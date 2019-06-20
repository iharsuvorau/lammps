/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "fix_femocs.h"
#include "comm.h"
#include "math_extra.h"
#include "atom.h"
#include "atom_masks.h"
#include "update.h"
#include "domain.h"
#include "respa.h"
#include "memory.h"

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
 * Read the arguments and initialize internals
 * ---------------------------------------------------------------------- */

FixFemocs::FixFemocs(LAMMPS *lmp, int narg, char **arg) :
                Fix(lmp, narg, arg),
                force_flag(0), ilevel_respa(0), maxatom(1), sforce(NULL),
                kin_energy(0), pot_energy(0), debug(0)
{
  // check for number of processors
  if (comm->nprocs > 1)
    error->all(FLERR,"Fix femocs can run only on single CPU core");

  // check for the existence of path to Femocs conf file
  if (narg < 4)
    error->all(FLERR,"Illegal fix femocs command");

  // read debug output flag
  if (narg > 4)
    debug = atoi(arg[4]);

  // read Femocs configuration parameters
  femocs.init(arg[3]);


  // flags related to force & velocity modification
  dynamic_group_allow = 1;
  scalar_flag = 1;
  vector_flag = 1;
  size_vector = 3;
  extscalar = 1;
  extvector = 1;
  respa_level_support = 1;
  virial_flag = 1;
  nevery = 1;
  global_freq = nevery;

  memory->create(sforce,maxatom,3,"femocs:sforce");
}

/* ---------------------------------------------------------------------- */

FixFemocs::~FixFemocs()
{
  memory->destroy(sforce);
}

/* ---------------------------------------------------------------------- */


void FixFemocs::print_msg(char* msg)
{
  if (debug > 0 && comm->me == 0) printf("FixFemocs: %s\n", msg);
}

/* ----------------------------------------------------------------------
 * Specify the routines LAMMPS engine is going to call
 * ---------------------------------------------------------------------- */

int FixFemocs::setmask()
{
  datamask_read = datamask_modify = 0;

  int mask = 0;
  mask |= FixConst::POST_FORCE;
  mask |= FixConst::THERMO_ENERGY;
  mask |= FixConst::POST_FORCE_RESPA;
  mask |= FixConst::MIN_POST_FORCE;
  mask |= FixConst::END_OF_STEP;
  return mask;
}

/* ----------------------------------------------------------------------
 * Continue initializing internals
 * ---------------------------------------------------------------------- */

void FixFemocs::init()
{
  if (strstr(update->integrate_style,"respa")) {
    ilevel_respa = ((Respa *) update->integrate)->nlevels-1;
    if (respa_level >= 0)
      ilevel_respa = MIN(respa_level,ilevel_respa);
  }
}

/* ----------------------------------------------------------------------
 * Main method to calculate forces.
 * ---------------------------------------------------------------------- */

void FixFemocs::post_force(int vflag)
{
  if (comm->nprocs > 1)
    error->all(FLERR,"Fix femocs can run only on single CPU core");

  // energy and virial setup
  if (vflag) v_setup(vflag);
  else evflag = 0;

  double virial[6];
  double unwrap[3];

  double **x = atom->x;
  double **f = atom->f;
//  int *tag = atom->tag;  // might be needed to perform properly RMSD check
  int *mask = atom->mask;
  imageint *image = atom->image;
  int nlocal = atom->nlocal;

  foriginal[0] = foriginal[1] = foriginal[2] = 0.0;
  force_flag = 0;

  // store total force before modification
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      foriginal[0] += f[i][0];
      foriginal[1] += f[i][1];
      foriginal[2] += f[i][2];
    }
  }

  // NB! Due to sorting, atom indices change between time steps
  //     Maybe storing atoms helps as described in Developer.pdf
  print_msg("importing atoms...");
  if (femocs.import_lammps(nlocal, x, NULL, mask, groupbit)) {
    print_msg("importing atoms failed!");
    return;
  }

  print_msg("solving equations...");
  if (femocs.run(update->ntimestep)) {
    print_msg("solving equations failed, using previous solution!");
  }

  // reallocate sforce array if necessary, and clean it
  if (atom->nmax > maxatom) {
    maxatom = atom->nmax;
    memory->destroy(sforce);
    memory->create(sforce,maxatom,3,"femocs:sforce");
  }
  for (int i = 0; i < nlocal; i++) {
    sforce[i][0] = sforce[i][1] = sforce[i][2] = 0;
  }

  print_msg("exporting forces...");
  if (femocs.export_data(&sforce[0][0], nlocal, "force")) {
    print_msg("exporting forces failed!");
    return;
  }

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      // add electrostatic force
      f[i][0] += sforce[i][0];
      f[i][1] += sforce[i][1];
      f[i][2] += sforce[i][2];

      // if requested, store atom-wise added energy
      if (evflag) {
        domain->unmap(x[i],image[i],unwrap);

        virial[0] = sforce[i][0] * unwrap[0];
        virial[1] = sforce[i][1] * unwrap[1];
        virial[2] = sforce[i][2] * unwrap[2];
        virial[3] = sforce[i][0] * unwrap[1];
        virial[4] = sforce[i][0] * unwrap[2];
        virial[5] = sforce[i][1] * unwrap[2];
        v_tally(i, virial);
      }
    }
  }

  print_msg("exporting potential energy...");
  if (femocs.export_data(&pot_energy, nlocal, "pot_energy")) {
    print_msg("exporting potential energy failed!");
    return;
  }
}

void FixFemocs::post_force_manycore(int vflag)
{
  /*

  // check if it's time to rerun field solver
  if (update->ntimestep % nevery) return;

  // energy and virial setup
  if (vflag) v_setup(vflag);
  else evflag = 0;

  if (lmp->kokkos)
    atom->sync_modify(Host, (unsigned int) (F_MASK | MASK_MASK),
        (unsigned int) F_MASK);

  double virial[6];
  double unwrap[3];

  double **x = atom->x;
  double **f = atom->f;
  double **v = atom->v;
  int *mask = atom->mask;
  imageint *image = atom->image;

  int nlocal = atom->nlocal;
  int nglobal = atom->natoms;
  int n_procs = comm->nprocs;

  // foriginal[0]     = "potential energy" for added force
  // foriginal[1,2,3] = force on atoms before extra force added
  foriginal[0] = foriginal[1] = foriginal[2] = foriginal[3] = 0.0;
  force_flag = 0;

  // store total force before modification
  print_msg("storing forces...");
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      foriginal[1] += f[i][0];
      foriginal[2] += f[i][1];
      foriginal[3] += f[i][2];
    }
  }

  // TODO make proper allocation
  int atoms_per_core[32];
  int displ[32];
  double x_all[1000];
  double local_forces[1000];

  MPI_Gather(nlocal, 1, MPI_INT, atoms_per_core, 1, MPI_INT, 0, world);

  // TODO check proper displacements

  for (int i = 0; i < n_procs; ++i)
    for (int j = 0; j < atoms_per_core[i]; ++j)
      displ[i] += atoms_per_core[j];

  for (int i = 0; i < n_procs; ++i) {
    atoms_per_core[i] *= 3;
    displ[i] *= 3;
  }

  // TODO check how to properly gather 2D arrays

  MPI_Gatherv(x, 3*nlocal, MPI_DOUBLE, x_all, atoms_per_core, displ, MPI_DOUBLE, 0, world);

  // TODO handle errors properly

  if (comm->me == 0) {
    if (run_femocs(nglobal, &x[0][0])) return;
    if (export_forces(nglobal)) return;
  }

  // TODO distribute pair potential between cores

  foriginal[0] = 0;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit)
      foriginal[0] += pair_pot[i];
  }

  // TODO distribute forces between cores

  print_msg("modifying forces & energies...");
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      // add electrostatic force
      f[i][0] += sforce[i][0];
      f[i][1] += sforce[i][1];
      f[i][2] += sforce[i][2];

      // if requested, store atom-wise added energy
      if (evflag) {
        domain->unmap(x[i],image[i],unwrap);

        virial[0] = sforce[i][0] * unwrap[0];
        virial[1] = sforce[i][1] * unwrap[1];
        virial[2] = sforce[i][2] * unwrap[2];
        virial[3] = sforce[i][0] * unwrap[1];
        virial[4] = sforce[i][0] * unwrap[2];
        virial[5] = sforce[i][1] * unwrap[2];
        v_tally(i, virial);
      }
    }
  }

  //*/
}

void FixFemocs::end_of_step()
{
  // NB! Velocity scaling must take place every timestep
  double **v = atom->v;
  int nlocal = atom->nlocal;

  print_msg("scaling velocities...");
  // import velocities here as there has been integration meanwhile
  femocs.import_lammps(nlocal, NULL, v, atom->mask, groupbit);
  if (femocs.export_data(&v[0][0], nlocal, "velocity")) {
    print_msg("scaling velocities failed!");
    return;
  }

  print_msg("exporting kinetic energy...");
  if (femocs.export_data(&kin_energy, nlocal, "kin_energy")) {
    print_msg("exporting kinetic energy failed!");
    return;
  }
}

/* ----------------------------------------------------------------------
 * Calculate forces before the first time step.
 * It is required by the velocity-Verlet method.
 * ---------------------------------------------------------------------- */

void FixFemocs::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet"))
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(ilevel_respa);
    post_force_respa(vflag,ilevel_respa,0);
    ((Respa *) update->integrate)->copy_f_flevel(ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
 * Setup in case the runs are continued and
 * information from the previous run is still valid.
 * ---------------------------------------------------------------------- */

void FixFemocs::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixFemocs::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == ilevel_respa) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixFemocs::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   potential energy of added force and kinetic energy of scaled velocity
------------------------------------------------------------------------- */

double FixFemocs::compute_scalar()
{
  return pot_energy + kin_energy;
}

/* ----------------------------------------------------------------------
   return components of total force on fix group before force was changed
------------------------------------------------------------------------- */

double FixFemocs::compute_vector(int n)
{
  // sum across procs only once

  if (force_flag == 0) {
    MPI_Allreduce(foriginal,foriginal_all,3,MPI_DOUBLE,MPI_SUM,world);
    force_flag = 1;
  }
  return foriginal_all[n];
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array and Femocs object
------------------------------------------------------------------------- */

double FixFemocs::memory_usage()
{
  // TODO! Implement femocs.size()
  return maxatom * 3 * sizeof(double) + sizeof(femocs);
}

