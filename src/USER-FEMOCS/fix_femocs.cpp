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

#include <cstring>
#include <cstdlib>
#include "fix_femocs.h"
#include "comm.h"
#include "math_extra.h"
#include "atom.h"
#include "atom_masks.h"
#include "update.h"
#include "modify.h"
#include "domain.h"
#include "region.h"
#include "respa.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"
#include "force.h"


using namespace LAMMPS_NS;

enum{NONE,CONSTANT,EQUAL,ATOM};

/* ---------------------------------------------------------------------- */

FixFemocs::FixFemocs(LAMMPS *lmp, int narg, char **arg) :
        Fix(lmp, narg, arg)
{
    if (narg != 4)
        error->all(FLERR,"Illegal fix femocs command");

    femocs.read_conf(arg[3]);
}

/* ---------------------------------------------------------------------- */

FixFemocs::~FixFemocs()
{}

/* ---------------------------------------------------------------------- */

void FixFemocs::print(char* msg)
{
    if (comm->me == 0) printf("FixFemocs: %s\n", msg);
}

/* ---------------------------------------------------------------------- */

int FixFemocs::setmask()
{
  datamask_read = datamask_modify = 0;

  int mask = 0;
  mask |= FixConst::POST_FORCE;
  mask |= FixConst::THERMO_ENERGY;
  mask |= FixConst::POST_FORCE_RESPA;
  mask |= FixConst::MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixFemocs::post_force(int vflag)
{
    require(comm->nprocs == 1, "fix femocs can run only on single CPU core!");

    double *xyz_1d = &atom->x[0][0];
    double *velocity_1d = &atom->v[0][0];
    double *force_1d = &atom->f[0][0];

    atom->nmax;

    int n_atoms = atom->nlocal;

    // import atomistic data to Femocs
    // NB! Atom indices change between time steps due to sorting
    //     Maybe storing atoms helps as described in Developer.pdf
    if (femocs.import_atoms(n_atoms, xyz_1d)) {
        print("importing atoms failed!");
        return;
    }

    // solve equations
    if (femocs.run()) {
        print("solving equations failed!");
        return;
    }

    // export forces
    if (femocs.export_data(force_1d, n_atoms, "force")) {
        print("exporting forces failed!");
        return;
    }

    // TODO handle energies

//    // export total potential energy of atoms without contribution of Lorentz interaction
//    success += femocs.export_data(Epair_all, 1, "pair_potential_sum");
//    Vpair = Epair_all[0];
//    // export potential energy per atom (pair-potential)
//    success += femocs.export_data(Epair_all, n_atoms, "pair_potential");
    //    // scale velocities
//    success += femocs.export_data(v, n_atoms, "velocity");
}

/* ---------------------------------------------------------------------- */

void FixFemocs::init()
{
  // TODO check what's going on here

}

/* ---------------------------------------------------------------------- */

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
   potential energy of added force
------------------------------------------------------------------------- */

double FixFemocs::compute_scalar()
{
  // only sum across procs one time

  if (force_flag == 0) {
    MPI_Allreduce(foriginal,foriginal_all,4,MPI_DOUBLE,MPI_SUM,world);
    force_flag = 1;
  }
  return foriginal_all[0];
}

/* ----------------------------------------------------------------------
   return components of total force on fix group before force was changed
------------------------------------------------------------------------- */

double FixFemocs::compute_vector(int n)
{
  // only sum across procs one time

  if (force_flag == 0) {
    MPI_Allreduce(foriginal,foriginal_all,4,MPI_DOUBLE,MPI_SUM,world);
    force_flag = 1;
  }
  return foriginal_all[n+1];
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixFemocs::memory_usage()
{
  double bytes = 0.0;
  if (varflag == ATOM) bytes = maxatom*4 * sizeof(double);
  return bytes;
}
