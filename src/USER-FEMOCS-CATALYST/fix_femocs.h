#ifdef FIX_CLASS

FixStyle(femocs,FixFemocs)

#else

#ifndef LMP_FIX_FEMOCS_H
#define LMP_FIX_FEMOCS_H

#include "fix.h"
#include "Femocs.h"

namespace LAMMPS_NS {

/**
  *  @class FixFemocs
  *  @brief Class for including effects of high electric field.
  */

class FixFemocs : public Fix {
public:
  FixFemocs(class LAMMPS *, int, char **);
  virtual ~FixFemocs();

  int setmask();
  void init();
  void setup(int);
  void min_setup(int);
  void post_force(int);
  void post_force_respa(int, int, int);
  void min_post_force(int);
  void end_of_step();

  double compute_scalar();
  double compute_vector(int);
  double memory_usage();

protected:
  femocs::Femocs femocs;

  double foriginal[3],foriginal_all[3];  ///< total force on atoms before extra force added
  int force_flag;
  int ilevel_respa;
  int maxatom;         ///< max # of owned+ghost in arrays on this proc
  double **sforce;     ///< extra force from FEMOCS

  double pot_energy;   ///< potential energy added/removed due to force scaling
  double kin_energy;   ///< kinetic energy added/removed due to temperature scaling
  int debug;           ///< print debug messages to the console

  void print_msg(char* msg);
};

}

#endif // LMP_FIX_FEMOCS_H
#endif // FIX_CLASS
