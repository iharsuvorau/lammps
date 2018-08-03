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

  double compute_scalar();
  double compute_vector(int);
  double memory_usage();

protected:
  femocs::Femocs femocs;

  double foriginal[4],foriginal_all[4];
  int force_flag;
  int ilevel_respa;
  int maxatom;
  double **sforce;
  double *pair_pot;

  void print_msg(char* msg);
};

}

#endif // LMP_FIX_FEMOCS_H
#endif // FIX_CLASS
