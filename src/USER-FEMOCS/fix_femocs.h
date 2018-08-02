#ifdef FIX_CLASS

FixStyle(femocs,FixFemocs)

#else

#ifndef LMP_FIX_FEMOCS_H
#define LMP_FIX_FEMOCS_H

#include "fix.h"
#include "Femocs.h"
#include <stdio.h>
#include <stdlib.h>

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

  // TODO clarify & clean them up

  double xvalue,yvalue,zvalue;
  int varflag,iregion;
  char *xstr,*ystr,*zstr,*estr;
  char *idregion;
  int xvar,yvar,zvar,evar,xstyle,ystyle,zstyle,estyle;
  double foriginal[4],foriginal_all[4];
  int force_flag;
  int ilevel_respa;
  int maxatom;
  double **sforce;

  void print(char* msg);
};

}

#endif // LMP_FIX_FEMOCS_H
#endif // FIX_CLASS
