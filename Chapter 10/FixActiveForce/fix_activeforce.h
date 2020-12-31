/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(activeforce,FixActiveForce)

#else

#ifndef LMP_FIX_ACTIVEFORCE_H
#define LMP_FIX_ACTIVEFORCE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixActiveForce : public Fix {
 public:
  FixActiveForce(class LAMMPS *, int, char **);
  ~FixActiveForce();
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

 private:
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
  
  // SM (21May20): integer to specify direction of active force
  int direction;	
  
  // SM (21May20): variable to store force magnitude 
  double fMagnitude;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Region ID for fix activeforce does not exist

Self-explanatory.

E: Variable name for fix activeforce does not exist

Self-explanatory.

E: Variable for fix activeforce is invalid style

Self-explanatory.

E: Cannot use variable energy with constant force in fix activeforce

This is because for constant force, LAMMPS can compute the change
in energy directly.

E: Must use variable energy with fix activeforce

Must define an energy variable when applying a dynamic
force during minimization.

*/
