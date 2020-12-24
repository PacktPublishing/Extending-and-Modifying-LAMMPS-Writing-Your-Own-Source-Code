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

FixStyle(temp/andersen,FixTempAndersen)

#else

#ifndef LMP_FIX_TEMP_ANDERSEN_H
#define LMP_FIX_TEMP_ANDERSEN_H

#include "fix.h"

namespace LAMMPS_NS {

class FixTempAndersen : public Fix {
 public:
  FixTempAndersen(class LAMMPS *, int, char **);
  int setmask();
  void end_of_step();

 
 private:
 
  double sigmaMB;			// SM
  int seedMB;				// SM
  int andFreq;				// SM
  class RanMars *random;	// SM
  int me;					// SM
  
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Fix temp/andersen period must be > 0.0

Self-explanatory.

E: Variable name for fix temp/andersen does not exist

Self-explanatory.

E: Variable for fix temp/andersen is invalid style

Only equal-style variables can be used.

E: Temperature ID for fix temp/andersen does not exist

Self-explanatory.

E: Computed temperature for fix temp/andersen cannot be 0.0

Self-explanatory.

E: Fix temp/andersen variable returned negative temperature

Self-explanatory.

E: Could not find fix_modify temperature ID

The compute ID for computing temperature does not exist.

E: Fix_modify temperature ID does not compute temperature

The compute ID assigned to the fix must compute temperature.

W: Group for fix_modify temp != fix group

The fix_modify command is specifying a temperature computation that
computes a temperature on a different group of atoms than the fix
itself operates on.  This is probably not what you want to do.

*/
