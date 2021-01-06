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

#include "string.h"
#include "stdlib.h"
#include "math.h"
#include "fix_temp_andersen.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "input.h"
#include "variable.h"
#include "group.h"
#include "update.h"
#include "modify.h"
#include "compute.h"
#include "error.h"
#include "random_mars.h"	// SM 

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NOBIAS,BIAS};
enum{CONSTANT,EQUAL};

/* ---------------------------------------------------------------------- */

FixTempAndersen::FixTempAndersen(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 6) error->all(FLERR,"Illegal fix temp/andersen command");	// SM

  nevery = 1;

// SM 
  andFreq = force->inumeric(FLERR,arg[3]);
  seedMB = force->inumeric(FLERR,arg[4]);
  sigmaMB = force->numeric(FLERR,arg[5]);

// SM: initialize Marsaglia RNG with processor-unique seed
  random = new RanMars(lmp,seedMB + comm->me);	
}

/* ---------------------------------------------------------------------- */

int FixTempAndersen::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  mask |= THERMO_ENERGY;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixTempAndersen::end_of_step()
{

  double **v = atom->v;
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

// SM: 
// 1. adjust velocities if currentTime % andFreq == 0
// 2. generate vx,vy,vz from random->gaussian() and scale by sigmaMB

  int currentTime = update->ntimestep;	// SM 
  if ( (currentTime > 0) && (currentTime % andFreq == 0) ) {
    for (int i = 0; i < nlocal; i++) {	
      if (mask[i] & groupbit) {
		// SM: assign vx,vy,vz
		v[i][0] = sigmaMB * random->gaussian();
		v[i][1] = sigmaMB * random->gaussian();
		v[i][2] = sigmaMB * random->gaussian();
      }       
    }	
  } 
}

/* ---------------------------------------------------------------------- */

