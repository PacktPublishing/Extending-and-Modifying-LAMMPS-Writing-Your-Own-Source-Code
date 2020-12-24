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
#include "fix_addboost.h"
#include "atom.h"
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

//SM
#include "group.h"	
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NONE,CONSTANT,EQUAL,ATOM};

/* ---------------------------------------------------------------------- */

FixAddBoost::FixAddBoost(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 10) error->all(FLERR,"Illegal fix addboost command");	// SM

  dynamic_group_allow = 1;
  scalar_flag = 1;
  vector_flag = 1;
  size_vector = 3;
  global_freq = 1;
  extscalar = 1;
  extvector = 1;

  // ------------------------
  // SM: Parse new parameters
  // ------------------------
  
  jgroupstr = atomstr = NULL;
  
  // SM: Set second group for interaction
  int n = strlen(arg[3]) + 1;
  jgroupstr = new char[n];
  strcpy(jgroupstr,arg[3]);
  jgroup = group->find(jgroupstr);
  if (jgroup == -1)
    error->all(FLERR,"FixAddboost second group ID does not exist");
  jgroupbit = group->bitmask[jgroup];
   
  // SM: Read atom ID 
  if (strstr(arg[4],"v_") == arg[4]) {
    int n = strlen(&arg[4][2]) + 1;
    atomstr = new char[n];
    strcpy(atomstr,&arg[4][2]);
  } else {
    atomID = force->inumeric(FLERR,arg[4]);
    atomstyle = CONSTANT;
  }  
  
  // SM: Set type of boost potential (LJ, Morse, harmonic)
  if (strstr(arg[5],"LJ") == arg[5]) potentialType = 1;
  else if (strstr(arg[5],"morse") == arg[5]) potentialType = 2;
  else if (strstr(arg[5],"harmonic") == arg[5]) potentialType = 3;
  else error->all(FLERR,"Invalid potential type for fix addboost");

  // SM: Read parameters and cutoff
  param1 = force->numeric(FLERR,arg[6]);
  param2 = force->numeric(FLERR,arg[7]);
  param3 = force->numeric(FLERR,arg[8]);
  cutoff = force->numeric(FLERR,arg[9]); 

  force_flag = 0;
  foriginal[0] = foriginal[1] = foriginal[2] = foriginal[3] = 0.0;
}

/* ---------------------------------------------------------------------- */

FixAddBoost::~FixAddBoost()
{
  delete [] jgroupstr;	// SM
  delete [] atomstr;	// SM
}

/* ---------------------------------------------------------------------- */

int FixAddBoost::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= THERMO_ENERGY;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAddBoost::init()
{
    // SM: deleted region check
   
	// SM: check variable for atomID
	if (atomstr) {
	  atomvar = input->variable->find(atomstr);
	  if (atomvar < 0)
		  error->all(FLERR,"Atom ID for fix addboost does not exist");
	  if (input->variable->equalstyle(atomvar)) 
		  atomstyle = EQUAL;
	  else error->all(FLERR,"Atom ID for fix addboost is invalid style");
	}
	
	// SM: need a full neighbor list, built whenever re-neighboring occurs 
	// (see fix_orient_fcc.cpp)
	int irequest = neighbor->request((void *) this);
	neighbor->requests[irequest]->pair = 0;	// not called by pair
	neighbor->requests[irequest]->fix = 1;	// called by fix
	neighbor->requests[irequest]->half = 0;	// half-list not requested
	neighbor->requests[irequest]->full = 1;	// full-list requested
}

/* ---------------------------------------------------------------------- */

// SM: init_list() method to get hold of the neighbour list pointer
void FixAddBoost::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixAddBoost::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet"))
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(nlevels_respa-1);
    post_force_respa(vflag,nlevels_respa-1,0);
    ((Respa *) update->integrate)->copy_f_flevel(nlevels_respa-1);
  }
}

/* ---------------------------------------------------------------------- */

void FixAddBoost::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixAddBoost::post_force(int vflag)
{
  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  imageint *image = atom->image;
  int nlocal = atom->nlocal;

  modify->clearstep_compute(); 
	force_flag = 0;
  	foriginal[0] = foriginal[1] = foriginal[2] = foriginal[3] = 0.0;
	boostV = 0.0;	// SM

	// SM: declare variables
	int i,maxZID;
	int j,jj,jnum;
	double xtmp,ytmp,ztmp,delx,dely,delz;
	double rsq,fpair;
	int *jlist,*numneigh,**firstneigh;
	
	// SM: read atom ID to apply boost
	if (atomstyle == CONSTANT) maxZID = atomID;
 	else if (atomstyle == EQUAL) 
		maxZID = (int) input->variable->compute_equal(atomvar);
	i = atom->map(maxZID);

	// SM: check if atom belongs to core and is not a ghost atom
	if ((i >= 0) && (i < nlocal)) {
		xtmp = x[i][0];
		ytmp = x[i][1];
		ztmp = x[i][2];
		
		// SM: read neighbors of boosted atom
		numneigh = list->numneigh;
		firstneigh = list->firstneigh;
		jlist = firstneigh[i];
		jnum = numneigh[i];
			
		// SM: loop over all neighbors of boosted atom
		for (jj = 0; jj < jnum; jj++) {	
		  j = jlist[jj];
		  j &= NEIGHMASK;

		  if (mask[j] & jgroupbit) {
			  delx = xtmp - x[j][0];
			  dely = ytmp - x[j][1];
			  delz = ztmp - x[j][2];
			  rsq = delx*delx + dely*dely + delz*delz;
			  
			  if (rsq <= cutoff*cutoff) {	
				// SM: apply boost pair potential		  
				if (potentialType == 1) {
					// LJ: param1 = epsilon; param2 = sigma
					double r2inv = 1.0/rsq;
					double r6inv = r2inv*r2inv*r2inv;
					double sigma6 = param2*param2*param2*param2*param2*param2;
					double sigma12 = sigma6*sigma6;
					boostV += 4.0 * param1 * r6inv * (sigma12*r6inv - sigma6);		
					double forcelj = 24.0 * param1 * r6inv * (2.0*sigma12*r6inv - sigma6);	
					fpair = forcelj*r2inv;
				}
				else if (potentialType == 2) {
					// MORSE: param1 = D0; param2 = alpha; param3 = R0
					double r = sqrt(rsq);
					double dr = r - param3;
					double dexp = exp(-param2 * dr);
					boostV += param1 * (dexp*dexp - 2.0*dexp);
					fpair = 2.0 * param1 * param2 * (dexp*dexp - dexp) / r;
				}
				else if (potentialType == 3) {
					// HARMONIC: param1 = k; param2 = well_depth
					boostV += 0.5*param1*rsq - param2;
					fpair = -1.0*param1;
				}

				fpair *= -1.0;	// SM: repulsive force

				f[i][0] += delx*fpair;
				f[i][1] += dely*fpair;
				f[i][2] += delz*fpair;
			
				f[j][0] -= delx*fpair;
				f[j][1] -= dely*fpair;
				f[j][2] -= delz*fpair;
				
				foriginal[1] += f[i][0];
				foriginal[2] += f[i][1];
				foriginal[3] += f[i][2];
			  }
		  }
	    }  
	}
 	foriginal[0] = -1.0*boostV;		// SM: repulsive potential
    modify->addstep_compute(update->ntimestep + 1);			// SM
    modify->addstep_compute_all(update->ntimestep + 1);		// SM	 
}

/* ---------------------------------------------------------------------- */

void FixAddBoost::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixAddBoost::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   potential energy of added force
------------------------------------------------------------------------- */

double FixAddBoost::compute_scalar()
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

double FixAddBoost::compute_vector(int n)
{
  // only sum across procs one time

  if (force_flag == 0) {
    MPI_Allreduce(foriginal,foriginal_all,4,MPI_DOUBLE,MPI_SUM,world);
    force_flag = 1;
  }
  return foriginal_all[n+1];
}
