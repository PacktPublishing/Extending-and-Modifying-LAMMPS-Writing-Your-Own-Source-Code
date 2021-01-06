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

#include "fix_activeforce.h"
#include <mpi.h>
#include <cstring>
#include <cstdlib>
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

// SM (22May20): import header files
#include "atom_vec.h"
#include "molecule.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NONE,CONSTANT,EQUAL,ATOM};

/* ---------------------------------------------------------------------- */

FixActiveForce::FixActiveForce(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  xstr(NULL), ystr(NULL), zstr(NULL), estr(NULL), idregion(NULL), sforce(NULL)

{
  if (narg < 5) error->all(FLERR,"Illegal fix activeforce command");

  dynamic_group_allow = 1;
  scalar_flag = 1;
  vector_flag = 1;
  size_vector = 3;
  global_freq = 1;
  extscalar = 1;
  extvector = 1;
  respa_level_support = 1;
  ilevel_respa = 0;
  virial_flag = 1;

  xstr = ystr = zstr = NULL;

  if (strstr(arg[3],"v_") == arg[3]) {
    int n = strlen(&arg[3][2]) + 1;
    xstr = new char[n];
    strcpy(xstr,&arg[3][2]);
  } else {
    fMagnitude = fabs( force->numeric(FLERR,arg[3]) );	// SM (21May20): changed xvalue to fMagnitude
    xstyle = CONSTANT;
  }

	// SM (21May20): Removed fy and fz entries

	// -----------------------------------------------------------------
	// SM (21May20): choose forward or reverse direction of active force
	// -----------------------------------------------------------------
	if (strcmp(arg[4],"forward") == 0) {
		direction = 1;
	} 
	else if (strcmp(arg[4],"reverse") == 0) {
		direction = -1;
	}
	else error->all(FLERR,"Illegal fix activeforce direction: must be forward or reverse");
	// -----------------------------------------------------------------


  // optional args

  nevery = 1;
  iregion = -1;

  int iarg = 5;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"every") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix activeforce command");
      nevery = atoi(arg[iarg+1]);
      if (nevery <= 0) error->all(FLERR,"Illegal fix activeforce command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"region") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix activeforce command");
      iregion = domain->find_region(arg[iarg+1]);
      if (iregion == -1)
        error->all(FLERR,"Region ID for fix activeforce does not exist");
      int n = strlen(arg[iarg+1]) + 1;
      idregion = new char[n];
      strcpy(idregion,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"energy") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix activeforce command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) {
        int n = strlen(&arg[iarg+1][2]) + 1;
        estr = new char[n];
        strcpy(estr,&arg[iarg+1][2]);
      } else error->all(FLERR,"Illegal fix activeforce command");
      iarg += 2;
    } else error->all(FLERR,"Illegal fix activeforce command");
  }

  force_flag = 0;
  foriginal[0] = foriginal[1] = foriginal[2] = foriginal[3] = 0.0;

  maxatom = 1;
  memory->create(sforce,maxatom,4,"activeforce:sforce");
}

/* ---------------------------------------------------------------------- */

FixActiveForce::~FixActiveForce()
{
  delete [] xstr;
  delete [] ystr;
  delete [] zstr;
  delete [] estr;
  delete [] idregion;
  memory->destroy(sforce);
}

/* ---------------------------------------------------------------------- */

int FixActiveForce::setmask()
{
  datamask_read = datamask_modify = 0;

  int mask = 0;
  mask |= POST_FORCE;
  mask |= THERMO_ENERGY;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixActiveForce::init()
{
  // check variables

  if (xstr) {
    xvar = input->variable->find(xstr);
    if (xvar < 0)
      error->all(FLERR,"Variable name for fix activeforce does not exist");
    if (input->variable->equalstyle(xvar)) xstyle = EQUAL;
    else if (input->variable->atomstyle(xvar)) xstyle = ATOM;
    else error->all(FLERR,"Variable for fix activeforce is invalid style");
  }
  if (ystr) {
    yvar = input->variable->find(ystr);
    if (yvar < 0)
      error->all(FLERR,"Variable name for fix activeforce does not exist");
    if (input->variable->equalstyle(yvar)) ystyle = EQUAL;
    else if (input->variable->atomstyle(yvar)) ystyle = ATOM;
    else error->all(FLERR,"Variable for fix activeforce is invalid style");
  }
  if (zstr) {
    zvar = input->variable->find(zstr);
    if (zvar < 0)
      error->all(FLERR,"Variable name for fix activeforce does not exist");
    if (input->variable->equalstyle(zvar)) zstyle = EQUAL;
    else if (input->variable->atomstyle(zvar)) zstyle = ATOM;
    else error->all(FLERR,"Variable for fix activeforce is invalid style");
  }
  if (estr) {
    evar = input->variable->find(estr);
    if (evar < 0)
      error->all(FLERR,"Variable name for fix activeforce does not exist");
    if (input->variable->atomstyle(evar)) estyle = ATOM;
    else error->all(FLERR,"Variable for fix activeforce is invalid style");
  } else estyle = NONE;

  // set index and check validity of region

  if (iregion >= 0) {
    iregion = domain->find_region(idregion);
    if (iregion == -1)
      error->all(FLERR,"Region ID for fix activeforce does not exist");
  }

  if (xstyle == ATOM || ystyle == ATOM || zstyle == ATOM)
    varflag = ATOM;
  else if (xstyle == EQUAL || ystyle == EQUAL || zstyle == EQUAL)
    varflag = EQUAL;
  else varflag = CONSTANT;

  if (varflag == CONSTANT && estyle != NONE)
    error->all(FLERR,"Cannot use variable energy with "
               "constant force in fix activeforce");
  if ((varflag == EQUAL || varflag == ATOM) &&
      update->whichflag == 2 && estyle == NONE)
    error->all(FLERR,"Must use variable energy with fix activeforce");

  if (strstr(update->integrate_style,"respa")) {
    ilevel_respa = ((Respa *) update->integrate)->nlevels-1;
    if (respa_level >= 0) ilevel_respa = MIN(respa_level,ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixActiveForce::setup(int vflag)
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

void FixActiveForce::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixActiveForce::post_force(int vflag)
{
  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  imageint *image = atom->image;
  double v[6];
  int nlocal = atom->nlocal;

  if (update->ntimestep % nevery) return;

  // energy and virial setup

  if (vflag) v_setup(vflag);
  else evflag = 0;

  if (lmp->kokkos)
    atom->sync_modify(Host, (unsigned int) (F_MASK | MASK_MASK),
                      (unsigned int) F_MASK);

  // update region if necessary

  Region *region = NULL;
  if (iregion >= 0) {
    region = domain->regions[iregion];
    region->prematch();
  }

  // reallocate sforce array if necessary

  if ((varflag == ATOM || estyle == ATOM) && atom->nmax > maxatom) {
    maxatom = atom->nmax;
    memory->destroy(sforce);
    memory->create(sforce,maxatom,4,"activeforce:sforce");
  }

  // foriginal[0] = "potential energy" for added force
  // foriginal[123] = force on atoms before extra force added

  foriginal[0] = foriginal[1] = foriginal[2] = foriginal[3] = 0.0;
  force_flag = 0;


	// ---------------------------------------------------------
	// SM (21May20): introduce atom pointer to read bonded atoms
	// -------------------------------------------------------------
	  tagint **bond_atom = atom->bond_atom;
	  tagint *tag = atom->tag;
	  tagint tagprev;
	  int molecular = atom->molecular;
	  int *num_bond = atom->num_bond;
	  int *molindex = atom->molindex;
	  int *molatom = atom->molatom;
	  Molecule **onemols = atom->avec->onemols;
	  int atom1, atom2;
	  int jmol,jatom,nb;           
	  double rsq,delx,dely,delz;
	  int atomForce;
	// -------------------------------------------------------------


  // constant force
  // potential energy = - x dot f in unwrapped coords

  if (varflag == CONSTANT) {
    double unwrap[3];
    for (int i = 0; i < nlocal; i++) {		
		if (mask[i] & groupbit) {
			if (region && !region->match(x[i][0],x[i][1],x[i][2])) continue;
		
			// -------------------------------------------------------------
			// SM (21May20): Identifying bonds and direction of active force 
			// (from compute_bond_local.cpp)
			// -------------------------------------------------------------
			atom1 = i;
			
			if (molecular == 1) {
				nb = num_bond[atom1];
				//fprintf(screen,"atom1 = %d \n", atom1);
			}
			else {
				if (molindex[atom1] < 0) continue;
				jmol = molindex[atom1];
				jatom = molatom[atom1];
				nb = onemols[jmol]->num_bond[jatom];
			}
			
			for (int j = 0; j < nb; j++) {
				if (molecular == 1) {
					atom2 = atom->map(bond_atom[atom1][j]);
					//fprintf(screen,"      atom2 = %d \n", atom2);
			    } else {
					tagprev = tag[atom1] - jatom - 1;
					atom2 = atom->map(onemols[jmol]->bond_atom[jatom][j]+tagprev);
				}

				if (atom2 < 0 || !(mask[atom2] & groupbit)) continue;		
				
				//fprintf(screen,"				tag[atom1] = %d \n", tag[atom1]);
				//fprintf(screen,"				tag[atom2] = %d \n", tag[atom2]);
			
				delx = x[atom2][0]-x[atom1][0];
				dely = x[atom2][1]-x[atom1][1];
				delz = x[atom2][2]-x[atom1][2];
				domain->minimum_image(delx,dely,delz);
				rsq = delx*delx + dely*dely + delz*delz;
				
				xvalue = fMagnitude * delx/sqrt(rsq) * direction;
				yvalue = fMagnitude * dely/sqrt(rsq) * direction;
				zvalue = fMagnitude * delz/sqrt(rsq) * direction;
				
				if (direction == 1) atomForce = atom1;
				if (direction == -1) atomForce = atom2;
			// -------------------------------------------------------------
			
				domain->unmap(x[atomForce],image[atomForce],unwrap);				
				foriginal[0] -= xvalue*unwrap[0] + yvalue*unwrap[1] + zvalue*unwrap[2];
				foriginal[1] += f[atomForce][0];
				foriginal[2] += f[atomForce][1];
				foriginal[3] += f[atomForce][2];
				f[atomForce][0] += xvalue;
				f[atomForce][1] += yvalue;
				f[atomForce][2] += zvalue;
				if (evflag) {
				  v[0] = xvalue * unwrap[0];
				  v[1] = yvalue * unwrap[1];
				  v[2] = zvalue * unwrap[2];
				  v[3] = xvalue * unwrap[1];
				  v[4] = xvalue * unwrap[2];
				  v[5] = yvalue * unwrap[2];
				  v_tally(atomForce,v);
				}
			}
		}	
	}
	
  // variable force, wrap with clear/add
  // potential energy = evar if defined, else 0.0
  // wrap with clear/add

  } else {
    double unwrap[3];

    modify->clearstep_compute();

    if (xstyle == EQUAL) fMagnitude = input->variable->compute_equal(xvar);
    else if (xstyle == ATOM)
      input->variable->compute_atom(xvar,igroup,&sforce[0][0],4,0);
    if (estyle == ATOM)
      input->variable->compute_atom(evar,igroup,&sforce[0][3],4,0);

    modify->addstep_compute(update->ntimestep + 1);

    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        if (region && !region->match(x[i][0],x[i][1],x[i][2])) continue;
		
		// -------------------------------------------------------------
		// SM (21May20): Identifying bonds and direction of active force 
		// (from compute_bond_local.cpp)
		// -------------------------------------------------------------
		atom1 = i;
		
		if (molecular == 1) {
			nb = num_bond[atom1];
		}
		else {
			if (molindex[atom1] < 0) continue;
			jmol = molindex[atom1];
			jatom = molatom[atom1];
			nb = onemols[jmol]->num_bond[jatom];
		}
		
		for (int j = 0; j < nb; j++) {
			if (molecular == 1) {
				atom2 = atom->map(bond_atom[atom1][j]);
			} else {
				tagprev = tag[atom1] - jatom - 1;
				atom2 = atom->map(onemols[jmol]->bond_atom[jatom][j]+tagprev);
			}

			if (atom2 < 0 || !(mask[atom2] & groupbit)) continue;		
		
			delx = x[atom2][0]-x[atom1][0];
			dely = x[atom2][1]-x[atom1][1];
			delz = x[atom2][2]-x[atom1][2];
			domain->minimum_image(delx,dely,delz);
			rsq = delx*delx + dely*dely + delz*delz;
			
			xvalue = fMagnitude * delx/sqrt(rsq) * direction;
			yvalue = fMagnitude * dely/sqrt(rsq) * direction;
			zvalue = fMagnitude * delz/sqrt(rsq) * direction;
			
			if (direction == 1) atomForce = atom1;
			if (direction == -1) atomForce = atom2;
		// -------------------------------------------------------------
		
			domain->unmap(x[atomForce],image[atomForce],unwrap);
			
			// SM (21May20): removed ystyle and zstyle
			if (xstyle == ATOM) {
				xvalue = sforce[atomForce][0];
				yvalue = sforce[atomForce][1];
				zvalue = sforce[atomForce][2];
			}

			// SM (21May20): removed ystyle and zstyle
			if (estyle == ATOM) {
			  foriginal[0] += sforce[atomForce][3];
			} else {
			  if (xstyle) {
				  foriginal[0] -= xvalue*unwrap[0];
				  foriginal[0] -= yvalue*unwrap[1];
				  foriginal[0] -= zvalue*unwrap[2];
			  }
			}
			
			foriginal[1] += f[atomForce][0];
			foriginal[2] += f[atomForce][1];
			foriginal[3] += f[atomForce][2];

			// SM (21May20): removed ystyle and zstyle
			if (xstyle) {
				f[atomForce][0] += xvalue;
				f[atomForce][1] += yvalue;
				f[atomForce][2] += zvalue;
			}
			
			if (evflag) {
			  // SM (21May20): changed ystyle (v[1],v[5]) and zstyle (v[2]) to xstyle
			  v[0] = xstyle ? xvalue*unwrap[0] : 0.0;
			  v[1] = xstyle ? yvalue*unwrap[1] : 0.0;
			  v[2] = xstyle ? zvalue*unwrap[2] : 0.0;
			  v[3] = xstyle ? xvalue*unwrap[1] : 0.0;
			  v[4] = xstyle ? xvalue*unwrap[2] : 0.0;
			  v[5] = xstyle ? yvalue*unwrap[2] : 0.0;
			  v_tally(atomForce, v);
			}
		}		
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixActiveForce::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == ilevel_respa) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixActiveForce::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   potential energy of added force
------------------------------------------------------------------------- */

double FixActiveForce::compute_scalar()
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

double FixActiveForce::compute_vector(int n)
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

double FixActiveForce::memory_usage()
{
  double bytes = 0.0;
  if (varflag == ATOM) bytes = maxatom*4 * sizeof(double);
  return bytes;
}
