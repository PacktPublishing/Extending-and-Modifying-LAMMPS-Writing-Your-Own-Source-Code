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

#include "fix_addforceXY.h"
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

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NONE,CONSTANT,EQUAL,ATOM};

/* ---------------------------------------------------------------------- */

FixAddForceXY::FixAddForceXY(LAMMPS *lmp, int narg, char **arg) :	// SM
  Fix(lmp, narg, arg),
  xstr(NULL), ystr(NULL), zstr(NULL), estr(NULL), idregion(NULL), sforce(NULL),
  kxstr(NULL), kystr(NULL), ework(NULL)

{
  if (narg < 8) error->all(FLERR,"Illegal fix AddForceXY command");	// SM

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
  kxstr = kystr = NULL;		// SM

  if (strstr(arg[3],"v_") == arg[3]) {
    int n = strlen(&arg[3][2]) + 1;
    xstr = new char[n];
    strcpy(xstr,&arg[3][2]);
  } else {
    x0 = force->numeric(FLERR,arg[3]);	// SM
    xstyle = CONSTANT;
  }
  if (strstr(arg[4],"v_") == arg[4]) {
    int n = strlen(&arg[4][2]) + 1;
    ystr = new char[n];
    strcpy(ystr,&arg[4][2]);
  } else {
    y0 = force->numeric(FLERR,arg[4]);	// SM
    ystyle = CONSTANT;
  }
  if (strstr(arg[5],"v_") == arg[5]) {
    int n = strlen(&arg[5][2]) + 1;
    zstr = new char[n];
    strcpy(zstr,&arg[5][2]);
  } else {
    zvalue = force->numeric(FLERR,arg[5]);
    zstyle = CONSTANT;
  }
  if (strstr(arg[6],"v_") == arg[6]) {	// SM
    int n = strlen(&arg[6][2]) + 1;
    kxstr = new char[n];
    strcpy(kxstr,&arg[6][2]);
  } else {
    kx = force->numeric(FLERR,arg[6]);
    kxstyle = CONSTANT;
  }
  if (strstr(arg[7],"v_") == arg[7]) {	// SM
    int n = strlen(&arg[7][2]) + 1;
    kystr = new char[n];
    strcpy(kystr,&arg[7][2]);
  } else {
    ky = force->numeric(FLERR,arg[7]);
    kystyle = CONSTANT;
  }

  // optional args

  nevery = 1;
  iregion = -1;

  int iarg = 8;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"every") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix AddForceXY command");
      nevery = atoi(arg[iarg+1]);
      if (nevery <= 0) error->all(FLERR,"Illegal fix AddForceXY command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"region") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix AddForceXY command");
      iregion = domain->find_region(arg[iarg+1]);
      if (iregion == -1)
        error->all(FLERR,"Region ID for fix AddForceXY does not exist");
      int n = strlen(arg[iarg+1]) + 1;
      idregion = new char[n];
      strcpy(idregion,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"energy") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix AddForceXY command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) {
        int n = strlen(&arg[iarg+1][2]) + 1;
        estr = new char[n];
        strcpy(estr,&arg[iarg+1][2]);
      } else error->all(FLERR,"Illegal fix AddForceXY command");
      iarg += 2;
    } else error->all(FLERR,"Illegal fix AddForceXY command");
  }

  force_flag = 0;
  foriginal[0] = foriginal[1] = foriginal[2] = foriginal[3] = 0.0;

  maxatom = 1;
  memory->create(sforce,maxatom,4,"addforceXY:sforce");
  memory->create(ework,1,"addforceXY::ework");	// SM
}

/* ---------------------------------------------------------------------- */

FixAddForceXY::~FixAddForceXY()
{
  delete [] xstr;
  delete [] ystr;
  delete [] zstr;
  delete [] estr;
  delete [] idregion;
  memory->destroy(sforce);
  memory->destroy(ework);	// SM
}

/* ---------------------------------------------------------------------- */

int FixAddForceXY::setmask()
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

void FixAddForceXY::init()
{
  // check variables

  if (xstr) {
    xvar = input->variable->find(xstr);
    if (xvar < 0)
      error->all(FLERR,"Variable name for fix AddForceXY does not exist");
    if (input->variable->equalstyle(xvar)) xstyle = EQUAL;
    else error->all(FLERR,"Variable for fix AddForceXY is invalid style");
  }
  if (ystr) {
    yvar = input->variable->find(ystr);
    if (yvar < 0)
      error->all(FLERR,"Variable name for fix AddForceXY does not exist");
    if (input->variable->equalstyle(yvar)) ystyle = EQUAL;
    else error->all(FLERR,"Variable for fix AddForceXY is invalid style");
  }
  if (zstr) {
    zvar = input->variable->find(zstr);
    if (zvar < 0)
      error->all(FLERR,"Variable name for fix AddForceXY does not exist");
    if (input->variable->equalstyle(zvar)) zstyle = EQUAL;
    else if (input->variable->atomstyle(zvar)) zstyle = ATOM;
    else error->all(FLERR,"Variable for fix AddForceXY is invalid style");
  }
  if (estr) {
    evar = input->variable->find(estr);
    if (evar < 0)
      error->all(FLERR,"Variable name for fix AddForceXY does not exist");
    if (input->variable->atomstyle(evar)) estyle = ATOM;
    else error->all(FLERR,"Variable for fix AddForceXY is invalid style");
  } else estyle = NONE;
  
  // SM
  if (kxstr) {
    kxvar = input->variable->find(kxstr);
    if (kxvar < 0)
      error->all(FLERR,"Variable name for fix AddForceXY does not exist");
    if (input->variable->equalstyle(kxvar)) kxstyle = EQUAL;
	else error->all(FLERR,"Variable kxvar for fix AddForceXY is invalid style");
  } else kxstyle = NONE;

  if (kystr) {
    kyvar = input->variable->find(kystr);
    if (kyvar < 0)
      error->all(FLERR,"Variable name for fix AddForceXY does not exist");
    if (input->variable->equalstyle(kyvar)) kystyle = EQUAL;
    else error->all(FLERR,"Variable kyvar for fix AddForceXY is invalid style");
  } else kystyle = NONE;

  // set index and check validity of region

  if (iregion >= 0) {
    iregion = domain->find_region(idregion);
    if (iregion == -1)
      error->all(FLERR,"Region ID for fix AddForceXY does not exist");
  }

  if (zstyle == ATOM)	// SM
    varflag = ATOM;
  else if (xstyle == EQUAL || ystyle == EQUAL || zstyle == EQUAL || kxstyle == EQUAL || kystyle == EQUAL)	// SM
    varflag = EQUAL;
  else varflag = CONSTANT;

  if (varflag == CONSTANT && estyle != NONE)
    error->all(FLERR,"Cannot use variable energy with "
               "constant force in fix AddForceXY");
  if ((varflag == EQUAL || varflag == ATOM) &&
      update->whichflag == 2 && estyle == NONE)
    error->all(FLERR,"Must use variable energy with fix AddForceXY");

  if (strstr(update->integrate_style,"respa")) {
    ilevel_respa = ((Respa *) update->integrate)->nlevels-1;
    if (respa_level >= 0) ilevel_respa = MIN(respa_level,ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixAddForceXY::setup(int vflag)
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

void FixAddForceXY::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixAddForceXY::post_force(int vflag)
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
	memory->destroy(ework);	// SM
    memory->create(sforce,maxatom,4,"addforceXY:sforce");
	memory->create(ework,1,"addforceXY::ework");	// SM
  }

  // foriginal[0] = "potential energy" for added force
  // foriginal[123] = force on atoms before extra force added

  foriginal[0] = foriginal[1] = foriginal[2] = foriginal[3] = 0.0;
  force_flag = 0;

  // constant force
  // potential energy = - x dot f in unwrapped coords

  if (varflag == CONSTANT) {
    double unwrap[3];
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        if (region && !region->match(x[i][0],x[i][1],x[i][2])) continue;
        domain->unmap(x[i],image[i],unwrap);
        foriginal[0] -= -kx*(x[i][0] - x0)*unwrap[0] - ky*(x[i][1] - y0)*unwrap[1];	// SM
        foriginal[1] += f[i][0];
        foriginal[2] += f[i][1];
        foriginal[3] += f[i][2];
		
		// SM
        f[i][0] += -kx*(x[i][0] - x0);
        f[i][1] += -ky*(x[i][1] - y0);
        f[i][2] += 0.0;
		
		// SM
        if (evflag) {
          v[0] = -kx*(x[i][0] - x0) * unwrap[0];
          v[1] = -ky*(x[i][1] - y0) * unwrap[1];
          v[2] = 0.0 * unwrap[2];
          v[3] = -kx*(x[i][0] - x0) * unwrap[1];
          v[4] = -kx*(x[i][0] - x0) * unwrap[2];
          v[5] = -ky*(x[i][1] - y0) * unwrap[2];
          v_tally(i,v);
        }
      }

  // variable force, wrap with clear/add
  // potential energy = evar if defined, else 0.0
  // wrap with clear/add

  } else {
    double unwrap[3];

    modify->clearstep_compute();

    if (xstyle == EQUAL) x0 = input->variable->compute_equal(xvar);
    if (ystyle == EQUAL) y0 = input->variable->compute_equal(yvar);
    if (zstyle == EQUAL) zvalue = input->variable->compute_equal(zvar);
    if (estyle == ATOM)
      input->variable->compute_atom(evar,igroup,&ework[0],1,0);

    // SM
    if (kxstyle == EQUAL) kx = input->variable->compute_equal(kxvar);
    if (kystyle == EQUAL) ky = input->variable->compute_equal(kyvar);

    modify->addstep_compute(update->ntimestep + 1);

    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        if (region && !region->match(x[i][0],x[i][1],x[i][2])) continue;
        domain->unmap(x[i],image[i],unwrap);
		
        if (estyle == ATOM) {
          foriginal[0] += ework[0];	// SM
        } else {
		  // SM
          if (xstyle) foriginal[0] -= -kx*(x[i][0] - x0)*unwrap[0];
          if (ystyle) foriginal[0] -= -ky*(x[i][1] - y0)*unwrap[1];
          if (zstyle) foriginal[0] -= 0.0;
        }
        foriginal[1] += f[i][0];
        foriginal[2] += f[i][1];
        foriginal[3] += f[i][2];

		// SM
        if (xstyle) f[i][0] += -kx*(x[i][0] - x0);
        if (ystyle) f[i][1] += -ky*(x[i][1] - y0);
        if (zstyle) f[i][2] += 0.0;
        if (evflag) {
          v[0] = xstyle ? (-kx*(x[i][0] - x0))*unwrap[0] : 0.0;
          v[1] = ystyle ? (-ky*(x[i][1] - y0))*unwrap[1] : 0.0;
          v[2] = zstyle ? 0.0*unwrap[2] : 0.0;
          v[3] = xstyle ? (-kx*(x[i][0] - x0))*unwrap[1] : 0.0;
          v[4] = xstyle ? (-kx*(x[i][0] - x0))*unwrap[2] : 0.0;
          v[5] = ystyle ? (-ky*(x[i][1] - y0))*unwrap[2] : 0.0;
          v_tally(i, v);
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixAddForceXY::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == ilevel_respa) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixAddForceXY::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   potential energy of added force
------------------------------------------------------------------------- */

double FixAddForceXY::compute_scalar()
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

double FixAddForceXY::compute_vector(int n)
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

double FixAddForceXY::memory_usage()
{
  double bytes = 0.0;
  if (varflag == ATOM) bytes = maxatom*4 * sizeof(double);
  return bytes;
}
