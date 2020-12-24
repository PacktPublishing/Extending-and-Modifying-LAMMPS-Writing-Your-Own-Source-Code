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

#include "fix_printEvaporate.h"
#include <mpi.h>
#include <cstring>
#include "update.h"
#include "input.h"
#include "modify.h"
#include "variable.h"
#include "memory.h"
#include "error.h"
#include "force.h"

#include "atom.h"	// SM
#include "group.h"  // SM
#include "region.h"	// SM
#include "domain.h" // SM

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixPrintEvaporate::FixPrintEvaporate(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  fp(NULL), string(NULL), copy(NULL), work(NULL), var_print(NULL), idregion(NULL) // SM
{
  if (narg < 5) error->all(FLERR,"Illegal fix printEvaporate command");
  if (strstr(arg[3],"v_") == arg[3]) {
    int n = strlen(&arg[3][2]) + 1;
    var_print = new char[n];
    strcpy(var_print,&arg[3][2]);
    nevery = 1;
  } else {
    nevery = force->inumeric(FLERR,arg[3]);
    if (nevery <= 0) error->all(FLERR,"Illegal fix printEvaporate command");
  }

  MPI_Comm_rank(world,&me);

  int n = strlen(arg[4]) + 1;
  string = new char[n];
  strcpy(string,arg[4]);

  copy = (char *) memory->smalloc(n*sizeof(char),"fix/printEvaporate:copy");
  work = (char *) memory->smalloc(n*sizeof(char),"fix/printEvaporate:work");
  maxcopy = maxwork = n;

  // parse optional args

  fp = NULL;
  screenflag = 1;
  char *title = NULL;

  iregion = -1; 	// SM

  int iarg = 5;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"file") == 0 || strcmp(arg[iarg],"append") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix printEvaporate command");
      if (me == 0) {
        if (strcmp(arg[iarg],"file") == 0) fp = fopen(arg[iarg+1],"w");
        else fp = fopen(arg[iarg+1],"a");
        if (fp == NULL) {
          char str[128];
          snprintf(str,128,"Cannot open fix printEvaporate file %s",arg[iarg+1]);
          error->one(FLERR,str);
        }
      }
      iarg += 2;
    } else if (strcmp(arg[iarg],"screen") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix printEvaporate command");
      if (strcmp(arg[iarg+1],"yes") == 0) screenflag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) screenflag = 0;
      else error->all(FLERR,"Illegal fix printEvaporate command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"title") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix printEvaporate command");
      delete [] title;
      int n = strlen(arg[iarg+1]) + 1;
      title = new char[n];
      strcpy(title,arg[iarg+1]);
      iarg += 2;
    } 

	// SM: parse region ID
	else if (strcmp(arg[iarg],"region") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix printEvaporate command");
      iregion = domain->find_region(arg[iarg+1]);
      if (iregion == -1)
        error->all(FLERR,"Region ID for fix printEvaporate does not exist");
      int n = strlen(arg[iarg+1]) + 1;
      idregion = new char[n];
      strcpy(idregion,arg[iarg+1]);
      iarg += 2;
    }
	
	else error->all(FLERR,"Illegal fix printEvaporate command");
  }

  // print file comment line

  if (fp && me == 0) {
    if (title) fprintf(fp,"%s\n",title);
    else fprintf(fp,"# Fix printEvaporate output for fix %s\n",id);
  }

  delete [] title;
}

/* ---------------------------------------------------------------------- */

FixPrintEvaporate::~FixPrintEvaporate()
{
  delete [] string;
  delete [] var_print;
  memory->sfree(copy);
  memory->sfree(work);
  delete [] idregion;	// SM

  if (fp && me == 0) fclose(fp);
}

/* ---------------------------------------------------------------------- */

int FixPrintEvaporate::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixPrintEvaporate::init()
{
  if (var_print) {
    ivar_print = input->variable->find(var_print);
    if (ivar_print < 0)
      error->all(FLERR,"Variable name for fix printEvaporate timestep does not exist");
    if (!input->variable->equalstyle(ivar_print))
      error->all(FLERR,"Variable for fix printEvaporate timestep is invalid style");
    next_print = static_cast<bigint>
      (input->variable->compute_equal(ivar_print));
    if (next_print <= update->ntimestep)
      error->all(FLERR,"Fix printEvaporate timestep variable returned a bad timestep");
  } else {
    if (update->ntimestep % nevery)
      next_print = (update->ntimestep/nevery)*nevery + nevery;
    else
      next_print = update->ntimestep;
  }
  
  // SM: check region validity
  if (iregion >= 0) {
    iregion = domain->find_region(idregion);
    if (iregion == -1)
      error->all(FLERR,"Region ID for fix printEvaporate does not exist");
  }

  // add next_print to all computes that store invocation times
  // since don't know a priori which are invoked via variables by this fix
  // once in end_of_step() can set timestep for ones actually invoked

  modify->addstep_compute_all(next_print);
}

/* ---------------------------------------------------------------------- */

void FixPrintEvaporate::setup(int /* vflag */)
{
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixPrintEvaporate::end_of_step()
{
  if (update->ntimestep != next_print) return;

  // make a copy of string to work on
  // substitute for $ variables (no printing)
  // append a newline and print final copy
  // variable evaluation may invoke computes so wrap with clear/add

  modify->clearstep_compute();

  strcpy(copy,string);
  input->substitute(copy,work,maxcopy,maxwork,0);

  if (var_print) {
    next_print = static_cast<bigint>
      (input->variable->compute_equal(ivar_print));
    if (next_print <= update->ntimestep)
      error->all(FLERR,"Fix printEvaporate timestep variable returned a bad timestep");
  } else {
    next_print = (update->ntimestep/nevery)*nevery + nevery;
  }

  modify->addstep_compute(next_print);

  // SM: update region if required
  Region *region = NULL;
  if (iregion >= 0) {
    region = domain->regions[iregion];
    region->prematch();
  }

  // SM: identify atoms located in region
  int nlocal = atom->nlocal;
  double **x = atom->x;
  int *mask = atom->mask;
  for (int i = 0; i < nlocal; i++)
	  if (mask[i] & groupbit) {
		if (region && !region->match(x[i][0],x[i][1],x[i][2])) continue;
		if (screenflag && screen) fprintf(screen,"%s\n",copy);
		if (screenflag && logfile) fprintf(logfile,"%s\n",copy);
		if (fp) {
		    fprintf(fp,"%s\n",copy);
		    fflush(fp);
		}
      }
  
}
