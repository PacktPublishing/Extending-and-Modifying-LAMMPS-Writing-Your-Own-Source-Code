------------------
FIX PRINTEVAPORATE
------------------

Source codes:
 - fix_printEvaporate.cpp
 - fix_printEvaporate.h
 
Syntax:

	fix FIXNAME GROUP printEvaporate N region R keywords

	- FIXNAME: name of fix 
	- GROUP: group on which fix acts
	- N: frequency of printing
	- R: region where atom must be located to print
	
Notes:	
 - modified from fix_print.cpp
 - prints output when atom in group is located in a given regionR
 
Code Snippets:

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