---------------
FIX WALL/REGION
---------------

Source codes:
 - fix_wall_region.cpp
 - fix_wall_region.h
 
Syntax:

	fix FIXNAME GROUP wall/region REGION lj126Expanded EPSILON SIGMA DELTALJ CUTOFF

	- FIXNAME: name of fix 
	- GROUP: group on which to apply wall force
	- REGION: region enclosed by wall
	- EPSILON: epsilon of expanded LJ force
	- SIGMA: sigma of expanded LJ force
	- DELTALJ: delta of expanded LJ force
	- CUTOFF: cutoff of expanded LJ force
	
Notes:	
 - modified from fix_wall_region.cpp
 - shifts the LJ potential by DELTALJ
 - the other wall forces described in the LAMMPS manual are unchanged
 
 
Code Snippets:

	// Shafat (5Jun20)
	else if (style == LJ126Expanded) {
		coeff1 = 48.0 * epsilon * pow(sigma,12.0);
		coeff2 = 24.0 * epsilon * pow(sigma,6.0);
		coeff3 = 4.0 * epsilon * pow(sigma,12.0);
		coeff4 = 4.0 * epsilon * pow(sigma,6.0);
		double r2inv = 1.0/((cutoff-deltaLJ)*(cutoff-deltaLJ));
		double r6inv = r2inv*r2inv*r2inv;
		offset = r6inv*(coeff3*r6inv - coeff4);
	}

	/* ----------------------------------------------------------------------
	   Shafat (5Jun20): LJ 12/6 Expanded interaction for particle with wall
	   compute eng and fwall = magnitude of wall force
	------------------------------------------------------------------------- */
	void FixWallRegion::lj126Expanded(double r)
	{
	  double rinv = 1.0/(r-deltaLJ);
	  double r2inv = rinv*rinv;
	  double r6inv = r2inv*r2inv*r2inv;
	  fwall = r6inv*(coeff1*r6inv - coeff2) * rinv;
	  eng = r6inv*(coeff3*r6inv - coeff4) - offset;
	}