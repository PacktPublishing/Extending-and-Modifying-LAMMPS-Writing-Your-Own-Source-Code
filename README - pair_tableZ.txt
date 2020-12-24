---------------
PAIR TABLEZ
---------------

Source codes:
 - pair_table_z.cpp
 - pair_table_z.h
 
Syntax:

	pair_style tableZ TABLE_STYLE N KEYWORD
	pair_coeff  TYPE1 TYPE2 FILE_NAME  VARIABLE  threshold1 threshold2 a0 a1 CUTOFF

	- TABLE_STYLE: table style lookup, linear, spline, bitmap
	- N: number of data points in table
	- KEYWORD: long-range solver type
	- TYPE1, TYPE2: atom types that interact via this potential
	- FILE_NAME: file name containing table
	- VARIABLE: identifying text from table file
	- threshold1: lower bound of decay region
	- threshold2: upper bound of decay region
	- a0, a1: coefficients of decay function (a0 + a1*zz^2)
	- CUTOFF: cutoff
	
Notes:	
 - applies a height-dependent quadratic decay function f(z) = a0 + a1*(z-z1)^2 to the potential between atoms within bounds 
 
 
Code Snippets:

		// -------------------------------------------------------------
		// SM: define polynomial decay for particle(s) within thresholds
		// no interaction within same molecule
		if (molecule[i] == molecule[j]) {	
			scalingZ = 0.0;
			deltaFz = 0.0;
		}
		
		// no decay interaction if both molecules are below threshold1
		else if (x[j][2]<=threshold1 && ztmp<=threshold1) { 	  
			scalingZ = a0;		
			deltaFz = 0.0;
		}
		
		// if atleast one molecule is above threshold2:
		//	- no decay force
		//  - set decay potential to the value at the upper threshold
		else if (x[j][2]>=threshold2 || ztmp>=threshold2) { 
			zz = threshold2 - threshold1;	// 28Mar15
			scalingZ = a0 + a1*zz*zz;
			deltaFz = 0.0;
		}
		
		// if atleast one molecule is located in between thresholds (and not above threshold2):
		//  - decay force  = -df/dz
		//  - decay potential = a0 + a1*zz^2
		else {				      
			if (ztmp >= x[j][2]) maxZ = ztmp;
			if (ztmp < x[j][2]) maxZ = x[j][2];	

			zz = maxZ - threshold1;		
			scalingZ = a0 + a1*zz*zz;	// f(z)
			dfdz = 2*a1*zz;				// df/dz

			// force_z = -df/dz for molecule located farther from surface
			// force_z = +df/dz for molecule located closer to surface
			if (ztmp >= x[j][2]) deltaFz = -dfdz; 	
			if (ztmp < x[j][2]) deltaFz = dfdz; 
		}	
		
		// -------------------------------------------------------------
		
		// SM: modify fx,fy,fz to incorporate scalingZ and deltaFz
        f[i][0] += scalingZ*delx*fpair;
        f[i][1] += scalingZ*dely*fpair;
        f[i][2] += scalingZ*delz*fpair + deltaFz;
        if (newton_pair || j < nlocal) {
          f[j][0] -= scalingZ*delx*fpair;
          f[j][1] -= scalingZ*dely*fpair;
          f[j][2] -= scalingZ*delz*fpair + deltaFz;
        }

