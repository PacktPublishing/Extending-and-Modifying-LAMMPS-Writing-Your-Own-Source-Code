---------------
PAIR HARMONIC
---------------

Source codes:
 - pair_harmonic.cpp
 - pair_harmonic.h
 
Syntax:

	pair_style harmonic GLOBAL_CUTOFF
	pair_coeff TYPE1 TYPE2 KSP R0 E0 LOCAL_CUTOFF

	- GLOBAL_CUTOFF: cutoff for all atom pairs interacting via this potential
	- TYPE1, TYPE2: atom types that interact via this potential
	- KSP: spring constant of potential
	- R0: location of energy minimum
	- E0: depth of potential well
	- LOCAL_CUTOFF: cutoff that applies to this potential for this atom pair specifically
	
Notes:	
 - applies a harmonic potential with the minimum located at r = r0 and depth of e0
 
 
Code Snippets:

      if (rsq < cutsq[itype][jtype]) {
        r = sqrt(rsq);
        dr = r - r0[itype][jtype];
        fpair = -factor_lj * ksp[itype][jtype] * dr / r;	// SM

	---------------------------------------------------------------

       if (eflag) {
          evdwl = 0.5 * ksp[itype][jtype] * (dr*dr) - e0[itype][jtype] -	
            offset[itype][jtype];	// SM
          evdwl *= factor_lj;
        }
