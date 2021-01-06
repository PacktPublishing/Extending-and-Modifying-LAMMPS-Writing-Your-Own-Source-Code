---------------
FIX ADDBOOST
---------------

Source codes:
 - fix_addboost.cpp
 - fix_addboost.h
 
Syntax:

	fix FIXNAME GROUP addboost GROUP2 ATOMID POTENTIAL P1 P2 P3 CUTOFF

	- FIXNAME: name of fix 
	- GROUP: placeholder
	- GROUP2: group with which the atom to be boosted interacts  
	- ATOMID: compute that returns the atom ID to be boosted
	- POTENTIAL: potential type LJ, morse or harmonic
	- P1, P2, P3: potential parameters
	- CUTOFF: cutoff of potential
	
Notes:	
 - modified from fix_addforce.cpp
 - identifies an atom from the compute and applies a boost potential of the form specified
 - applies equal and opposite forces to the neighbouring atoms as well
 
Code Snippets:
