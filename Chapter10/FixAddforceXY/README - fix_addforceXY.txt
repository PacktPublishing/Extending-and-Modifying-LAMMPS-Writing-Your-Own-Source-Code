------------------
FIX ADDFORCEXY
------------------

Source codes:
 - fix_addforceXY.cpp
 - fix_addforceXY.h
 
Syntax:

	fix FIXNAME GROUP addforceXY X0 Y0 Z0 KX KY 
	
	- FIXNAME: name of fix
	- GROUP: group on which fix acts
	- X0: x-coordinate of point towards which force is directed
	- Y0: y-coordinate of point towards which force is directed
	- Z0: placeholder quantity, ignored in fix
	- KX: spring constant of force along x-direction
	- KY: spring constant of force along y-direction

Notes:	
 - Applies a 2D force [-KX*(x-X0), -KY*(y-Y0)] on group atoms
 
 
