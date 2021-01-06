---------------
FIX ACTIVEFORCE
---------------

Source codes:
 - fix_activeforce.cpp
 - fix_activeforce.h
 
Syntax:

	fix FIXNAME GROUP activeforce FORCE DIRECTION

	- FIXNAME: name of fix 
	- GROUP: group on which to apply active force
	- FORCE: absolute value of force applied on group; can be constant or variable (i.e. v_...)
	- DIRECTION: forward (direction of increasing atom ID's) or reverse (direction of decreasing atom ID's)
	
Notes:	
 - modified from fix_addforce.cpp
 - The absolute value of the entry "FORCE" is used
 - applies force on each bonded atom in molecule except the last atom in chain
 - if direction is "forward", force is applied in the same order that bonds are listed (e.g. bond 1 2 => force on atom 2)
 - if direction is "reverse", force is applied in the reversed order that bonds are listed (e.g. bond 1 2 => force on atom 1)
 - direction of force points from first atom to second atom  
 
 
Code Snippet:

			// -------------------------------------------------------------
			// SM: Identifying bonds and direction of active force 
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
				
				if (direction == 1) atomIndexI = atom1;
				if (direction == -1) atomIndexI = atom2;
			// -------------------------------------------------------------