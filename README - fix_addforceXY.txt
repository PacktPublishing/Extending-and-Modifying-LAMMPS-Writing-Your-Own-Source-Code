------------------
PAIR CUNDALLSTRACK
------------------

Source codes:
 - pair_CundallStrack.cpp
 - pair_CundallStrack.h
 
Syntax:

	pair_style CundallStrack Fn Kt MuS A MuK DAMPFLAG
	pair_coeff TYPE1 TYPE2 
	
	- Fn: normal force between spheres
	- Kt: spring constant of tangential spring that applies friction
	- MuS: coefficient fo static friction
	- A: drag coefficient of velocity dependent tangential force
	- MuK: coefficient of kinetic friction that determines upper limit of tangential force
	- DAMPFLAG: 0 or 1 activate tangential force
	- TYPE1, TYPE2: atom types on which potential acts

Notes:	
 - applies a Cundall-Strack scheme tangential friction on rotating spheres 
 
 
Code Snippets:

        fsmax = MuK*FN;		// SM [8Nov20]

        if (fs > MuS*FN) {			// SM [8Nov20]
          if (shrmag != 0.0) {
			  
			shrmaginv = 1.0/shrmag;		// SM [8Nov20]  
			  
			sgn = 1;
			if (fs1<0) sgn=-1;
            shear[0] = -(shear[0]*shrmaginv) * fsmax/kt * sgn;		// SM [8Nov20]

			sgn = 1;
			if (fs2<0) sgn=-1;
            shear[1] = -(shear[1]*shrmaginv) * fsmax/kt * sgn;		// SM [8Nov20]
			
			sgn = 1;
			if (fs3<0) sgn=-1;
			shear[2] = -(shear[2]*shrmaginv) * fsmax/kt * sgn;		// SM [8Nov20]
            
			fs1 *= fsmax/fs;	// SM [8Nov20]
            fs2 *= fsmax/fs;	// SM [8Nov20]
            fs3 *= fsmax/fs;	// SM [8Nov20]
			
          } else fs1 = fs2 = fs3 = 0.0;
        }

		fprintf(screen,"%f %f %f \n", vtr2, fs2, shear[1]);	// SM [8Nov20]

        // forces & torques

        fx = delx*rinv*FN + fs1;	// SM [8Nov20]
        fy = dely*rinv*FN + fs2;	// SM [8Nov20]
        fz = delz*rinv*FN + fs3;	// SM [8Nov20]
        f[i][0] += fx;
        f[i][1] += fy;
        f[i][2] += fz;

