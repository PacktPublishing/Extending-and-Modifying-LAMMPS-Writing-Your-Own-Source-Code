---------------
FIX NH
---------------

Source codes:
 - fix_nh.cpp
 - fix_nh.h
 
Syntax:

	fix FIXNAME GROUP nvt tempExp TSTART TPERIOD BETA BOOST OUTPUTFLAG

	- FIXNAME: name of fix 
	- GROUP: group on which thermostat acts
	- TSTART: starting temperature
	- TPERIOD: temperature damping parameter
	- BETA: heating rate
	- BOOST: boost potential (can be a variable)
	- OUTPUTFLAG: 0 or 1 to indicate output of energy or target temperature
	
Notes:	
 - modified from fix_nh.cpp
 - increments target temperature exponentially based on boost potential and temperature 
 
Code Snippets:

  if (delta != 0.0) delta /= update->endstep - update->beginstep;

  // SM: calculate t_target without boost
  if (tempExpFlag == 0)
	t_target = t_start + delta * (t_stop-t_start);
 
  // SM: calculate t_target with boost
  if (tempExpFlag == 1) {
	modify->clearstep_compute();
	if (boostStyle == EQUAL) 
		boostV = input->variable->compute_equal(boostVar);
	
	if (delta == 0) t_target = t_start;
	double kBT_Vb = boostV/(boltz*t_target);
	t_target += betaTemp * dtv * pow(M_E,kBT_Vb);	
		
	modify->addstep_compute(update->ntimestep + 1);
  }