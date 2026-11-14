TITLE decay of internal calcium concentration
:MODIFIED BY Basak Kocaoglu

: This mechanism was published in:  Destexhe, A. Babloyantz, A. and 
: Sejnowski, TJ.  Ionic mechanisms for intrinsic slow oscillations in
: thalamic relay neurons. Biophys. J. 65: 1538-1552, 1993)
:
: Written by Alain Destexhe, Salk Institute, Nov 12, 1992
:

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX cabuf
	USEION ca READ ica, cai WRITE cai
	GLOBAL cainf,taur
}

UNITS {
	(molar) = (1/liter)			: moles do not appear in units
	(mM)	= (millimolar)
	(um)	= (micron)
	(mA)	= (milliamp)
	(msM)	= (ms mM)
	FARADAY = (faraday) (coulomb)
}


PARAMETER {
	taur	= 5	(ms)		: rate of calcium removal
	cainf	= 2e-4	(mM)
}

STATE {
	cai		(mM) <1e-6>
}

INITIAL {
	cai = cainf	:really should be done in hoc (initialization
			: of states is consistent but not currents from
			: fcurrent() since eca comes from nernst with
			: cai0_ca_ion and cao0_ca_ion first and then
			: this initial block is called.

}

ASSIGNED {
	ica		(mA/cm2)
	drive_channel	(mM/ms)
	A (mM/cm2)
}
	
BREAKPOINT {
	SOLVE state METHOD derivimplicit
}

DERIVATIVE state { 

	A = 10000/(2*FARADAY) :10/

	drive_channel = - A * ica

	if (drive_channel <= 0.) { drive_channel = 0. }	: cannot pump inward

	cai' = drive_channel + (cainf-cai)/taur
}
