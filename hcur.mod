TITLE h current for Varela Lab
:
: MODIFIED by Basak Kocaoglu
: Written by Destexhe, 1996
: 
: 

NEURON {
	SUFFIX hcur
	USEION h READ eh WRITE ih VALENCE 1
	USEION ca READ cai :WRITE cai
        	RANGE ghbar, ih_inf, tau_s, ihshift
	GLOBAL k2, cac, k4, Pc, nca, nexp, ginc, taum, qm, qh
		RANGE poneic, openic, openlockedic, alpha, beta
}

UNITS {
	(molar) = (1/liter)		
	(mM)	= (millimolar)
	(um)	= (micron)
	(mA)	= (milliamp)
	(msM)	= (ms mM)
	FARADAY = (faraday) (coulomb)
}

PARAMETER {
	qm	= 5		: q10's for activation and inactivation
	qh	= 3		: from Coulter et al., J Physiol 414: 587, 1989
	shift	= 2 	(mV)	: corresponds to 2mM ext Ca++
	actshift = 0 	(mV)	: shift of activation curve (towards hyperpol)
	v		(mV)
:below is for Ih
	eh	(mV)
	ghbar	= 2 (S/cm2)
	cac	= 0.002 (mM)		: half-activation of calcium dependence
	k2	= 0.0004 (1/ms)		: inverse of time constant
	Pc	= 0.007			: half-activation of CB protein dependence
	k4	= 0.001	(1/ms)		: backward binding on Ih
	nca	= 4			: number of binding sites of ca++
	nexp	= 1			: number of binding sites on Ih channels
	ginc	= 2			: augmentation of conductance with Ca++
	taum	= 20	(ms)		: min value of tau
	ihshift = 0
	poneic= 0.0
	openic = 0.6
	openlockedic = 0.0
	cai (mM)
}
ASSIGNED {
	ih_inf 
	tau_s	(ms)
	alpha	(1/ms)
	beta	(1/ms)
	ih	(mV)
	celsius (degC)
	ncaexp
	nexpexp
}
STATE {
	pone
	open
	openlocked
}
BREAKPOINT {
	SOLVE state METHOD derivimplicit
	eval(v)
	ih = ghbar * (open + ginc * openlocked) * (v-eh) : Varela Dynasim Model implementation
	:ih = ghbar * (open + ginc *(1 - openlocked - open)) * (v-eh) : Soplata 2017 original equation
}
DERIVATIVE state {
	ncaexp = (cai/cac) * (cai/cac) * (cai/cac) * (cai/cac) : in parameter block it is assigned as 4
	pone' = k2* ncaexp * ((1-pone) - k2 * pone) 
	open' = alpha * (1 - open - openlocked) - (beta * open)
	nexpexp = nexp : in parameter block it is assigned as 1
	openlocked' = (k4 * (pone/Pc) * nexpexp) * open - k4 * openlocked
}


INITIAL {
	eval(v)
	ih = ih_inf 

	pone = poneic
	open = openic
	openlocked = openlockedic

}

UNITSOFF
PROCEDURE eval(v (mV)) {
	ih_inf = 1 / ( 1 + exp((v+75-ihshift)/5.5) )
	tau_s = (taum + 1000 / ( exp((v+71.5-ihshift)/14.2) + exp(-(v+89-ihshift)/11.6) ) )
	alpha = ih_inf / tau_s
	beta  = ( 1 - ih_inf ) / tau_s
}

UNITSON

