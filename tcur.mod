TITLE T channel for Varela Lab
:
:
: Modified from Destexhe 1995 and Hines 2001
: by Basak Kocaoglu
: 

NEURON {
	SUFFIX tcur
	USEION ca READ eca, cai WRITE ica
		RANGE  m_inf, h_inf, tau_h, shift, actshift, tau_m
		RANGE cabuffer
	GLOBAL qm, qh
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
	gT	= 0.00033 (S/cm2) :0.33	(mS/cm2) : T soma taken as default here
	eca (mV)
	cai (mM)
}
ASSIGNED {
	ica	(mA/cm2)
	m_inf
	tau_m	(ms)
	h_inf
	tau_h	(ms)
	phi_h
	cabuffer
	celsius	(degC)
	phi_m
}
STATE {
	m
	h
}
BREAKPOINT {
	SOLVE state METHOD derivimplicit 
	eval(v)
	cabuffer = 1000 * 8.31441 * (273.15 + celsius) / (2 * FARADAY) * log(2 / cai) :((exp(1/2)) / (exp(1/cai)))
	ica = gT* m*m*h * ((v-cabuffer))
	
}
DERIVATIVE state {
	m' = m_inf :  (m_inf - m) / tau_m
	h' = (h_inf - h) / tau_h
}
INITIAL {
	m = m_inf
	h = h_inf
}
UNITSOFF
PROCEDURE eval(v(mV)) {
	m_inf = 1.0 / ( 1 + exp(-(v+shift+actshift+57)/6.2) )
	h_inf = 1.0 / ( 1 + exp((v+shift+81)/4.0) )
	:phi_h = qh ^ ((celsius-24)/10)
	phi_h = 3.73
	tau_h = (30.8 + (211.4 + exp(v+shift+113.2)/qm))/(1+exp((v+shift+84)/3.2))/phi_h
	
	phi_m = qm ^ ((celsius-24)/10)
	tau_m = ( 3 + 1.0 / ( exp((v+shift+25)/10) + exp(-(v+shift+100)/15) ) ) / phi_m
}
UNITSON
