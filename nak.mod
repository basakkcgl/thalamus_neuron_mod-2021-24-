TITLE Na+ and K+ channel for Varela Lab
:
: Modified from Destexhe, Houweling by Basak Kocaoglu.
:
: 
: Update = vtraub changed as vshiftna


INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX nak
	USEION na READ ena WRITE ina
	USEION k READ ek WRITE ik
	RANGE gnabar, vshiftna
	RANGE m_inf, h_inf, am,bm,ah,bh,vdiff, n_inf
	RANGE ik,ina 
	RANGE gkbar, vshiftk, gKleak
	RANGE tau_m, tau_h, tau_n
	THREADSAFE : assigned GLOBALs will be per thread
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
        (S) = (siemens)
	:(mS) = (millisiemens) : NEURON seems to have an issue of overriding default parameter units. Make changes at interpreter level. 
}

PARAMETER {
	gnabar	= 0.09 	(S/cm2)
	:gnabar	= 90 	(mS/cm2) : NEURON seems to have an issue of overriding default parameter units. Make changes at interpreter level. 
	ena		(mV)
	v               (mV)
	vshiftna	= -40	(mV)
	vshiftk	= -25	(mV)
	gkbar	= 0.01	(S/cm2)
	:gkbar	= 10	(mS/cm2) : NEURON seems to have an issue of overriding default parameter units. Make changes at interpreter level. 
	ek		(mV)
	gKleak = 0.000014 (S/cm2)
	eKleak = -105 (mV)
	
}

STATE {
	m h n
}

ASSIGNED {
	ina	(mA/cm2)
	ik	(mA/cm2)
	m_inf
	h_inf
	n_inf
	am
	bm
	ah
	bh
	an
	bn
	vdiffna
	vdiffk
	tau_m
	tau_h
	tadj
	tau_n
	celsius (degC)
}


BREAKPOINT {

	SOLVE states METHOD derivimplicit
	ina = gnabar * m*m*m*h * (v - ena)
	ik  = gkbar * n*n*n*n * (v - ek) - gKleak * (v- eKleak) : Potassium leak subtracted 
}

DERIVATIVE states {   : use this for exact Hodgkin-Huxley equations
	evaluate_fct(v)
	m' = (m_inf - m)/ tau_m
	h' = (h_inf - h)/ tau_h
	n' = (n_inf - n)/ tau_n
	
	:m'= am* - am * m - bm
	:h' = ah * - ah * h - bh
}
UNITSOFF
INITIAL {
	tadj = 3.0 ^ ((celsius-6.3)/ 10 )
	evaluate_fct(v)
	m= m_inf
	h= h_inf
	n= n_inf
}
PROCEDURE evaluate_fct(v(mV)) {
	vdiffna = v - vshiftna
	vdiffk = v - vshiftk
	am = 0.32 * (13-vdiffna) / ( exp((13-vdiffna)/4) - 1)
	bm = 0.28 * (vdiffna-40) / ( exp((vdiffna-40)/5) - 1)
	ah = 0.128 * exp((17-vdiffna)/18)
	bh = 4 / ( 1 + exp((40-vdiffna)/5) )
	an = 0.032 * (15-vdiffk) / ( exp((15-vdiffk)/5) - 1)
	bn = 0.5 * exp((10-vdiffk)/40)
	n_inf = an / (an + bn)
	m_inf = am / (am + bm)
	h_inf = ah / (ah + bh)
	tau_m = 1 / (am + bm) / tadj
	tau_n = 1 / (an + bn) / tadj
	tau_h = 1 / (ah + bh) / tadj
}
UNITSON
