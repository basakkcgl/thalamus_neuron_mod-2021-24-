TITLE leaks 
 

:Leaks for HH implementation.
:Modified from HH built-in mod file in NEURON by Basak Kocaoglu.


UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
	(S) = (siemens)
}

NEURON {
        SUFFIX leaks
        NONSPECIFIC_CURRENT il
        RANGE el, gl
        :GLOBAL 
	THREADSAFE : assigned GLOBALs will be per thread
}

PARAMETER {
        gl = 0.0001 (S/cm2)	<0,1e9> :.0003 (S/cm2)	<0,1e9>
        el = -60 (mV) :-54.3 (mV)
}

ASSIGNED {
        v (mV)
        celsius (degC)
        il (mA/cm2)
}

BREAKPOINT {    
        il = gl*(v - el)
}


