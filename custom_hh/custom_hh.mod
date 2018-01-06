:SOMA

: vogdb & Marco Capogrosso & Emanuele Formento
:
:
: This model has been adapted and is described in detail in:
:
: McIntyre CC and Grill WM. Extracellular Stimulation of Central Neurons:
: Influence of Stimulus Waveform and Frequency on Neuronal Output
: Journal of Neurophysiology 88:1592-1604, 2002.

TITLE Motor Axon Soma
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX custom_hh
	NONSPECIFIC_CURRENT ina
	NONSPECIFIC_CURRENT ikrect
	NONSPECIFIC_CURRENT il
	RANGE gnabar, gl, ena, ek, el, gkrect
	RANGE m_inf, h_inf, n_inf
	RANGE tau_m, tau_h, tau_n
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	:SOMA PARAMETERS
	gnabar	= 0.05	(mho/cm2)
	gl	= 0.002 (mho/cm2)
	gkrect = 0.3  (mho/cm2)

	ena     = 50.0  (mV)
	ek      = -80.0 (mV)
	el	= -70.0 (mV)

	dt              (ms)
	v               (mV)

	amA = 0.4
	amB = 66
	amC = 5
	bmA = 0.4
	bmB = 32
	bmC = 5
	R=8.314472
	F=96485.34
}

STATE {
	 m h n
}

ASSIGNED {
	ina	 (mA/cm2)
	il      (mA/cm2)
	ikrect    (mA/cm2)
	m_inf
	h_inf
	n_inf
	tau_m
	tau_h
	tau_n
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ina = gnabar * m*m*m*h*(v - ena)
	ikrect   = gkrect *n*n*n*n*(v - ek)
	il   = gl * (v - el)
}

DERIVATIVE states {  
	 : exact Hodgkin-Huxley equations
    evaluate_fct(v)
	m' = (m_inf - m) / tau_m
	h' = (h_inf - h) / tau_h
	n' = (n_inf - n) / tau_n
}

UNITSOFF

INITIAL {
	evaluate_fct(v)
	m = m_inf
	h = h_inf
	n = n_inf
}

PROCEDURE evaluate_fct(v(mV)) { LOCAL a,b,v2

	:FAST SODIUM
	:m
	a = alpham(v)
	b = betam(v)
	tau_m = 1 / (a + b)
	m_inf = a / (a + b)
	:h
	tau_h = 30 / (Exp((v+60)/15) + Exp(-(v+60)/16))
	h_inf = 1 / (1 + Exp((v+65)/7))

	
	:DELAYED RECTIFIER POTASSIUM 
	tau_n = 5 / (Exp((v+50)/40) + Exp(-(v+50)/50))
	n_inf = 1 / (1 + Exp(-(v+38)/15))
}


FUNCTION alpham(x) {
	if (fabs((x+amB)/amC) < 1e-6) {
		alpham = amA*amC
	}else{
		alpham = (amA*(x+amB)) / (1 - Exp(-(x+amB)/amC))
	}
}



FUNCTION betam(x) {
	if (fabs((x+bmB)/bmC) < 1e-6) {
		betam = -bmA*bmC
	}else{
		betam = (bmA*(-(x+bmB))) / (1 - Exp((x+bmB)/bmC))
	}
}

FUNCTION Exp(x) {
	if (x < -100) {
		Exp = 0
	}else{
		Exp = exp(x)
	}
}

UNITSON
