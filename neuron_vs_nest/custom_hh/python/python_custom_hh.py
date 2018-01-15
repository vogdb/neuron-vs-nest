import numpy
import scipy as sp
import pylab as pylab
from scipy.integrate import odeint
from neuron_vs_nest import common_util

# Constants
C_m  = 200.0 # membrane capacitance, in pF
g_Na = 5000.0 # maximum conductances, in nS
g_K  = 30000.0
g_L  = 200
g_CaN  = 5000.0
g_CaL  = 10.0
g_CaK  = 30000.0
E_Na = 50.0 # Nernst reversal potentials, in mV
E_K  = -80.0
E_L  = -70.0
R = 8.314472
F = 96485.34
# E_Ca = 40.0
Temperature = 309.15
Ca_out=2.0
Ca_in_init=0.0001

# Channel gating kinetics
# Functions of membrane voltage
def alpha_m(V): return 0.4 * (V + 66.0) / (1.0 - sp.exp(-(V + 66.0)/5.0))
def beta_m(V): return 0.4 * (-(V + 32.0)) / (1.0 - sp.exp((V + 32.0)/5.0))
def h_inf(V): return 1.0 / (1.0 + sp.exp((V + 65.0)/7.0))
def h_tau(V): return 30.0 / (sp.exp((V + 60.0)/15.0) + sp.exp(-(V + 60.0)/16.0))
def n_inf(V): return 1.0 / (1.0 + sp.exp(-(V + 38.0)/15.0))
def n_tau(V): return 5.0/ (sp.exp((V + 50.0)/40.0) + sp.exp(-(V + 50.0)/50.0))
def p_inf(V): return 1.0 / (1.0 + sp.exp(-(V + 55.8)/3.7))
def p_tau(): return 400.0
def mc_inf(V): return 1.0 / (1.0 + sp.exp(-(V + 32.0)/5.0))
def mc_tau(): return 15.0
def hc_inf(V): return 1.0 / (1.0 + sp.exp((V + 50.0)/5.0))
def hc_tau(): return 50.0

#  Sodium (Na = element name)
def I_Na(V,m,h):return g_Na * m**3 * h * (V - E_Na)
#  Potassium (K = element name)
def I_K(V, n):  return g_K  * n**4     * (V - E_K)
#  Leak
def I_L(V):     return g_L             * (V - E_L)
#  Calcium N
def I_CaN(V, mc, hc, E_Ca):     return g_CaN * mc**2 * hc * (V - E_Ca)
#  Calcium L
def I_CaL(V, p, E_Ca):     return g_CaL * p * (V - E_Ca)
#  Calcium K
def I_CaK(V, Ca_in):     return 0.6 * g_CaK * ((Ca_in * Ca_in)/(Ca_in * Ca_in + 0.014 * 0.014))*(V - E_K)
#  E_Ca
def E_Ca_func(Ca_in):     return ((1000.0 * R * Temperature)/(2.0 * F))*sp.log(Ca_out/Ca_in)
# External current
def I_inj(t): return 700.0

# The time to integrate over
t = sp.arange(0.0, 150.0, 0.1)

# Integrate!
def dALLdt(X, t):
    V, m, h, n, p, mc, hc, Ca_in = X
    E_Ca = ((1000.0 * R * Temperature)/(2.0 * F))*sp.log(Ca_out/Ca_in)
    #calculate membrane potential & activation variables
    dVdt = (I_inj(t) - I_Na(V, m, h) - I_K(V, n) - I_L(V) - I_CaN(V, mc, hc, E_Ca) - I_CaL(V, p, E_Ca) - I_CaK(V, Ca_in)) / C_m
    dmdt = alpha_m(V)*(1.0-m) - beta_m(V)*m
    dhdt = (h_inf(V) - h) / h_tau(V)
    dndt = (n_inf(V) - n) / n_tau(V)
    dpdt = (p_inf(V) - p) / p_tau()
    dmcdt = (mc_inf(V) - mc) / mc_tau()
    dhcdt = (hc_inf(V) - hc) / hc_tau()
    dca_indt = 0.01 * (-0.00001 * (I_CaN(V, mc, hc, E_Ca) + I_CaL(V, p, E_Ca)) - 4.0*Ca_in)
    return dVdt, dmdt, dhdt, dndt, dpdt, dmcdt, dhcdt, dca_indt

V_Init = -65
X = odeint(dALLdt, [V_Init, alpha_m(V_Init)/(alpha_m(V_Init) + beta_m(V_Init)), h_inf(V_Init), n_inf(V_Init), p_inf(V_Init), mc_inf(V_Init), hc_inf(V_Init), Ca_in_init], t)
V = X[:,0]
m = X[:,1]
h = X[:,2]
n = X[:,3]
p = X[:,4]
mc = X[:,5]
hc = X[:,6]
Ca_in = X[:,7]
ina = I_Na(V,m,h)
ik = I_K(V, n)
il = I_L(V)
E_Ca = E_Ca_func(Ca_in)
ican = I_CaN(V,mc,hc, E_Ca)
ical = I_CaL(V,p, E_Ca)
icak = I_CaK(V,Ca_in)


pylab.figure()
pylab.title('Hodgkin-Huxley Neuron')

pylab.subplot(4,1,1)
pylab.ylabel('V (mV)')
pylab.plot(t, V, label="V_m")
pylab.legend()

pylab.subplot(4,1,2)
pylab.ylabel('Ca Concentration')
pylab.yticks(numpy.arange(0.0001, 0.0010, 0.0002))
pylab.plot(t, Ca_in, label='Ca_in')
pylab.legend()

pylab.subplot(4,1,3)
pylab.ylim(0, 1)
pylab.ylabel('m, n, h particles')
pylab.plot(t, m, 'r', label='m')
pylab.plot(t, n, 'g', label='n')
pylab.plot(t, h, 'b', label='h')
pylab.legend()

pylab.subplot(4,1,4)
pylab.ylim(0, 1)
pylab.ylabel('p, mc, hc particles')
pylab.plot(t, p, 'r', label='p')
pylab.plot(t, mc, 'g', label='mc')
pylab.plot(t, hc, 'b', label='hc')
pylab.legend()

pylab.show()