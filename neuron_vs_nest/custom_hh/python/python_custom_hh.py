import scipy as sp
import pylab as pylab
from scipy.integrate import odeint
from neuron_vs_nest import common_util

# Constants
C_m  = 200.0 # membrane capacitance, in pF
g_Na = 5000.0 # maximum conductances, in nS
g_K  = 30000.0
g_L  = 200
E_Na = 50.0 # Nernst reversal potentials, in mV
E_K  = -80.0
E_L  = -70.0

# Channel gating kinetics
# Functions of membrane voltage
def alpha_m(V): return 0.4 * (V + 66.0) / (1.0 - sp.exp(-(V + 66.0)/5.0))
def beta_m(V): return 0.4 * (-(V + 32.0)) / (1.0 - sp.exp((V + 32.0)/5.0))
def h_inf(V): return 1.0 / (1.0 + sp.exp((V + 65.0)/7.0))
def h_tau(V): return 30.0 / (sp.exp((V + 60.0)/15.0) + sp.exp(-(V + 60.0)/16.0))
def n_inf(V): return 1.0 / (1.0 + sp.exp(-(V + 38.0)/15.0))
def n_tau(V): return 5.0/ (sp.exp((V + 50.0)/40.0) + sp.exp(-(V + 50.0)/50.0))

#  Sodium (Na = element name)
def I_Na(V,m,h):return g_Na * m**3 * h * (V - E_Na)
#  Potassium (K = element name)
def I_K(V, n):  return g_K  * n**4     * (V - E_K)
#  Leak
def I_L(V):     return g_L             * (V - E_L)
# External current
def I_inj(t): return 700.0

# The time to integrate over
t = sp.arange(0.0, common_util.RUN_TIME, 0.1)

# Integrate!
def dALLdt(X, t):
    V, m, h, n = X
    #calculate membrane potential & activation variables
    dVdt = (I_inj(t) - I_Na(V, m, h) - I_K(V, n) - I_L(V)) / C_m
    dmdt = alpha_m(V)*(1.0-m) - beta_m(V)*m
    dhdt = (h_inf(V) - h) / h_tau(V)
    dndt = (n_inf(V) - n) / n_tau(V)
    return dVdt, dmdt, dhdt, dndt

V_Init = -65
X = odeint(dALLdt, [V_Init, alpha_m(V_Init)/(alpha_m(V_Init) + beta_m(V_Init)), h_inf(V_Init), n_inf(V_Init)], t)
V = X[:,0]
m = X[:,1]
h = X[:,2]
n = X[:,3]
ina = I_Na(V,m,h)
ik = I_K(V, n)
il = I_L(V)


pylab.figure("python V_m")
pylab.plot(t, V)

pylab.figure("python m")
pylab.plot(t, m)
pylab.figure("python n")
pylab.plot(t, n)
pylab.figure("python h")
pylab.plot(t, h)

pylab.show()