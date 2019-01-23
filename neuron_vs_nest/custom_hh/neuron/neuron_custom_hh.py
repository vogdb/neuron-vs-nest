import numpy
import pylab
from neuron import h
from neuron import gui  # VERY IMPORTANT to include 'gui'! Despite the fact that it is unused.

from neuron_vs_nest import common_util

soma = h.Section(name='soma')
soma.L = 200.
soma.diam = 50. / numpy.pi
soma.nseg = 1
h.define_shape()  # Translate into 3D points.
soma.cm = 2  # Membrane capacitance in micro Farads / cm^2
soma.insert('custom_hh')

iclamp = h.IClamp(soma(0.5))
iclamp.delay = 0
iclamp.dur = common_util.RUN_TIME
iclamp.amp = 0.7  # nA

t_vec = h.Vector()
v_vec = h.Vector()
m_vec = h.Vector()
n_vec = h.Vector()
h_vec = h.Vector()
t_vec.record(h._ref_t)
v_vec.record(soma(0.5)._ref_v)
m_vec.record(soma(0.5)._ref_m_custom_hh)
n_vec.record(soma(0.5)._ref_n_custom_hh)
h_vec.record(soma(0.5)._ref_h_custom_hh)

duration = common_util.RUN_TIME
h.tstop = duration
h.run()

pylab.figure("custom_hh V_m")
pylab.plot(t_vec, v_vec)
pylab.figure("custom_hh m")
pylab.plot(t_vec, m_vec)
pylab.figure("custom_hh n")
pylab.plot(t_vec, n_vec)
pylab.figure("custom_hh h")
pylab.plot(t_vec, h_vec)

pylab.show()
