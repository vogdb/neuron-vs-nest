from neuron import h
from neuron import gui  # VERY IMPORTANT to include 'gui'! Despite the fact that it is unused.
import pylab
import numpy
import common_util

soma = h.Section(name='soma')
soma.L = 200.
soma.diam = 50. / numpy.pi
soma.nseg = 1
h.define_shape()  # Translate into 3D points.
soma.cm = 2  # Membrane capacitance in micro Farads / cm^2
soma.insert('hh')

iclamp = h.IClamp(soma(0.5))
iclamp.delay = 0
iclamp.dur = common_util.RUN_TIME
iclamp.amp = 0.7  # nA

v_vec = h.Vector()
t_vec = h.Vector()
v_vec.record(soma(0.5)._ref_v)
t_vec.record(h._ref_t)

duration = common_util.RUN_TIME
h.tstop = duration
h.run()

pylab.figure("neuron iclamp")
pylab.plot(t_vec, v_vec)

pylab.show()
