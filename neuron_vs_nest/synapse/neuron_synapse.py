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
soma.insert('hh')

syn_exc = h.ExpSyn(0.5, sec=soma)
syn_exc.e = 0.0
syn_exc.tau = 0.8

generator = h.NetStim()
generator.number = 1000000 # pool of available spike to emit.
generator.start = common_util.SPIKES_START
generator.interval = common_util.SPIKES_INTERVAL
generator.noise = 0

conn = h.NetCon(generator, syn_exc)
conn.delay = 0
conn.weight[0] = 1

v_vec = h.Vector()
t_vec = h.Vector()
v_vec.record(soma(0.5)._ref_v)
t_vec.record(h._ref_t)

duration = common_util.RUN_TIME
h.tstop = duration
h.run()

pylab.figure("neuron synapse")
pylab.plot(t_vec, v_vec)

pylab.show()
