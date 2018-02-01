import pylab
from neuron import h
from neuron import gui  # VERY IMPORTANT to include 'gui'! Despite the fact that it is unused.

RUN_TIME = 200
SPIKES_START = 1
SPIKES_INTERVAL = 20
SPIKES_NUMBER = 10000
SPIKE_TIMES = [i for i in range(SPIKES_START, SPIKES_NUMBER, SPIKES_INTERVAL)]

soma = h.Section(name='soma')

iaf_neuron = h.IntFire1()
iaf_neuron.tau = 1
iaf_neuron.refrac = 2

generator = h.NetStim()
generator.number = SPIKES_NUMBER # pool of available spike to emit.
generator.start = SPIKES_START
generator.interval = SPIKES_INTERVAL
generator.noise = 0

conn = h.NetCon(generator, iaf_neuron)
conn.delay = 0
conn.weight[0] = 2

v_vec = h.Vector()
t_vec = h.Vector()
t_vec.record(h._ref_t)
v_vec.record(iaf_neuron._ref_m)

duration = RUN_TIME
h.tstop = duration
h.run()

pylab.figure("Neuron IntFire")

pylab.subplot(2, 1, 1)
pylab.xlim(0, RUN_TIME)
pylab.plot(t_vec, v_vec, 'r.')
pylab.ylabel("V_m")

pylab.subplot(2, 1, 2)
pylab.xlim(0, RUN_TIME)
pylab.plot(SPIKE_TIMES, [1] * len(SPIKE_TIMES), ".")
pylab.ylabel("spikes")

pylab.show()
