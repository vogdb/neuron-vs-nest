import nest
import pylab
import common_util

neuron = nest.Create(
    "hh_psc_alpha", params={
        "C_m": 200.0,  # pF
        "t_ref": 0.0,
        'tau_syn_ex':2.0,
        # 'tau_syn_in':1.,
    })
multimeter = nest.Create(
    'multimeter',
    params={"record_from": ["V_m"], "withtime": True}
)

spike_times = [common_util.SPIKES_INTERVAL * i + common_util.SPIKES_START + 1 for i in range(common_util.SPIKES_NUMBER)]
generator = nest.Create(
    'spike_generator',
    params={"start": common_util.SPIKES_START, "spike_times": spike_times}
)

nest.CopyModel("static_synapse","excitatory",{"weight":100.,  "delay":0.5})
# nest.CopyModel("static_synapse","inhibitory",{"weight":200., "delay":0.5})

nest.Connect(multimeter, neuron)
# nest.Connect(generator, neuron, syn_spec={
#     "delay": 0.0,
#     "weight": 1000.0,
#     "model": "static_synapse"
# })
nest.Connect(generator, neuron, syn_spec="excitatory")
nest.Simulate(common_util.RUN_TIME)

status = nest.GetStatus(multimeter)[0]
events = status['events']
times = events['times']
pylab.figure("nest synapse")
pylab.plot(times, events["V_m"])

pylab.show()

