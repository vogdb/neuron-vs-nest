import nest
import pylab
import common_util

neuron = nest.Create(
    "hh_psc_alpha", params={
        "C_m": 200.0,  # pF
        "t_ref": 0.0,
        "I_e": 700.0,
    })
multimeter = nest.Create(
    'multimeter',
    params={"record_from": ["V_m"], "withtime": True}
)

nest.Connect(multimeter, neuron)
nest.Simulate(common_util.RUN_TIME)

status = nest.GetStatus(multimeter)[0]
events = status['events']
times = events['times']
pylab.figure("Nest iclamp")
pylab.plot(times, events["V_m"])

pylab.show()

