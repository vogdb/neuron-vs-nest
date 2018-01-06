import pylab
import nest
from neuron_vs_nest import common_util

nest.Install("custom_hh_model")
# nest.SetKernelStatus(dict(resolution=0.5))

neuron = nest.Create(
    'custom_hh', params={
        "I_e": 700.0,  # pA
        "C_m": 200.0,  # pF
        "t_ref": 0.0,
    }
)
multimeter = nest.Create(
    'multimeter',
    params={"record_from": ["V_m", "Act_m", "Act_h", "Inact_n"], "withtime": True}
)

nest.Connect(multimeter, neuron)
nest.Simulate(common_util.RUN_TIME)


def plot_parameter(plot_name, device, param_to_display):
    status = nest.GetStatus(device)[0]
    events = status['events']
    times = events['times']
    pylab.figure(plot_name)
    pylab.plot(times, events[param_to_display])


plot_parameter('nest iclamp V_m', multimeter, 'V_m')
plot_parameter('nest iclamp Act_m', multimeter, 'Act_m')
plot_parameter('nest iclamp Act_h', multimeter, 'Act_h')
plot_parameter('nest iclamp Inact_n', multimeter, 'Inact_n')

pylab.show()

