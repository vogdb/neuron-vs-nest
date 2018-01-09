/*
*  custom_hh.h
*
*  This file is part of NEST.
*
*  Copyright (C) 2004 The NEST Initiative
*
*  NEST is free software: you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation, either version 2 of the License, or
*  (at your option) any later version.
*
*  NEST is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with NEST.  If not, see <http://www.gnu.org/licenses/>.
*
*/
#ifndef CUSTOM_HH
#define CUSTOM_HH

#include "config.h"

#ifdef HAVE_GSL

// External includes:
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

// forwards the declaration of the function
/**
 * Function computing right-hand side of ODE for GSL solver.
 * @note Must be declared here so we can befriend it in class.
 * @note Must have C-linkage for passing to GSL. Internally, it is
 *       a first-class C++ function, but cannot be a member function
 *       because of the C-linkage.
 * @note No point in declaring it inline, since it is called
 *       through a function pointer.
 * @param void* Pointer to model neuron instance.
 */
extern "C" inline int custom_hh_dynamics( double, const double y[], double f[], void* pnode );

// Includes from nestkernel:
#include "archiving_node.h"
#include "connection.h"
#include "event.h"
#include "nest_types.h"
#include "ring_buffer.h"
#include "universal_data_logger.h"


// Includes from sli:
#include "dictdatum.h"

/* BeginDocumentation
  Name: custom_hh.

Description:

Parameters:
  The following parameters can be set in the status dictionary.
  t_ref      [ms   ] Refractory period
  g_Na       [1 / GOhm] Sodium peak conductance
  g_L        [1 / GOhm] Leak conductance
  g_K        [1 / GOhm] Potassium peak conductance
  C_m        [pF   ] Membrane capacitance
  tau_syn_ex [ms   ] Rise time of the excitatory synaptic alpha function i
  tau_syn_in [ms   ] Rise time of the inhibitory synaptic alpha function
  I_e        [pA   ] Constant Current in pA

State variables:
  These variables can be also set in the status dictionary.
  r          [integer] number of steps in the current refractory phase

  V_m        [mV   ] Membrane potential

References: Empty

Sends: nest::SpikeEvent

Receives: Spike, Current,  DataLoggingRequest
*/
class custom_hh : public nest::Archiving_Node
{
public:
  /**
  * The constructor is only used to create the model prototype in the model manager.
  */
  custom_hh();

  /**
  * The copy constructor is used to create model copies and instances of the model.
  * @node The copy constructor needs to initialize the parameters and the state.
  *       Initialization of buffers and interal variables is deferred to
  *       @c init_buffers_() and @c calibrate().
  */
  custom_hh(const custom_hh&);

  /**
  * Releases resources.
  */
  ~custom_hh();

  /**
   * Import sets of overloaded virtual functions.
   * @see Technical Issues / Virtual Functions: Overriding, Overloading, and
   * Hiding
   */
  using nest::Node::handles_test_event;
  using nest::Node::handle;

  /**
  * Used to validate that we can send nest::SpikeEvent to desired target:port.
  */
  nest::port send_test_event(nest::Node& target, nest::rport receptor_type, nest::synindex, bool);

  /**
  * @defgroup mynest_handle Functions handling incoming events.
  * We tell nest that we can handle incoming events of various types by
  * defining @c handle() and @c connect_sender() for the given event.
  * @{
  */
  void handle(nest::SpikeEvent &);        /** accept spikes */
  void handle(nest::CurrentEvent &);      /** accept input current */
  void handle(nest::DataLoggingRequest &);/** allow recording with multimeter */

  nest::port handles_test_event(nest::SpikeEvent&, nest::port);
  nest::port handles_test_event(nest::CurrentEvent&, nest::port);
  nest::port handles_test_event(nest::DataLoggingRequest&, nest::port);
  /** @} */

  // SLI communication functions:
  void get_status(DictionaryDatum &) const;
  void set_status(const DictionaryDatum &);

private:

  /** Reset parameters and state of neuron. */
  void init_state_(const Node& proto);

  /** Reset internal buffers of neuron. */
  void init_buffers_();

  /** Initialize auxiliary quantities, leave parameters and state untouched. */
  void calibrate();

  /** Take neuron through given time interval */
  void update(nest::Time const &, const long, const long);

  // The next two classes need to be friends to access the State_ class/member
  friend class nest::RecordablesMap<custom_hh>;
  friend class nest::UniversalDataLogger<custom_hh>;

  /**
  * Free parameters of the neuron.
  *
  * 
  *
  * These are the parameters that can be set by the user through @c SetStatus.
  * They are initialized from the model prototype when the node is created.
  * Parameters do not change during calls to @c update() and are not reset by
  * @c ResetNetwork.
  *
  * @note Parameters_ need neither copy constructor nor @c operator=(), since
  *       all its members are copied properly by the default copy constructor
  *       and assignment operator. Important:
  *       - If Parameters_ contained @c Time members, you need to define the
  *         assignment operator to recalibrate all members of type @c Time . You
  *         may also want to define the assignment operator.
  *       - If Parameters_ contained members that cannot copy themselves, such
  *         as C-style arrays, you need to define the copy constructor and
  *         assignment operator to copy those members.
  */
  struct Parameters_
  {
      
double t_ref; //! Refractory period
      
double g_Na; //! Sodium peak conductance
      
double g_L; //! Leak conductance
      
double g_K; //! Potassium peak conductance
      
double C_m; //! Membrane capacitance
      
double E_Na; 
      
double E_K; 
      
double E_L; 
      
double tau_syn_ex; //! Rise time of the excitatory synaptic alpha function i
      
double tau_syn_in; //! Rise time of the inhibitory synaptic alpha function
      
double I_e; //! Constant Current in pA

    double __gsl_error_tol;

    /** Initialize parameters to their default values. */
    Parameters_();
  };

  /**
  * Dynamic state of the neuron.
  *
  * 
  *
  * These are the state variables that are advanced in time by calls to
  * @c update(). In many models, some or all of them can be set by the user
  * through @c SetStatus. The state variables are initialized from the model
  * prototype when the node is created. State variables are reset by @c ResetNetwork.
  *
  * @note State_ need neither copy constructor nor @c operator=(), since
  *       all its members are copied properly by the default copy constructor
  *       and assignment operator. Important:
  *       - If State_ contained @c Time members, you need to define the
  *         assignment operator to recalibrate all members of type @c Time . You
  *         may also want to define the assignment operator.
  *       - If State_ contained members that cannot copy themselves, such
  *         as C-style arrays, you need to define the copy constructor and
  *         assignment operator to copy those members.
  */
  struct State_
  {
      /** Symbolic indices to the elements of the state vector y */
      enum StateVecElems
      {
          /** Membrane potential */
          V_m,
          
          Act_m,
          
          Act_h,
          
          Inact_n,
          
          I_syn_in__1,
          
          I_syn_in,
          
          I_syn_ex__1,
          
          I_syn_ex,
        STATE_VEC_SIZE
      };
      /** state vector, must be C-array for GSL solver */
      double ode_state[ STATE_VEC_SIZE ];
        
long r; //! number of steps in the current refractory phase
        
double alpha_m_init; 
        
double beta_m_init; 
        
double h_inf_init; 
        
double n_inf_init; 

      
double I_syn_exc; 
      
double I_syn_inh; 
      
double I_Na; 
      
double I_K; 
      
double I_L; 
      
double n_inf; 
      
double n_tau; 
      
double alpha_m; 
      
double beta_m; 
      
double h_inf; 
      
double h_tau; 

    State_();
  };

  /**
  * Internal variables of the neuron.
  *
  * 
  *
  * These variables must be initialized by @c calibrate, which is called before
  * the first call to @c update() upon each call to @c Simulate.
  * @node Variables_ needs neither constructor, copy constructor or assignment operator,
  *       since it is initialized by @c calibrate(). If Variables_ has members that
  *       cannot destroy themselves, Variables_ will need a destructor.
  */
  struct Variables_ {
      
long RefractoryCounts; //! refractory time in steps
  };

  /**
    * Buffers of the neuron.
    * Usually buffers for incoming spikes and data logged for analog recorders.
    * Buffers must be initialized by @c init_buffers_(), which is called before
    * @c calibrate() on the first call to @c Simulate after the start of NEST,
    * ResetKernel or ResetNetwork.
    * @node Buffers_ needs neither constructor, copy constructor or assignment operator,
    *       since it is initialized by @c init_nodes_(). If Buffers_ has members that
    *       cannot destroy themselves, Buffers_ will need a destructor.
    */
  struct Buffers_ {
    Buffers_(custom_hh&);
    Buffers_(const Buffers_ &, custom_hh&);

    /** Logger for all analog data */
    nest::UniversalDataLogger<custom_hh> logger_;

        // spike buffers
        inline nest::RingBuffer& get_spikeInh() {return spikeInh; }
        //!< Buffer incoming [0, 0, 0, 0, 0, 0, 1, -12]is through delay, as sum
nest::RingBuffer spikeInh;
        double spikeInh_grid_sum_;
        // spike buffers
        inline nest::RingBuffer& get_spikeExc() {return spikeExc; }
        //!< Buffer incoming [0, 0, 0, 0, 0, 0, 1, -12]is through delay, as sum
nest::RingBuffer spikeExc;
        double spikeExc_grid_sum_;


      // current buffers
      //!< Buffer incoming [0, 0, 0, 0, 0, 0, 1, -12]is through delay, as sum
nest::RingBuffer currents;
      inline nest::RingBuffer& get_currents() {return currents; }
      double currents_grid_sum_;

      // GSL ODE stuff

      /** stepping function */
      gsl_odeiv_step* __s;
      /** adaptive stepsize control function */
      gsl_odeiv_control* __c;
      /** evolution function */
      gsl_odeiv_evolve* __e;
      /** struct describing system */
      gsl_odeiv_system __sys;

      // IntergrationStep_ should be reset with the neuron on ResetNetwork,
      // but remain unchanged during calibration. Since it is initialized with
      // step_, and the resolution cannot change after nodes have been created,
      // it is safe to place both here.
      /** step size in ms */
      double __step;
      /** current integration time step, updated by GSL */
      double __integration_step;
  };

    

  /** returns number of steps in the current refractory phase in integer */
  inline long get_r() const {
    return S_. r;
  }

  inline void set_r(const long __v) {
    S_. r = __v;
  }


    

  /** returns Membrane potential in mV */
  inline double get_V_m() const {
    return S_. ode_state[State_::V_m];
  }

  inline void set_V_m(const double __v) {
    S_. ode_state[State_::V_m] = __v;
  }

    

  /**  */
  inline double get_alpha_m_init() const {
    return (0.4*(S_.ode_state[State_::V_m]/1.0+66.))/(1.-std::exp(-(S_.ode_state[State_::V_m]/1.0+66.)/5.));
  }

    

  /**  */
  inline double get_beta_m_init() const {
    return (0.4*(-(S_.ode_state[State_::V_m]/1.0+32.)))/(1.-std::exp((S_.ode_state[State_::V_m]/1.0+32.)/5.));
  }

    

  /**  */
  inline double get_Act_m() const {
    return S_. ode_state[State_::Act_m];
  }

  inline void set_Act_m(const double __v) {
    S_. ode_state[State_::Act_m] = __v;
  }

    

  /**  */
  inline double get_h_inf_init() const {
    return 1./(1.+std::exp((S_.ode_state[State_::V_m]/1.0+65.)/7.));
  }

    

  /**  */
  inline double get_Act_h() const {
    return S_. ode_state[State_::Act_h];
  }

  inline void set_Act_h(const double __v) {
    S_. ode_state[State_::Act_h] = __v;
  }

    

  /**  */
  inline double get_n_inf_init() const {
    return 1./(1.+std::exp(-(S_.ode_state[State_::V_m]/1.0+38.)/15.));
  }

    

  /**  */
  inline double get_Inact_n() const {
    return S_. ode_state[State_::Inact_n];
  }

  inline void set_Inact_n(const double __v) {
    S_. ode_state[State_::Inact_n] = __v;
  }

    

  /**  */
  inline double get_I_syn_in__1() const {
    return S_. ode_state[State_::I_syn_in__1];
  }

  inline void set_I_syn_in__1(const double __v) {
    S_. ode_state[State_::I_syn_in__1] = __v;
  }

    

  /**  */
  inline double get_I_syn_in() const {
    return S_. ode_state[State_::I_syn_in];
  }

  inline void set_I_syn_in(const double __v) {
    S_. ode_state[State_::I_syn_in] = __v;
  }

    

  /**  */
  inline double get_I_syn_ex__1() const {
    return S_. ode_state[State_::I_syn_ex__1];
  }

  inline void set_I_syn_ex__1(const double __v) {
    S_. ode_state[State_::I_syn_ex__1] = __v;
  }

    

  /**  */
  inline double get_I_syn_ex() const {
    return S_. ode_state[State_::I_syn_ex];
  }

  inline void set_I_syn_ex(const double __v) {
    S_. ode_state[State_::I_syn_ex] = __v;
  }


    

  /** returns Refractory period in ms */
  inline double get_t_ref() const {
    return P_. t_ref;
  }

  inline void set_t_ref(const double __v) {
    P_. t_ref = __v;
  }

    

  /** returns Sodium peak conductance in 1 / GOhm */
  inline double get_g_Na() const {
    return P_. g_Na;
  }

  inline void set_g_Na(const double __v) {
    P_. g_Na = __v;
  }

    

  /** returns Leak conductance in 1 / GOhm */
  inline double get_g_L() const {
    return P_. g_L;
  }

  inline void set_g_L(const double __v) {
    P_. g_L = __v;
  }

    

  /** returns Potassium peak conductance in 1 / GOhm */
  inline double get_g_K() const {
    return P_. g_K;
  }

  inline void set_g_K(const double __v) {
    P_. g_K = __v;
  }

    

  /** returns Membrane capacitance in pF */
  inline double get_C_m() const {
    return P_. C_m;
  }

  inline void set_C_m(const double __v) {
    P_. C_m = __v;
  }

    

  /**  */
  inline double get_E_Na() const {
    return P_. E_Na;
  }

  inline void set_E_Na(const double __v) {
    P_. E_Na = __v;
  }

    

  /**  */
  inline double get_E_K() const {
    return P_. E_K;
  }

  inline void set_E_K(const double __v) {
    P_. E_K = __v;
  }

    

  /**  */
  inline double get_E_L() const {
    return P_. E_L;
  }

  inline void set_E_L(const double __v) {
    P_. E_L = __v;
  }

    

  /** returns Rise time of the excitatory synaptic alpha function i in ms */
  inline double get_tau_syn_ex() const {
    return P_. tau_syn_ex;
  }

  inline void set_tau_syn_ex(const double __v) {
    P_. tau_syn_ex = __v;
  }

    

  /** returns Rise time of the inhibitory synaptic alpha function in ms */
  inline double get_tau_syn_in() const {
    return P_. tau_syn_in;
  }

  inline void set_tau_syn_in(const double __v) {
    P_. tau_syn_in = __v;
  }

    

  /** returns Constant Current in pA in pA */
  inline double get_I_e() const {
    return P_. I_e;
  }

  inline void set_I_e(const double __v) {
    P_. I_e = __v;
  }


    

  /** returns refractory time in steps in integer */
  inline long get_RefractoryCounts() const {
    return V_. RefractoryCounts;
  }

  inline void set_RefractoryCounts(const long __v) {
    V_. RefractoryCounts = __v;
  }


      

  /**  */
  inline double get_I_syn_exc() const {
    return S_. I_syn_exc;
  }

  inline void set_I_syn_exc(const double __v) {
    S_. I_syn_exc = __v;
  }

      

  /**  */
  inline double get_I_syn_inh() const {
    return S_. I_syn_inh;
  }

  inline void set_I_syn_inh(const double __v) {
    S_. I_syn_inh = __v;
  }

      

  /**  */
  inline double get_I_Na() const {
    return P_.g_Na*S_.ode_state[State_::Act_m]*S_.ode_state[State_::Act_m]*S_.ode_state[State_::Act_m]*S_.ode_state[State_::Act_h]*(S_.ode_state[State_::V_m]-P_.E_Na);
  }

      

  /**  */
  inline double get_I_K() const {
    return P_.g_K*S_.ode_state[State_::Inact_n]*S_.ode_state[State_::Inact_n]*S_.ode_state[State_::Inact_n]*S_.ode_state[State_::Inact_n]*(S_.ode_state[State_::V_m]-P_.E_K);
  }

      

  /**  */
  inline double get_I_L() const {
    return P_.g_L*(S_.ode_state[State_::V_m]-P_.E_L);
  }

      

  /**  */
  inline double get_n_inf() const {
    return 1./(1.+std::exp(-(S_.ode_state[State_::V_m]/1.0+38.)/15.));
  }

      

  /**  */
  inline double get_n_tau() const {
    return 5./(std::exp((S_.ode_state[State_::V_m]/1.0+50.)/40.)+std::exp(-(S_.ode_state[State_::V_m]/1.0+50.)/50.));
  }

      

  /**  */
  inline double get_alpha_m() const {
    return (0.4*(S_.ode_state[State_::V_m]/1.0+66.))/(1.-std::exp(-(S_.ode_state[State_::V_m]/1.0+66.)/5.));
  }

      

  /**  */
  inline double get_beta_m() const {
    return (0.4*(-(S_.ode_state[State_::V_m]/1.0+32.)))/(1.-std::exp((S_.ode_state[State_::V_m]/1.0+32.)/5.));
  }

      

  /**  */
  inline double get_h_inf() const {
    return 1./(1.+std::exp((S_.ode_state[State_::V_m]/1.0+65.)/7.));
  }

      

  /**  */
  inline double get_h_tau() const {
    return 30./(std::exp((S_.ode_state[State_::V_m]/1.0+60.)/15.)+std::exp(-(S_.ode_state[State_::V_m]/1.0+60.)/16.));
  }


    inline nest::RingBuffer& get_spikeInh() {return B_.get_spikeInh(); };
    inline nest::RingBuffer& get_spikeExc() {return B_.get_spikeExc(); };
    inline nest::RingBuffer& get_currents() {return B_.get_currents(); };

  // Generate function header
  /**
   * @defgroup iaf_psc_alpha_data
   * Instances of private data structures for the different types
   * of data pertaining to the model.
   * @note The order of definitions is important for speed.
   * @{
   */
  Parameters_ P_;  // Free parameters.
  State_      S_;  // Dynamic state.
  Variables_  V_;  // Internal Variables
  Buffers_    B_;  // Buffers.
  /** @} */

  //! Mapping of recordables names to access functions
  static nest::RecordablesMap<custom_hh> recordablesMap_;

    friend int custom_hh_dynamics( double, const double y[], double f[], void* pnode );
}; /* neuron custom_hh */

inline
nest::port custom_hh::send_test_event(nest::Node& target, nest::rport receptor_type, nest::synindex, bool)
{
  // You should usually not change the code in this function.
  // It confirms that the target of connection @c c accepts @c nest::SpikeEvent on
  // the given @c receptor_type.
  nest::SpikeEvent e;
  e.set_sender(*this);

  return target.handles_test_event(e, receptor_type);
}

inline
nest::port custom_hh::handles_test_event(nest::SpikeEvent&, nest::port receptor_type)
{
    // You should usually not change the code in this function.
    // It confirms to the connection management system that we are able
    // to handle @c SpikeEvent on port 0. You need to extend the function
    // if you want to differentiate between input ports.
    if (receptor_type != 0)
      throw nest::UnknownReceptorType(receptor_type, get_name());
    return 0;

}

inline
nest::port custom_hh::handles_test_event(nest::CurrentEvent&, nest::port receptor_type)
{
  // You should usually not change the code in this function.
  // It confirms to the connection management system that we are able
  // to handle @c CurrentEvent on port 0. You need to extend the function
  // if you want to differentiate between input ports.
  if (receptor_type != 0)
  throw nest::UnknownReceptorType(receptor_type, get_name());
  return 0;
}
inline
nest::port custom_hh::handles_test_event(nest::DataLoggingRequest& dlr,
nest::port receptor_type)
{
  // You should usually not change the code in this function.
  // It confirms to the connection management system that we are able
  // to handle @c DataLoggingRequest on port 0.
  // The function also tells the built-in UniversalDataLogger that this node
  // is recorded from and that it thus needs to collect data during simulation.
  if (receptor_type != 0)
  throw nest::UnknownReceptorType(receptor_type, get_name());

  return B_.logger_.connect_logging_device(dlr, recordablesMap_);
}

// TODO call get_status on used or internal components
inline
void custom_hh::get_status(DictionaryDatum &__d) const
{
  
  def< double >(__d, "t_ref", get_t_ref());

  
  def< double >(__d, "g_Na", get_g_Na());

  
  def< double >(__d, "g_L", get_g_L());

  
  def< double >(__d, "g_K", get_g_K());

  
  def< double >(__d, "C_m", get_C_m());

  
  def< double >(__d, "E_Na", get_E_Na());

  
  def< double >(__d, "E_K", get_E_K());

  
  def< double >(__d, "E_L", get_E_L());

  
  def< double >(__d, "tau_syn_ex", get_tau_syn_ex());

  
  def< double >(__d, "tau_syn_in", get_tau_syn_in());

  
  def< double >(__d, "I_e", get_I_e());

    
  def< long >(__d, "r", get_r());



  (*__d)[nest::names::recordables] = recordablesMap_.get_list();
    def< double >(__d, nest::names::gsl_error_tol, P_.__gsl_error_tol);
    if ( P_.__gsl_error_tol <= 0. )
    {
      throw nest::BadProperty( "The gsl_error_tol must be strictly positive." );
    }

}

inline
void custom_hh::set_status(const DictionaryDatum &__d)
{
  
  double tmp_t_ref = get_t_ref();
  updateValue<double>(__d, "t_ref", tmp_t_ref);

  
  double tmp_g_Na = get_g_Na();
  updateValue<double>(__d, "g_Na", tmp_g_Na);

  
  double tmp_g_L = get_g_L();
  updateValue<double>(__d, "g_L", tmp_g_L);

  
  double tmp_g_K = get_g_K();
  updateValue<double>(__d, "g_K", tmp_g_K);

  
  double tmp_C_m = get_C_m();
  updateValue<double>(__d, "C_m", tmp_C_m);

  
  double tmp_E_Na = get_E_Na();
  updateValue<double>(__d, "E_Na", tmp_E_Na);

  
  double tmp_E_K = get_E_K();
  updateValue<double>(__d, "E_K", tmp_E_K);

  
  double tmp_E_L = get_E_L();
  updateValue<double>(__d, "E_L", tmp_E_L);

  
  double tmp_tau_syn_ex = get_tau_syn_ex();
  updateValue<double>(__d, "tau_syn_ex", tmp_tau_syn_ex);

  
  double tmp_tau_syn_in = get_tau_syn_in();
  updateValue<double>(__d, "tau_syn_in", tmp_tau_syn_in);

  
  double tmp_I_e = get_I_e();
  updateValue<double>(__d, "I_e", tmp_I_e);


  
  long tmp_r = get_r();
  updateValue<long>(__d, "r", tmp_r);


  // We now know that (ptmp, stmp) are consistent. We do not
  // write them back to (P_, S_) before we are also sure that
  // the properties to be set in the parent class are internally
  // consistent.
  Archiving_Node::set_status(__d);

  // if we get here, temporaries contain consistent set of properties
    
  set_t_ref(tmp_t_ref);

    
  set_g_Na(tmp_g_Na);

    
  set_g_L(tmp_g_L);

    
  set_g_K(tmp_g_K);

    
  set_C_m(tmp_C_m);

    
  set_E_Na(tmp_E_Na);

    
  set_E_K(tmp_E_K);

    
  set_E_L(tmp_E_L);

    
  set_tau_syn_ex(tmp_tau_syn_ex);

    
  set_tau_syn_in(tmp_tau_syn_in);

    
  set_I_e(tmp_I_e);


    
  set_r(tmp_r);


  updateValue< double >(__d, nest::names::gsl_error_tol, P_.__gsl_error_tol);
  if ( P_.__gsl_error_tol <= 0. )
  {
    throw nest::BadProperty( "The gsl_error_tol must be strictly positive." );
  }
};

#endif /* #ifndef CUSTOM_HH */
#endif /* HAVE GSL */
