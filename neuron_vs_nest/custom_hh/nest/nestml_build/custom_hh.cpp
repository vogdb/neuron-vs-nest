/*
*  custom_hh.cpp
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

// C++ includes:
#include <limits>

// Includes from libnestutil:
#include "numerics.h"

// Includes from nestkernel:
#include "exceptions.h"
#include "kernel_manager.h"
#include "universal_data_logger_impl.h"

// Includes from sli:
#include "dict.h"
#include "dictutils.h"
#include "doubledatum.h"
#include "integerdatum.h"
#include "lockptrdatum.h"

#include "custom_hh.h"

/* ----------------------------------------------------------------
* Recordables map
* ---------------------------------------------------------------- */
nest::RecordablesMap<custom_hh> custom_hh::recordablesMap_;

namespace nest
{
  // Override the create() method with one call to RecordablesMap::insert_()
  // for each quantity to be recorded.
  template <>
  void RecordablesMap<custom_hh>::create()
  {
    // use standard names whereever you can for consistency!
      



      



      



      



      



      



      



      



      



      



      



      



      



      


insert_("V_m", &custom_hh::get_V_m);

      


insert_("alpha_m_init", &custom_hh::get_alpha_m_init);

      


insert_("beta_m_init", &custom_hh::get_beta_m_init);

      


insert_("Act_m", &custom_hh::get_Act_m);

      


insert_("h_inf_init", &custom_hh::get_h_inf_init);

      


insert_("Act_h", &custom_hh::get_Act_h);

      


insert_("n_inf_init", &custom_hh::get_n_inf_init);

      


insert_("Inact_n", &custom_hh::get_Inact_n);

      


insert_("I_syn_in__1", &custom_hh::get_I_syn_in__1);

      


insert_("I_syn_in", &custom_hh::get_I_syn_in);

      


insert_("I_syn_ex__1", &custom_hh::get_I_syn_ex__1);

      


insert_("I_syn_ex", &custom_hh::get_I_syn_ex);

      



      



      



      



      



      



      



      



      



      



      



  }

}

/* ----------------------------------------------------------------
 * Default constructors defining default parameters and state
 * Note: the implementation is empty. The initialization is of variables
 * is a part of the custom_hh's constructor.
 * ---------------------------------------------------------------- */
custom_hh::Parameters_::Parameters_()
{
}

custom_hh::State_::State_()
{
}

/* ----------------------------------------------------------------
* Parameter and state extractions and manipulation functions
* ---------------------------------------------------------------- */

custom_hh::Buffers_::Buffers_(custom_hh &n): logger_(n)
  , __s( 0 )
  , __c( 0 )
  , __e( 0 )
{
  // Initialization of the remaining members is deferred to
  // init_buffers_().
}

custom_hh::Buffers_::Buffers_(const Buffers_ &, custom_hh &n): logger_(n)
  , __s( 0 )
  , __c( 0 )
  , __e( 0 )
{
  // Initialization of the remaining members is deferred to
  // init_buffers_().
}

/* ----------------------------------------------------------------
 * Default and copy constructor for node, and destructor
 * ---------------------------------------------------------------- */
custom_hh::custom_hh():Archiving_Node(), P_(), S_(), B_(*this)
{
  recordablesMap_.create();
     // use a default `good` enough value for the absolute error.
     // it cab be adjusted via `SetStatus`
     P_.__gsl_error_tol = 1e-3;

    
    P_.t_ref = 2.0*1.0; // as ms

    
    P_.g_Na = 5000.0*1.0; // as 1 / GOhm

    
    P_.g_L = 200.0*1.0; // as 1 / GOhm

    
    P_.g_K_rect = 30000.0*1.0; // as 1 / GOhm

    
    P_.C_m = 100.0*1.0; // as pF

    
    P_.E_Na = 50.0*1.0; // as mV

    
    P_.E_K = -80.0*1.0; // as mV

    
    P_.E_L = -70.0*1.0; // as mV

    
    P_.tau_syn_ex = 0.2*1.0; // as ms

    
    P_.tau_syn_in = 2.0*1.0; // as ms

    
    P_.I_e = 0*1.0; // as pA


    
    S_.r = 0; // as integer


    
    S_.ode_state[State_::V_m] = -70.*1.0; // as mV

    
    S_.alpha_m_init = (0.4*(S_.ode_state[State_::V_m]/1.0+66.))/(1.-std::exp(-(S_.ode_state[State_::V_m]/1.0+66.)/5.)); // as real

    
    S_.beta_m_init = (0.4*(-(S_.ode_state[State_::V_m]/1.0+32.)))/(1.-std::exp((S_.ode_state[State_::V_m]/1.0+32.)/5.)); // as real

    
    S_.ode_state[State_::Act_m] = get_alpha_m_init()/(get_alpha_m_init()+get_beta_m_init()); // as real

    
    S_.h_inf_init = 1./(1.+std::exp((S_.ode_state[State_::V_m]/1.0+65.)/7.)); // as real

    
    S_.ode_state[State_::Act_h] = get_h_inf_init(); // as real

    
    S_.n_inf_init = 1./(1.+std::exp(-(S_.ode_state[State_::V_m]/1.0+38.)/15.)); // as real

    
    S_.ode_state[State_::Inact_n] = get_n_inf_init(); // as real

    
    S_.ode_state[State_::I_syn_in__1] = 0; // as real

    
    S_.ode_state[State_::I_syn_in] = 0; // as real

    
    S_.ode_state[State_::I_syn_ex__1] = 0; // as real

    
    S_.ode_state[State_::I_syn_ex] = 0; // as real



}

custom_hh::custom_hh(const custom_hh& __n): Archiving_Node(), P_(__n.P_), S_(__n.S_), B_(__n.B_, *this)
{
    P_.t_ref = __n.P_.t_ref;
    P_.g_Na = __n.P_.g_Na;
    P_.g_L = __n.P_.g_L;
    P_.g_K_rect = __n.P_.g_K_rect;
    P_.C_m = __n.P_.C_m;
    P_.E_Na = __n.P_.E_Na;
    P_.E_K = __n.P_.E_K;
    P_.E_L = __n.P_.E_L;
    P_.tau_syn_ex = __n.P_.tau_syn_ex;
    P_.tau_syn_in = __n.P_.tau_syn_in;
    P_.I_e = __n.P_.I_e;

    S_.r = __n.S_.r;

    V_.RefractoryCounts = __n.V_.RefractoryCounts;
}

custom_hh::~custom_hh()
{
    // GSL structs may not have been allocated, so we need to protect destruction
    if ( B_.__s )
      gsl_odeiv_step_free( B_.__s );
    if ( B_.__c )
      gsl_odeiv_control_free( B_.__c );
    if ( B_.__e )
      gsl_odeiv_evolve_free( B_.__e );
}

/* ----------------------------------------------------------------
* Node initialization functions
* ---------------------------------------------------------------- */

void
custom_hh::init_state_(const Node& proto)
{
  const custom_hh& pr = downcast<custom_hh>(proto);
  S_ = pr.S_;
}

extern "C" inline int
custom_hh_dynamics( double, const double ode_state[], double f[], void* pnode )
{
  typedef custom_hh::State_ State_;
  // get access to node so we can almost work as in a member function
  assert( pnode );
  const custom_hh& node = *( reinterpret_cast< custom_hh* >( pnode ) );

  // ode_state[] here is---and must be---the state vector supplied by the integrator,
  // not the state vector in the node, node.S_.ode_state[].


    double I_syn_exc = ode_state[State_::I_syn_ex];
    double I_syn_inh = ode_state[State_::I_syn_in];
    double I_Na = node.get_g_Na()*ode_state[State_::Act_m]*ode_state[State_::Act_m]*ode_state[State_::Act_m]*ode_state[State_::Act_h]*(ode_state[State_::V_m]-node.get_E_Na());
    double I_K = node.get_g_K_rect()*ode_state[State_::Inact_n]*ode_state[State_::Inact_n]*ode_state[State_::Inact_n]*ode_state[State_::Inact_n]*(ode_state[State_::V_m]-node.get_E_K());
    double I_L = node.get_g_L()*(ode_state[State_::V_m]-node.get_E_L());
    double n_inf = 1./(1.+std::exp(-(ode_state[State_::V_m]/1.0+38.)/15.));
    double n_tau = 5./(std::exp((ode_state[State_::V_m]/1.0+50.)/40.)+std::exp(-(ode_state[State_::V_m]/1.0+50.)/50.));
    double alpha_m = (0.4*(ode_state[State_::V_m]/1.0+66.))/(1.-std::exp(-(ode_state[State_::V_m]/1.0+66.)/5.));
    double beta_m = (0.4*(-(ode_state[State_::V_m]/1.0+32.)))/(1.-std::exp((ode_state[State_::V_m]/1.0+32.)/5.));
    double h_inf = 1./(1.+std::exp((ode_state[State_::V_m]/1.0+65.)/7.));
    double h_tau = 30./(std::exp((ode_state[State_::V_m]/1.0+60.)/15.)+std::exp(-(ode_state[State_::V_m]/1.0+60.)/16.));

    f[ State_::V_m ] = (-(I_Na+I_K+I_L)+node.B_.currents_grid_sum_+node.get_I_e()+I_syn_inh+I_syn_exc)/node.get_C_m();
    f[ State_::Act_m ] = (alpha_m*(1.-ode_state[State_::Act_m])-beta_m*ode_state[State_::Act_m])/1.0;
    f[ State_::Act_h ] = (h_inf-ode_state[State_::Act_h])/h_tau;
    f[ State_::Inact_n ] = (n_inf-ode_state[State_::Inact_n])/n_tau;
    f[ State_::I_syn_in__1 ] = -1/pow((node.get_tau_syn_in()), (2))*ode_state[State_::I_syn_in]+-2/node.get_tau_syn_in()*ode_state[State_::I_syn_in__1];
    f[ State_::I_syn_in ] = ode_state[State_::I_syn_in__1];
    f[ State_::I_syn_ex__1 ] = -1/pow((node.get_tau_syn_ex()), (2))*ode_state[State_::I_syn_ex]+-2/node.get_tau_syn_ex()*ode_state[State_::I_syn_ex__1];
    f[ State_::I_syn_ex ] = ode_state[State_::I_syn_ex__1];

  return GSL_SUCCESS;
}




void
custom_hh::init_buffers_()
{
  get_spikeInh().clear(); //includes resize
  get_spikeExc().clear(); //includes resize
  get_currents().clear(); //includes resize
  B_.logger_.reset(); // includes resize
  Archiving_Node::clear_history();
    if ( B_.__s == 0 )
    {
      B_.__s = gsl_odeiv_step_alloc( gsl_odeiv_step_rkf45, 8 );
    }
    else
    {
      gsl_odeiv_step_reset( B_.__s );
    }

    if ( B_.__c == 0 )
    {
      B_.__c = gsl_odeiv_control_y_new( P_.__gsl_error_tol, 0.0 );
    }
    else
    {
      gsl_odeiv_control_init( B_.__c, P_.__gsl_error_tol, 0.0, 1.0, 0.0 );
    }

    if ( B_.__e == 0 )
    {
      B_.__e = gsl_odeiv_evolve_alloc( 8 );
    }
    else
    {
      gsl_odeiv_evolve_reset( B_.__e );
    }

    B_.__sys.function = custom_hh_dynamics;
    B_.__sys.jacobian = NULL;
    B_.__sys.dimension = 8;
    B_.__sys.params = reinterpret_cast< void* >( this );
    B_.__step = nest::Time::get_resolution().get_ms();
    B_.__integration_step = nest::Time::get_resolution().get_ms();

}

void
custom_hh::calibrate()
{
  B_.logger_.init();

    

V_. RefractoryCounts =
    nest::Time(nest::Time::ms((double) P_.t_ref)).get_steps()
;






}

/* ----------------------------------------------------------------
* Update and spike handling functions
* ---------------------------------------------------------------- */

/*
 
 */
void
custom_hh::update(
        nest::Time const & origin,
        const long from, const long to)
{
    double __t = 0;

  for ( long lag = from ; lag < to ; ++lag ) {
          B_.spikeInh_grid_sum_ = get_spikeInh().get_value( lag );
          B_.spikeExc_grid_sum_ = get_spikeExc().get_value( lag );
          B_.currents_grid_sum_ = get_currents().get_value( lag );

          
      double U_old = S_.ode_state[State_::V_m];



      __t = 0;
// numerical integration with adaptive step size control:
// ------------------------------------------------------
// gsl_odeiv_evolve_apply performs only a single numerical
// integration step, starting from t and bounded by step;
// the while-loop ensures integration over the whole simulation
// step (0, step] if more than one integration step is needed due
// to a small integration step size;
// note that (t+IntegrationStep > step) leads to integration over
// (t, step] and afterwards setting t to step, but it does not
// enforce setting IntegrationStep to step-t; this is of advantage
// for a consistent and efficient integration across subsequent
// simulation intervals
while ( __t < B_.__step )
{
  const int status = gsl_odeiv_evolve_apply(B_.__e,
                                            B_.__c,
                                            B_.__s,
                                            &B_.__sys,              // system of ODE
                                            &__t,                   // from t
                                            B_.__step,              // to t <= step
                                            &B_.__integration_step, // integration step size
                                            S_.ode_state);          // neuronal state

  if ( status != GSL_SUCCESS ) {
    throw nest::GSLSolverFailure( get_name(), status );
  }
}



  // # sending spikes: crossing 0 mV, pseudo-refractoriness and local maximum...
    if (S_.r>0) {
        S_.r
  -=
  1;




} else if((((S_.ode_state[State_::V_m]>0*1.0))) && (((U_old>S_.ode_state[State_::V_m])))) {
        S_.r
  =
  V_.RefractoryCounts;



      set_spiketime(nest::Time::step(origin.get_steps()+lag+1));
nest::SpikeEvent se;
nest::kernel().event_delivery_manager.send(*this, se, lag);;




} /* if end */




        S_.ode_state[State_::I_syn_ex__1]
  +=
  B_.spikeExc_grid_sum_*numerics::e/P_.tau_syn_ex;



        S_.ode_state[State_::I_syn_ex]
  +=
  B_.spikeExc_grid_sum_*0;



        S_.ode_state[State_::I_syn_in__1]
  +=
  B_.spikeInh_grid_sum_*numerics::e/P_.tau_syn_in;



        S_.ode_state[State_::I_syn_in]
  +=
  B_.spikeInh_grid_sum_*0;





    // voltage logging
    B_.logger_.record_data(origin.get_steps()+lag);
  }

}

// Do not move this function as inline to h-file. It depends on
// universal_data_logger_impl.h being included here.
void
custom_hh::handle(nest::DataLoggingRequest& e)
{
  B_.logger_.handle(e);
}


void
custom_hh::handle(nest::SpikeEvent &e)
{
  assert(e.get_delay() > 0);

      const double weight = e.get_weight();
      const double multiplicity = e.get_multiplicity();
        if ( weight < 0.0 ) // inhibitory
        {
          get_spikeInh().add_value(e.get_rel_delivery_steps( nest::kernel().simulation_manager.get_slice_origin()),
                      
                       weight * multiplicity );
        }
        if ( weight >= 0.0 ) // excitatory
        {
          get_spikeExc().add_value(e.get_rel_delivery_steps( nest::kernel().simulation_manager.get_slice_origin()),
                       weight * multiplicity );
        }

}

void
custom_hh::handle(nest::CurrentEvent& e)
{
  assert(e.get_delay() > 0);

  const double current=e.get_current();
  const double weight=e.get_weight();

  // add weighted current; HEP 2002-10-04
    get_currents().add_value(
               e.get_rel_delivery_steps( nest::kernel().simulation_manager.get_slice_origin()),
               weight * current );
}
