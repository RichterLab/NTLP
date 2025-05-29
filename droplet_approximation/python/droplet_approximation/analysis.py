import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

from .data import create_droplet_batch
from .models import do_inference
from .physics import dydt
from .physics import timed_solve_ivp


def analyze_model_iterative_performance( model, input_parameters = None, dt = 0.05, final_time=10.0, figure_size=None ):
    """
    Creates a figure with four plots to qualitatively assess the supplied model's
    performance on a single droplet's parameters.  For both the radius and
    temperature variables the ODE outputs are plotted against the model's estimates.
    Additionally, the relative and absolute differences of each variable are
    plotted so fine-grained differences in the solutions can be reviewed.
    
    XXX: This should return the normalized RMSE for radius and temperature.
    
    Takes 3 arguments:
    
      model            - PyTorch model to compare against the ODEs.
      input_parameters - Optional NumPy array of input parameters to evaluate performance
                         on.  If omitted, a random set of input parameters are sampled.
      figure_size      - Optional sequence, of length 2, containing the width
                         and height of the figure created.  If omitted, defaults
                         to something big enough to comfortably assess a single
                         model's outputs.
    
    Returns nothing.
    
    """

    # Default to something that fits a set of four plots (in a 2x2 grid) comfortably.
    if figure_size is None:
        figure_size = (9, 8)
    
    # Smoothly evaluate the solution for 10 seconds into the future.
    NUMBER_TIME_POINTS = int(np.floor(final_time/dt))

    # Specify the colors/styles in our plots.
    TRUTH_COLOR    = "b"
    MODEL_COLOR    = "g."
    RELATIVE_COLOR = "darkmagenta"
    ABSOLUTE_COLOR = "dodgerblue"
    
    # Sample the times in log-space so we get good coverage in both the DNS and
    # LES regions.
    t_eval = np.linspace(0, final_time, NUMBER_TIME_POINTS)
    
    # Get a random droplet if we weren't provided one.
    if input_parameters is None:
        # XXX: This is wasteful.  Just call np.random.uniform()
        (input_parameters,
         _,
         _,
         _,
         _) = create_droplet_batch( 1 )

    # Report the droplet we're evaluating in case we randomly sampled it.
    print( "Inputs: {}".format( input_parameters ) )
        
    # Get truth from the ODEs.
    y0           = (input_parameters[0], input_parameters[1])    
    solution     = solve_ivp( dydt, [0, final_time], y0, method="BDF", t_eval=t_eval, args=(input_parameters[2:],) )
    truth_output = np.vstack( (solution.y[0][:], solution.y[1][:]) ).T

    # Get the model's estimate.
    model_output = do_iterative_inference( input_parameters,
                                 t_eval,
                                 model,
                                 "cpu" )

    # Compute the normalized RMSE across the entire time scale.
    rmse_0 = np.sqrt( np.mean( (truth_output[:, 0] - model_output[:, 0])**2 ) ) / np.mean( truth_output[:, 0] ) 
    rmse_1 = np.sqrt( np.mean( (truth_output[:, 1] - model_output[:, 1])**2 ) ) / np.mean( truth_output[:, 1] ) 
    print( "NRMSE: {:g}%, {:g}%".format( rmse_0 * 100, rmse_1 * 100) )
    
    # Removed fined-grained (DNS vs. LES) analysis since dt is constant
    
    # Create our figure and embed the parameters that were evaluated.
    fig_h, ax_h = plt.subplots( 2, 2, figsize=figure_size )
    fig_h.suptitle( "Droplet Size and Temperature\nRadius={:g}, Temperature={:g}, m_s={:g}, Air Temp={:g}, RH={:g}, rhoa={:g}".format(
        input_parameters[0],
        input_parameters[1],
        input_parameters[2],
        input_parameters[3],
        input_parameters[4],
        input_parameters[5]
        ) )

    # Truth vs model predictions.
    ax_h[0][0].plot( t_eval, truth_output[:, 0], TRUTH_COLOR,
                     t_eval, model_output[:, 0], MODEL_COLOR )
    ax_h[0][1].plot( t_eval, truth_output[:, 1], TRUTH_COLOR,
                     t_eval, model_output[:, 1], MODEL_COLOR )

    # Relative difference between truth and model.
    ax_h[1][0].plot( t_eval, 
                     np.abs( truth_output[:, 0] - model_output[:, 0] ) / truth_output[:, 0] * 100,
                     color=RELATIVE_COLOR )
    ax_h[1][0].tick_params( axis="y", labelcolor=RELATIVE_COLOR )
    ax_h[1][1].plot( t_eval, 
                     np.abs( truth_output[:, 1] - model_output[:, 1] ) / truth_output[:, 1] * 100,
                     color=RELATIVE_COLOR )
    ax_h[1][1].tick_params( axis="y", labelcolor=RELATIVE_COLOR )

    # Aboslute difference between truth and model.
    ax_h_twin_radius = ax_h[1][0].twinx()
    ax_h_twin_radius.plot( t_eval, np.abs( truth_output[:, 0] - model_output[:, 0] ), color=ABSOLUTE_COLOR )
    ax_h_twin_radius.set_ylabel( "Absolute Difference (m)" )
    ax_h_twin_radius.tick_params( axis="y", labelcolor=ABSOLUTE_COLOR )

    ax_h_twin_temperature = ax_h[1][1].twinx()
    ax_h_twin_temperature.plot( t_eval, np.abs( truth_output[:, 1] - model_output[:, 1] ), color=ABSOLUTE_COLOR )
    ax_h_twin_temperature.set_ylabel( "Absolute Difference (K)" )
    ax_h_twin_temperature.tick_params( axis="y", labelcolor=ABSOLUTE_COLOR )

    # Label the comparison plots' lines.
    ax_h[0][0].legend( ["Truth", "Model"] )
    ax_h[0][1].legend( ["Truth", "Model"] )

    # Label the columns of plots.
    ax_h[0][0].set_title( "Radius" )
    ax_h[0][1].set_title( "Temperature" )

    ax_h[0][0].set_ylabel( "Radius (m)" )
    ax_h[0][1].set_ylabel( "Temperature (K)" )
    ax_h[1][0].set_ylabel( "Relative Difference (%)" )
    ax_h[1][0].set_xlabel( "Time (s)" )
    ax_h[1][1].set_ylabel( "Relative Difference (%)" )
    ax_h[1][1].set_xlabel( "Time (s)" )

    # Show time in log-scale as well as the droplets' radius.
    ax_h[0][0].set_xscale( "log" )
    ax_h[0][0].set_yscale( "log" )
    ax_h[0][1].set_xscale( "log" )
    ax_h[0][1].set_yscale( "log" )
    ax_h[1][0].set_xscale( "log" )
    ax_h[1][1].set_xscale( "log" )

    fig_h.tight_layout()




def plot_droplet_size_temperature( size_temperatures, times ):
    """
    Plots a droplet's size and temperature as a function of time.

    Takes 2 arguments:

      size_temperatures - NumPy array, sized 2 x number_times, containing
                          droplet radii and temperatures in the first and second
                          columns, respectively.
      times             - NumPy vector, of length number_times, containing the
                          times corresponding to each entry in size_temperatures.

    Returns 2 values:

      fig_h - Figure handle created.
      ax_h  - Sequence of two axes handles, one for the size plot and another
              for the temperature plot.

    """

    fig_h, ax_h = plt.subplots( 1, 2, figsize=(9, 4) )

    fig_h.suptitle( "Droplet Size and Temperature (Truth)" )

    ax_h[0].plot( times, size_temperatures[0, :], label="radius" )
    ax_h[0].set_xlabel( "Time (s)" )
    ax_h[0].set_ylabel( "Radius (m)" )
    ax_h[0].set_xscale( "log" )
    ax_h[0].set_yscale( "log" )

    ax_h[1].plot( times, size_temperatures[1, :], label="temp" )
    ax_h[1].set_xlabel( "Time (s)" )
    ax_h[1].set_ylabel( "Temperature (K)" )
    ax_h[1].set_xscale( "log" )

    return fig_h, ax_h

def analyze_model_performance( model, input_parameters=None, figure_size=None ):
    """
    Creates a figure with four plots to qualitatively assess the supplied model's
    performance on a single droplet's parameters.  For both the radius and
    temperature variables the ODE outputs are plotted against the model's estimates.
    Additionally, the relative and absolute differences of each variable are
    plotted so fine-grained differences in the solutions can be reviewed.

    XXX: This should return the normalized RMSE for radius and temperature.

    Takes 3 arguments:

      model            - PyTorch model to compare against the ODEs.
      input_parameters - Optional NumPy array of input parameters to evaluate performance
                         on.  If omitted, a random set of input parameters are sampled.
      figure_size      - Optional sequence, of length 2, containing the width
                         and height of the figure created.  If omitted, defaults
                         to something big enough to comfortably assess a single
                         model's outputs.

    Returns nothing.

    """

    # Default to something that fits a set of four plots (in a 2x2 grid) comfortably.
    if figure_size is None:
        figure_size = (9, 8)

    # Smoothly evaluate the solution for 10 seconds into the future.
    FINAL_TIME         = 10.0
    NUMBER_TIME_POINTS = 1000

    # Specify the colors/styles in our plots.
    TRUTH_COLOR    = "b"
    MODEL_COLOR    = "r."
    RELATIVE_COLOR = "darkmagenta"
    ABSOLUTE_COLOR = "dodgerblue"

    # Sample the times in log-space so we get good coverage in both the DNS and
    # LES regions.
    t_eval = np.logspace( -3, np.log10( FINAL_TIME ), NUMBER_TIME_POINTS )

    # Get a random droplet if we weren't provided one.
    if input_parameters is None:
        # XXX: This is wasteful.  Just call np.random.uniform()
        (input_parameters,
         _,
         _,
         _,
         _) = create_droplet_batch( 1 )

    # Report the droplet we're evaluating in case we randomly sampled it.
    print( "Inputs: {}".format( input_parameters ) )

    # Get truth from the ODEs.
    y0           = (input_parameters[0, 0], input_parameters[0, 1])
    solution     = solve_ivp( dydt, [0, FINAL_TIME], y0, method="Radau", t_eval=t_eval, args=(input_parameters[0, 2:],) )
    truth_output = np.vstack( (solution.y[0][:], solution.y[1][:]) ).T

    # Get the model's estimate.
    model_output = do_inference( np.tile( input_parameters, (t_eval.shape[0], 1) ),
                                 t_eval,
                                 model,
                                 "cpu" )

    # Compute the normalized RMSE across the entire time scale.
    rmse_0 = np.sqrt( np.mean( (truth_output[:, 0] - model_output[:, 0])**2 ) ) / np.mean( truth_output[:, 0] )
    rmse_1 = np.sqrt( np.mean( (truth_output[:, 1] - model_output[:, 1])**2 ) ) / np.mean( truth_output[:, 1] )
    print( "NRMSE: {:g}%, {:g}%".format( rmse_0 * 100, rmse_1 * 100) )

    # Compute the normalized RMSE for different regions (DNS and LES) for
    # finer-grained performance assessment.
    t_eval_mask = np.empty( (NUMBER_TIME_POINTS, 4),
                            dtype=np.bool_ )
    t_eval_mask[:, 0] = (t_eval  < 1e-2)
    t_eval_mask[:, 1] = (t_eval >= 1e-2) & (t_eval < 1e-1)
    t_eval_mask[:, 2] = (t_eval >= 1e-1) & (t_eval < 1e0)
    t_eval_mask[:, 3] = (t_eval >= 1e0)

    rmse_0 = np.empty( (4,), dtype=np.float32 )
    rmse_1 = np.empty( (4,), dtype=np.float32 )

    for scale_index in np.arange( t_eval_mask.shape[1] ):
        rmse_0[scale_index] = np.sqrt( np.mean( (truth_output[t_eval_mask[:, scale_index], 0] - model_output[t_eval_mask[:, scale_index], 0])**2 ) ) / np.mean( truth_output[t_eval_mask[:, scale_index], 0] )
        rmse_1[scale_index] = np.sqrt( np.mean( (truth_output[t_eval_mask[:, scale_index], 1] - model_output[t_eval_mask[:, scale_index], 1])**2 ) ) / np.mean( truth_output[t_eval_mask[:, scale_index], 1] )

    print( "NRMSE: {}%, {}%".format( rmse_0 * 100, rmse_1 * 100 ) )

    # Create our figure and embed the parameters that were evaluated.
    fig_h, ax_h = plt.subplots( 2, 2, figsize=figure_size )
    fig_h.suptitle( "Droplet Size and Temperature\nRadius={:g}, Temperature={:g}, m_s={:g}, Air Temp={:g}, RH={:g}, rhoa={:g}".format(
        input_parameters[0, 0],
        input_parameters[0, 1],
        input_parameters[0, 2],
        input_parameters[0, 3],
        input_parameters[0, 4],
        input_parameters[0, 5]
        ) )

    # Truth vs model predictions.
    ax_h[0][0].plot( t_eval, truth_output[:, 0], TRUTH_COLOR,
                     t_eval, model_output[:, 0], MODEL_COLOR )
    ax_h[0][1].plot( t_eval, truth_output[:, 1], TRUTH_COLOR,
                     t_eval, model_output[:, 1], MODEL_COLOR )

    # Relative difference between truth and model.
    ax_h[1][0].plot( t_eval,
                     np.abs( truth_output[:, 0] - model_output[:, 0] ) / truth_output[:, 0] * 100,
                     color=RELATIVE_COLOR )
    ax_h[1][0].tick_params( axis="y", labelcolor=RELATIVE_COLOR )
    ax_h[1][1].plot( t_eval,
                     np.abs( truth_output[:, 1] - model_output[:, 1] ) / truth_output[:, 1] * 100,
                     color=RELATIVE_COLOR )
    ax_h[1][1].tick_params( axis="y", labelcolor=RELATIVE_COLOR )

    # Aboslute difference between truth and model.
    ax_h_twin_radius = ax_h[1][0].twinx()
    ax_h_twin_radius.plot( t_eval, np.abs( truth_output[:, 0] - model_output[:, 0] ), color=ABSOLUTE_COLOR )
    ax_h_twin_radius.set_ylabel( "Absolute Difference (m)" )
    ax_h_twin_radius.tick_params( axis="y", labelcolor=ABSOLUTE_COLOR )

    ax_h_twin_temperature = ax_h[1][1].twinx()
    ax_h_twin_temperature.plot( t_eval, np.abs( truth_output[:, 1] - model_output[:, 1] ), color=ABSOLUTE_COLOR )
    ax_h_twin_temperature.set_ylabel( "Absolute Difference (K)" )
    ax_h_twin_temperature.tick_params( axis="y", labelcolor=ABSOLUTE_COLOR )

    # Label the comparison plots' lines.
    ax_h[0][0].legend( ["Truth", "Model"] )
    ax_h[0][1].legend( ["Truth", "Model"] )

    # Label the columns of plots.
    ax_h[0][0].set_title( "Radius" )
    ax_h[0][1].set_title( "Temperature" )

    ax_h[0][0].set_ylabel( "Radius (m)" )
    ax_h[0][1].set_ylabel( "Temperature (K)" )
    ax_h[1][0].set_ylabel( "Relative Difference (%)" )
    ax_h[1][0].set_xlabel( "Time (s)" )
    ax_h[1][1].set_ylabel( "Relative Difference (%)" )
    ax_h[1][1].set_xlabel( "Time (s)" )

    # Show time in log-scale as well as the droplets' radius.
    ax_h[0][0].set_xscale( "log" )
    ax_h[0][0].set_yscale( "log" )
    ax_h[0][1].set_xscale( "log" )
    ax_h[0][1].set_yscale( "log" )
    ax_h[1][0].set_xscale( "log" )
    ax_h[1][1].set_xscale( "log" )

    fig_h.tight_layout()
