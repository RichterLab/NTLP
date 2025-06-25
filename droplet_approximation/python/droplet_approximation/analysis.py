import multiprocessing

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

from .data import create_droplet_batch, read_training_file
from .models import do_inference, do_iterative_inference
from .physics import dydt, scale_droplet_parameters
from .physics import DROPLET_TIME_LOG_RANGE, normalize_droplet_parameters

def analyze_model_iterative_performance( model, input_parameters=None, dt=0.05, final_time=10.0, figure_size=None ):
    """
    Creates a figure with four plots to qualitatively assess the supplied model's
    iterative performance on a single droplet's parameters. For both the radius and
    temperature variables the ODE outputs are plotted against the model's estimates.
    Additionally, the relative and absolute differences of each variable are
    plotted so fine-grained differences in the solutions can be reviewed.

    The model output is calculated iteratively. For instance, if `dt=0.5`, the model
    output radius/temperature from the input parameters for `dt=0.5` will be used as
    the input to calculate the radius/temperature for `t=1` and so on.

    Takes 5 arguments:

      model            - PyTorch model to compare against the ODEs.
      input_parameters - Optional NumPy array of input parameters to evaluate performance
                         on.  If omitted, a random set of input parameters are sampled.
                         The last four parameters will be fixed throughout the integration
                         while the first two will be replace with the model output
                         at each time step.
      dt               - Optional float specifying the timestep to evaluate successive
                         particles at.  If omitted, defaults to 0.05 seconds.
      final_time       - Optional float fixing how long to integrate out to in seconds.
                         If omitted, defaults to 10 seconds.
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
    NUMBER_TIME_POINTS = int( np.floor( final_time / dt ) )

    # Specify the colors/styles in our plots.
    TRUTH_COLOR    = "b"
    MODEL_COLOR    = "g."
    RELATIVE_COLOR = "darkmagenta"
    ABSOLUTE_COLOR = "dodgerblue"

    # Sample the times in log-space so we get good coverage in both the DNS and
    # LES regions.
    t_eval = np.linspace( 0, final_time, NUMBER_TIME_POINTS )

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
        rmse_0[scale_index] = (np.sqrt( np.mean( (truth_output[t_eval_mask[:, scale_index], 0] -
                                                 model_output[t_eval_mask[:, scale_index], 0])**2 ) ) /
                               np.mean( truth_output[t_eval_mask[:, scale_index], 0] ))
        rmse_1[scale_index] = (np.sqrt( np.mean( (truth_output[t_eval_mask[:, scale_index], 1] -
                                                  model_output[t_eval_mask[:, scale_index], 1])**2 ) ) /
                               np.mean( truth_output[t_eval_mask[:, scale_index], 1] ))

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

def analyze_model_particle_performance( times, truth_output, model_output, distances, title_string=None, figure_size=None, time_range=None ):
    """
    Creates a figure with four plots to qualitatively assess the supplied model's
    performance on a particular particle's trajectory from an NTLP dump. For both
    the radius and temperature variables the ODE outputs are plotted against the
    model's estimates.

    The NRMSE is calculated to assess overall model performance. Additionally, the
    relative and absolute differences of each variable are plotted so fine-grained
    differences in the solutions can be reviewed.

    Requires a dataframe with both the mlp output and distance already calculated.

    Takes 7 arguments:

      times           - NumPy array, contains an array of the simulation
                        times being analyzed
      truth_ouput     - NumPy array, shaped 2 x len( times ). Contains an
                        array of the particle's actual radius and temperature
                        at index 0 and 1 respectively.
      model_output    - NumPy array, shaped 2 x len( times ). Contains an
                        array of the particle's estimated radius and temperature
                        at index 0 and 1 respectively.
      distances       - NumPy array, shaped 2 x len( times ). Contains
                        an array of the distance between the particle's true
                        output and the MLP output.
      title_string    - Optional string to append to plot title
      figure_size     - Optional sequence, of length 2, containing the width
                        and height of the figure created.  If omitted, defaults
                        to something big enough to comfortably assess a single
                        model's outputs.
      time_range      - Optional sequence, of length 2, specifying the start and
                        end time to display performance for.  If omitted, defaults
                        to None and the times input is used for the X axis.

    Returns nothing.

    """

    # Default to something that fits a set of four plots (in a 2x2 grid) comfortably.
    if figure_size is None:
        figure_size = (9, 8)

    # Specify the colors/styles in our plots.
    TRUTH_COLOR    = "b"
    MODEL_COLOR    = "r."
    RELATIVE_COLOR = "darkmagenta"
    ABSOLUTE_COLOR = "dodgerblue"

    # Sort the dataframe by time
    t_eval = times

    # Compute the normalized RMSE across the entire time scale.
    rmse_0 = np.sqrt( np.mean( distances[:, 0]**2 ) ) / np.abs( np.mean( np.log10( truth_output[:, 0] ) ) ) # TODO fix hard-coded log
    rmse_1 = np.sqrt( np.mean( distances[:, 1]**2 ) ) / np.abs( np.mean( truth_output[:, 1] ) )

    # Create our figure and embed the parameters that were evaluated.
    fig_h, ax_h = plt.subplots( 2, 2, figsize=figure_size )
    fig_h.suptitle( "Droplet Size and Temperature\nRadius NRMSE={:g}%, Temperature NRMSE={:g}%\n{}".format( rmse_0 * 100, rmse_1 * 100, title_string) )

    # Truth vs model predictions.
    ax_h[0][0].plot( t_eval, truth_output[:, 0], TRUTH_COLOR,
                     t_eval, model_output[:, 0], MODEL_COLOR,
                     ms=2 )
    ax_h[0][1].plot( t_eval, truth_output[:, 1], TRUTH_COLOR,
                     t_eval, model_output[:, 1], MODEL_COLOR,
                     ms=2 )

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

    if time_range is not None:
        ax_h[0][0].set_xlim( time_range )
        ax_h[0][1].set_xlim( time_range )
        ax_h[1][0].set_xlim( time_range )
        ax_h[1][1].set_xlim( time_range )

    fig_h.tight_layout()

def calculate_nrmse( truth_output, model_output, distance_function ):
    """
    Calculates the normalized root-mean squared error (NRMSE) between the
    provided truth and model outputs using a user-supplied error function.
    Returns the NRMSE across all observations (rows).

    Takes 3 arguments:

      truth_output      - NumPy array, sized number_observations x 2, containing
                          the truth values for radii and temperatures.
      model_output      - NumPy array, sized number_observations x 2, containing
                          the model values for radii and temparatures.
      distance_function - The distance function to compute the error between
                          truth and model outputs.  Takes two arguments,
                          truth_output and model_output, and returns a NumPy
                          array, sized number_observations x 2, containing the
                          error (distance) between the inputs.

    Returns 2 values:

      radii_nrmse        - NRMSE of the radii.
      temperatures_nrmse - NRMSE of the temperatures.

    """

    distances = distance_function( truth_output, model_output )

    return (np.sqrt( np.mean( distances[:, 0]**2 ) ) /
            np.abs( np.mean( np.log10( truth_output[:, 0] ) ) ) +
            np.sqrt( np.mean( distances[:, 1]**2 ) ) /
            np.abs( np.mean( truth_output[:, 1] ) ))

def mse_score_models( models, file_name, device, weighted=False, normalized=False ):
    """
    DEPRECATED.

    Will be rewritten after merge with new methods written locally.
    Calculates the mean square error on an array of models for a dataset.

    Takes 5 arguments:

      models          - Sequence of PyTorch models to be evaluated.
      file_name       - Path to training data.
      device          - Device string to perform model evaluation on.
      weighted        - Optional boolean, defaults to False, if True,
                        weights MSE loss by the reciprocal of the
                        integration time.
      normalized      - Optional boolean, defaults to False, if True,
                        normalizes droplet parameters before
                        calculating MSE loss.

    Returns 1 value:

      losses          - NumPy array, length of models, MSE loss for
                        each model on the provided data set.

    """

    input_parameters, output_parameters, integration_times = read_training_file( file_name )

    # Normalizes the radii' reciprocal logarithmically since it is harder
    # to learn.
    weights = 10**((DROPLET_TIME_LOG_RANGE[0] + DROPLET_TIME_LOG_RANGE[1]) / 2.0) * np.reciprocal( integration_times )
    weights = np.stack( (weights, weights), axis=-1 )

    BATCH_SIZE = 1024 * 10

    number_droplets = len( input_parameters )
    number_batches  = (number_droplets + BATCH_SIZE + 1) // BATCH_SIZE

    losses = np.zeros( shape=(len( models ), 2), dtype=np.float64 )
    for batch_index in range( number_batches ):
        if batch_index != (number_batches - 1):
            batch_size = BATCH_SIZE
        else:
            batch_size = number_droplets % BATCH_SIZE

        start_index = BATCH_SIZE * batch_index
        end_index   = start_index + batch_size

        inputs          = input_parameters[start_index:end_index]
        target_outputs  = output_parameters[start_index:end_index]
        times           = integration_times[start_index:end_index]
        current_weights = weights[start_index:end_index]

        # If this gets too large could average over batches, however
        # not all batches are the same size. Could multiply by
        # batch_size/BATCH_SIZE and then average by batch size if
        # need be
        for model_index in range( len( models ) ):
            inferred_outputs = do_inference( inputs,
                                             times,
                                             models[model_index],
                                             device )
            if normalized:
                error = (normalize_droplet_parameters( inferred_outputs ) -
                         normalize_droplet_parameters( target_outputs ))
            else:
                error = np.array( [np.log10( inferred_outputs[:, 0] /
                                             target_outputs[:, 0] ),
                                   inferred_outputs[:, 1] - target_outputs[:, 1]] )
            if weighted:
                losses[model_index] += ((current_weights * error)**2).sum( axis=1 )
            else:
                losses[model_index] += (error**2).sum( axis=1 )

    losses /= number_droplets

    return losses

def parallel_analyze_model_iterative_performance_NTLP_data( model, df, iterations, device, cores=64 ):
    """
    NEEDS TO BE REWORKED - DO NOT USE

    Calculates the model output by integrating along with each particle's
    background conditions for `iterations` time steps. Groups the dataset
    by processor that spawned each particle and then pools these jobs
    onto `cores` workers.

    Takes 5 arguments:

      model            - PyTorch model to compare against the ODEs.
      df               - Pandas DataFrame to evaluate particles from.
      iterations       - Integer, number of time steps to iterate through
                         for each input.
      device           - Device to evaluate model on.
      cores            - Integer, number of workers to parallelize across.  If omitted
                         defaults to 64 cores.

    Returns 1 value:

      output           - NumPy array of output radii and temperatures
                         after `iterations` time steps, shaped 2 x len(df).

    """

    print( "Sorting data" )
    df.sort_values( by=["time", "particle id", "processor"], ascending=[True, True, True], inplace=True )
    processor_groups = df.groupby( "processor" )

    print( "Integrating data" )
    inputs = [(dataset, iterations, model, device) for _, dataset in processor_groups]

    with multiprocessing.Pool( processes=cores ) as pool:
        results = np.vstack( pool.starmap( do_iterative_inference, inputs ) )

    # TODO trim

    print( "Scaling data" )
    scaled_results = scale_droplet_parameters( results )

    return scaled_results

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

def standard_distance( truth_output, model_output ):
    """
    Calculates the "distance" between the `mlp` output and the true output
    on a dataframe. Uses the absolute value of the difference between the log
    of both output radii and the absolute value of the difference between the
    two output temperatures. Requires a dataframe with the mlp output.

    Takes 2 argumenst:

      truth_output        - NumPy array, shaped number 2 x time steps containing:
                            the backwards euler radius output at index 0
                            the backwards euler temperature output at index 1
      model_output        - NumPy array, shaped number 2 x time steps containing:
                            the MLP radius output at index 0
                            the MLP temperature output at index 1

    Returns 1 value:

      NTLP_distance_data  - NumPy array, shaped data 2 x time steps, containing
                            the absolute value of the difference of the log of
                            the output radii at index 0 the absolute value
                            difference of of the output temperatures at index 1.

    """

    return np.array( [np.abs( np.log10( model_output[:, 0] / truth_output[:, 0] ) ),
                      np.abs( model_output[:, 1] - truth_output[:, 1] )] ).T
