import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

from .data import create_droplet_batch
from .models import do_inference, do_iterative_inference
from .physics import dydt
from .scoring import calculate_nrmse, identity_norm

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
    fig_h, ax_h = plt.subplots( 2, 2, figsize=figure_size, sharex=True )
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
    fig_h, ax_h = plt.subplots( 2, 2, figsize=figure_size, sharex=True )
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

def analyze_model_particle_performance( times, truth_output, model_output, norm=None, title_string=None, figure_size=None, time_range=None ):
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
      norm            - Optional user provided norm to apply to all
                        data before calculating nrmse. Accepts an array
                        sized number_observations x 2 and returns a normed
                        array of the same length. Defaults to `identity_norm`.
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

    # Default to identity norm
    if norm is None:
        norm = identity_norm

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

    normed_truth_output = norm( truth_output )
    normed_model_output = norm( model_output )

    # Compute the normalized RMSE across the entire time scale.
    nrmse = calculate_nrmse( normed_truth_output, normed_model_output )

    # Create our figure and embed the parameters that were evaluated.
    fig_h, ax_h = plt.subplots( 2, 2, figsize=figure_size, sharex=True )
    fig_h.suptitle( "Droplet Size and Temperature\n NRMSE={:g}%\n{}".format( nrmse * 100, title_string) )

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
    ax_h[0][0].set_yscale( "log" )
    ax_h[0][1].set_yscale( "log" )

    if time_range is not None:
        ax_h[0][0].set_xlim( time_range )
        ax_h[0][1].set_xlim( time_range )
        ax_h[1][0].set_xlim( time_range )
        ax_h[1][1].set_xlim( time_range )

    fig_h.tight_layout()

    return fig_h, ax_h

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

    fig_h, ax_h = plt.subplots( 1, 2, figsize=(9, 4), sharex=True )

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

def plot_particles( particles_df, force_flag=False, time_range=[-np.inf, np.inf] ):
    """
    Plots one or more particles' characteristics along their lifetimes.  Creates
    a figure with 7 sub-plots, arranged in a vertical configuration, visualizing
    the following characteristics:

      1. Integration time, in seconds
      2. Output radius, in meters
      3. Output particle temperature, in Kelvin
      4. Background air temperature, in Kelvin
      5. Background relative humidity, as a percentage
      6. Salt mass, in kilograms
      7. Air density, in kg/m^3

    Takes 2 arguments:

      particles_df - Particles DataFrame whose particles should be visualized.
      force_flag   - Optional flag indicating that large numbers of particles
                     should be visualized.  If omitted, defaults to False and
                     ValueError is raised when too many particles are present
                     in particles_df.
      time_range   - Optional tuple containing the time bounds, lower and upper,
                     to plot particles.  If omitted, defaults to [-np.inf, np.inf]
                     and all particles are plotted.  If a particle in particles_df
                     does not have any observations in time_range then it is
                     not plotted.

    Returns 2 values:

      fig_h - The Matplotlib figure handle created.
      ax_h  - The Matplotlib axes handle containing the line plots.

    """

    # 10 particles in the same time window can get busy so set this as our soft
    # limit for plotting.
    MAXIMUM_NUMBER_OF_PARTICLES = 10

    if isinstance( particles_df, pd.Series ):
        particles_df = particles_df.to_frame().T

    # Blow up if we were provided too many particles unless we're forced to.
    # Things get *really* busy really with more than a handful of particles so
    # we don't want to plot them if it's not imperative.  Additionally, for very
    # large numbers of particles (hundreds or thousands) this can take quite a
    # long time to complete and we want to protect against accidentally plotting
    # an entire DataFrame and not just a subset of rows.
    if len( particles_df ) > MAXIMUM_NUMBER_OF_PARTICLES and not force_flag:
        raise ValueError( "Too many particles to plot ({:d})!".format(
            len( particles_df ) ) )

    fig_h, ax_h = plt.subplots( 7, 1, figsize=(8, 15), sharex=True )

    fig_h.suptitle( "{:d} Particle{:s}\n".format(
        len( particles_df ),
        "" if len( particles_df ) == 1 else "s"
    ) )

    # Initialize our temperature range to an infinite range that will get
    # necked down with each particle processed.
    temperature_min =  np.inf
    temperature_max = -np.inf

    # Plot each of the particles sequentially.
    for particle_index in range( len( particles_df ) ):
        particle = particles_df.iloc[particle_index]

        # Build this particle's timeline.
        times = particle["birth time"] + np.cumsum( particle["integration times"] ) - particle["integration times"][0]

        # Identify the observations of interest.
        times_mask = (times >= time_range[0]) & (times <= time_range[1])

        # Skip this particle if it's observations aren't within the requested
        # range.
        if times_mask.sum() == 0:
            continue

        # Keep our particle and air temperatures on the same scale.
        temperature_min = min( temperature_min,
                               min( particle["output temperatures"][times_mask].min(),
                                    particle["air temperatures"][times_mask].min() ) )
        temperature_max = max( temperature_max,
                               max( particle["output temperatures"][times_mask].max(),
                                    particle["air temperatures"][times_mask].max() ) )

        # Create label for this particle
        particle_label = str( particle.name )

        # Plot our quantities with labels
        ax_h[0].plot( times[times_mask], particle["integration times"][times_mask], label=particle_label )
        ax_h[0].set_title( "Integration Time" )
        ax_h[0].set_ylabel( "dt (s)" )

        ax_h[1].plot( times[times_mask], particle["output radii"][times_mask], label=particle_label )
        ax_h[1].set_title( "Output Radius" )
        ax_h[1].set_ylabel( "Size (m)" )

        ax_h[2].plot( times[times_mask], particle["output temperatures"][times_mask], label=particle_label )
        ax_h[2].set_title( "Particle Temperature" )
        ax_h[2].set_ylabel( "Kelvin" )

        ax_h[3].plot( times[times_mask], particle["air temperatures"][times_mask], label=particle_label )
        ax_h[3].set_title( "Air Temperature" )
        ax_h[3].set_ylabel( "Kelvin" )

        ax_h[4].plot( times[times_mask], particle["relative humidities"][times_mask] * 100, label=particle_label )
        ax_h[4].set_title( "Relative Humidity" )
        ax_h[4].set_ylabel( "Percentage (%)" )

        ax_h[5].plot( times[times_mask], particle["salt masses"][times_mask], label=particle_label )
        ax_h[5].set_title( "Salt Mass" )
        ax_h[5].set_ylabel( "kg" )

        ax_h[6].plot( times[times_mask], particle["air densities"][times_mask], label=particle_label )
        ax_h[6].set_title( "Air Density" )
        ax_h[6].set_ylabel( "kg/m^3" )
        ax_h[6].set_xlabel( "Time (s)" )

    # Set the particle and air temperatures to the same scale.
    ax_h[2].set_ylim( [temperature_min, temperature_max] )
    ax_h[3].set_ylim( [temperature_min, temperature_max] )

    # Turn on minor ticks on both the X- and Y-axis so its easier to interpret
    # the data.
    for axes_index in range( len( ax_h ) ):
        ax_h[axes_index].minorticks_on()

    # Create a single legend positioned outside the plot area on the right,
    # centered on the middle axis.
    ax_h[3].legend( bbox_to_anchor=(1.01, 0.5),
                    loc="center left",
                    title="Particles" )

    # Adjust layout to prevent legend from being cut off.  This makes room for
    # the legend on the right.
    plt.subplots_adjust( right=0.85 )
    fig_h.tight_layout()

    return fig_h, ax_h

def plot_particle_history( history_nc ):
    """
    Plots the particle-related variables in an NTLP history file.  Creates a
    figure with two horizontal subplots containing:

      1. Particle Events: The counts for new, destroyed, failed Gauss-Newton,
                          and failed iterative solve (BE), as a function of
                          time, line and scatter plots.
      2. Population Counts: The total particle count and breakdown between
                            aerosol and droplets, as a function of time, line
                            plots.

    NOTE: This does not visualize the number of activated/deactivated particles.
          This wasn't needed during verification.

    Takes 1 argument:

      history_nc - netCDF4 Dataset object containing the NTLP history variables.

    Returns 2 values:

      fig_h - The Matplotlib figure handle created.
      ax_h  - The Matplotlib axes handle containing the line plots.

    """

    simulation_time      = history_nc.variables["time"][:]
    number_particles     = history_nc.variables["tnumpart"][:]
    number_droplets      = history_nc.variables["tnumdrop"][:]
    number_aerosols      = history_nc.variables["tnumaerosol"][:]
    number_destroyed     = history_nc.variables["tnum_destroy"][:]
    number_gn_failed     = history_nc.variables["tnum100"][:]
    number_impossible    = history_nc.variables["tnumimpos"][:]
    number_new_particles = history_nc.variables["tot_reintro"][:]

    fig_h, ax_h = plt.subplots( 1, 2, figsize=(10, 6) )

    fig_h.suptitle( "Particle Measurements in '{:s}".format( history_nc.filepath() ) )

    #
    # NOTE: The number of new and failed Gauss-Newton (GN) particles are scatter
    #       plots so as to not occlude important data with their regular swings
    #       in value.  Failed GN will always be larger than failed iterative
    #       solves but we don't want times where there weren't any failures
    #       (zero value) to clutter the plot with a dip.  Similarly for new
    #       particles as injections happen at a fixed cadence resulting in
    #       regular dips covering everything else in the plot (as this is the
    #       largest magnitude value visualized).
    #
    ax_h[0].plot( simulation_time, number_new_particles, ".", label="New Particles" )
    ax_h[0].plot( simulation_time, number_destroyed,          label="Destroyed Particles" )
    ax_h[0].plot( simulation_time, number_gn_failed, ".",     label="Gauss-Newton Failed" )
    ax_h[0].plot( simulation_time, number_impossible,         label="Failed Iterative Solves" )

    # Put the legend in the top-right corner.
    #
    # NOTE: There isn't a great place to put it and this is the least worst
    #       option.  This likely occludes the new particles plot though its
    #       trend should be understandable from the visible portions.
    #
    ax_h[0].legend( loc="upper right" )
    ax_h[0].set_ylabel( "Count" )
    ax_h[0].set_xlabel( "Time (s)" )
    ax_h[0].set_title( "Particle Events" )

    ax_h[1].plot( simulation_time, number_droplets,  label="Droplet Count" )
    ax_h[1].plot( simulation_time, number_aerosols,  label="Aerosol Count" )
    ax_h[1].plot( simulation_time, number_particles, label="Particle Count" )

    ax_h[1].legend( loc="center right" )
    ax_h[1].set_xlabel( "Time (s)" )
    ax_h[1].set_title( "Population Counts" )

    return fig_h, ax_h
