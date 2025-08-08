import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .data import BE_TAG_NAME, \
                  get_evaluation_column_names
from .models import do_iterative_inference
from .physics import BDF_TOLERANCE_ABSOLUTE, \
                     BDF_TOLERANCE_RELATIVE, \
                     dydt, \
                     get_parameter_ranges, \
                     timed_solve_ivp
from .scoring import calculate_nrmse

from itertools import islice

def plot_droplet_size_temperatures( times, size_temperatures, background_parameters={},
                                    compare=None, ax_h=None, title_string=None ):
    """
    Generic function for plotting radius/temperature data alongside background parameters.
    Can graph one or more time series of radius/temperature data and compare them for
    relative/absolute difference. The first entry is `size_temperatures` is treated as
    the baseline for comparison.

    Takes 6 Arguments:

      times                 - Array, contains the times corresponding to the provided
                              time series data.
      size_temperatures     - Dictionary contains key:value pairs of
                              time_series_label:time_series_data to graph/compare.
      background_parameters - Optional dictionary, contains key:value pairs of
                              data_label:time_series_data for additional background
                              parameters to plot.
      compare               - Optional Boolean, determines whether to generate
                              absolute/relative difference plots between radius/temperature
                              data provided. Defaults to False if one radius/temperature
                              time series is provided and True if two or more are provided.
      ax_h                  - Optional axes array to graph onto.
      title_string          - Optional String to title the plot. Defaults to "Droplet Size
                              and Temperature."

    Returns 2 Values:

      fig_h - Figure generated for the plots. Equals none if ax_h is provided
              since no new plot is generated.
      ax_h  - Array of axes that were used for plotting.

    """

    if compare is None:
        compare = len( size_temperatures ) > 1
    if title_string is None:
        title_string = "Droplet Size and Temperature"

    time_series_count = len( size_temperatures )

    # Determine programmatically how many rows are needed for the plot.
    subplot_height = 2 if compare else 1
    if background_parameters is not None:
        subplot_height += (len( background_parameters ) + 1) // 2

    if ax_h is None:
        fig_h, ax_h = plt.subplots( subplot_height, 2,
                                    figsize=(10, 3.2*subplot_height + 0.5),
                                    sharex=True,
                                    layout="constrained" )
        fig_h.suptitle( title_string )
    else:
        fig_h = None

    # If there's only one row of axes, make it a singleton list so it can be
    # handled like a multi-row figure.
    if subplot_height == 1:
        ax_h = np.array( [ax_h] )

    # Plot the radius and temperature data.  Use a qualitative colormap with 9
    # darker colors so each of the evaluations is distinct from each other.
    cmap   = plt.get_cmap( "Set1" )
    colors = cmap( np.linspace( 0.0, 1.0, time_series_count ) )
    for color, (label, time_series_data) in zip( colors, size_temperatures.items() ):
        ax_h[0][0].plot( times, time_series_data[..., 0], label=label, color=color )
        ax_h[0][1].plot( times, time_series_data[..., 1], label=label, color=color )

    ax_h[0][0].set_ylabel( "Radius (m)" )
    ax_h[0][0].set_yscale( "log" )
    ax_h[0][1].set_ylabel( "Temperature (K)" )

    ax_h[0][0].legend()
    ax_h[0][1].legend()

    if compare:
        if time_series_count == 1:
            raise( Exception( "Error: compare flag true but no other time series to compare!" ) )

        # Setup twin axes so we can see both relative and absolute differences
        # between the reference and each comparison with a single plot.
        ax_h_twin_radius      = ax_h[1][0].twinx()
        ax_h_twin_temperature = ax_h[1][1].twinx()

        # Get the first label and datapoint to use a reference in the
        # comparisons.  Plot the remaining relative to it.
        reference_label, reference_data = next( iter( size_temperatures.items() ) )
        for color, (label, comparison_data) in islice( zip( colors, size_temperatures.items() ),
                                                       1,
                                                       None ):
            ax_h[1][0].plot( times,
                             (np.abs( reference_data[:, 0] - comparison_data[:, 0] ) /
                              reference_data[:, 0] * 100),
                             color=color,
                             label="{:s} relative error".format( label ) )
            ax_h[1][1].plot( times,
                             (np.abs( reference_data[:, 1] - comparison_data[:, 1] ) /
                              reference_data[:, 1] * 100),
                             color=color,
                             label="{:s} error".format( label ) )
            ax_h_twin_radius.plot( times,
                                   np.abs( reference_data[:, 0] - comparison_data[:, 0] ),
                                   c=color,
                                   label="{:s} absolute error".format( label ),
                                   linestyle="dashed" )

        # Label the axes.
        ax_h[1][0].set_title( "Radius Error against {:s}".format( reference_label ) )
        ax_h[1][1].set_title( "Temperature Error against {:s}".format( reference_label ) )
        ax_h[1][0].set_ylabel( "Relative Difference (%)" )
        ax_h[1][1].set_ylabel( "Relative Difference (%)" )
        ax_h_twin_temperature.set_ylabel( "Absolute Difference (K)" )
        ax_h_twin_radius.set_ylabel( "Absolute Difference (m)" )

        # Set the temperature graph's absolute difference axes to match the
        # relative axes since we didn't plot anything.
        #
        # NOTE: We take care to compute the extrema without NaN's by using masks
        #       as we're not interested in .nanmin()/.nanmax()'s RuntimeWarning.
        #       It's also worth noting that we're guaranteed a finite range here
        #       as all NaN's will have already exploded above.  As a result we
        #       provide a bogus initial argument as it won't be used.
        #
        temperature_relative_limits = ax_h[1][1].get_ylim()
        reference_minimum           = np.amin( reference_data[:, 1], where=~np.isnan( reference_data[:, 1] ), initial=0.0 )
        reference_maximum           = np.amax( reference_data[:, 1], where=~np.isnan( reference_data[:, 1] ), initial=0.0 )
        temperature_absolute_limits = (temperature_relative_limits[0] / 100.0 * reference_minimum,
                                       temperature_relative_limits[1] / 100.0 * reference_maximum)
        ax_h_twin_temperature.set_ylim( temperature_absolute_limits )

        ax_h[1][0].legend( loc=(0.05, 0.75) )
        ax_h[1][1].legend()
        ax_h_twin_radius.legend( loc=(0.05, 0.85) )

        # Turn on minor ticks for the twins to match what we do for the other
        # axes below.
        ax_h_twin_radius.minorticks_on()
        ax_h_twin_temperature.minorticks_on()

    # Plot background parameters in the remaining subplots.
    starting_index = 4 if compare else 2
    for index, (label, time_series) in enumerate( background_parameters.items() ):
       axis_row_index    = (starting_index + index) // 2
       axis_column_index = (starting_index + index) % 2
       current_axis      = ax_h[axis_row_index][axis_column_index]

       current_axis.set_title( label )
       current_axis.plot( times, time_series )
       current_axis.set_ylabel( label )

    # Format all subplots.
    #
    # NOTE: This does not update the comparisons' twin axes!
    #
    for current_axis in ax_h.flat:
        current_axis.minorticks_on()
        current_axis.grid( color="b", alpha=0.1 )
        current_axis.set_xlabel( "Time (s)" )

    return fig_h, ax_h

def plot_droplet_size_temperatures_dataframe( particle_dataframe, evaluation_tags, **kwargs ):
    """
    Wrapper for plot_droplet_size_temperatures. Plots radius/temperatures from a dataframe
    for specified evaluation tags.

    Takes 2 Arguments:

      particle_dataframe - Pandas DataFrame containing the row for the particle to plot.
      evaluation_tags     - String List or String containing the evaluation tags to compare.
      **kwargs            - Additional arguments to pass to plot_droplet_size_temperatures.

    Returns 2 Values:

      fig_h - Figure generated for the plots. Equals none if ax_h is provided
              since no new plot is generated.
      ax_h  - Array of axes that were used for plotting.

    """

    # If just one tag was provided, wrap it into an array
    if type( evaluation_tags ) == str:
        evaluation_tags = [evaluation_tags]

    # Add a default title to the figure if the user did not provide one.
    if "title_string" not in kwargs:
        evaluation_string = ", ".join( [evaluation_tag for evaluation_tag in evaluation_tags] )
        particle_id       = particle_dataframe.name
        kwargs["title_string"] = ("Droplet size/temperatures for {:s}\n" +
                                  "Particle {:d}").format( evaluation_string, particle_id )

    # Get the data to plot as series.
    times             = particle_dataframe["times"]
    size_temperatures = {
        evaluation_tag: np.stack( particle_dataframe[["output {:s} radii".format( evaluation_tag ),
                                                      "output {:s} temperatures".format( evaluation_tag )]],
                                  axis=-1 )
        for evaluation_tag in evaluation_tags
    }

    fig_h, ax_h = plot_droplet_size_temperatures( times, size_temperatures, **kwargs )

    return fig_h, ax_h

def plot_droplet_size_temperatures_domain( input_parameters, model=None, dt=None, final_time=10.0, **kwargs ):
    """
    Wrapper for plot_droplet_size_temperatures. Evaluates the BDF radius/temperature solution
    for the given background conditions until final_time. If a model is provided, evalutes the model
    iteratively with a time step of dt. Plots the resulting radius/temperature data.

    Takes 5 Arugments:

      input_parameters - Array sized 6 containing the initial droplet parameters.
      model            - Optional PyTorch model to compare with BDF.
      dt               - Optional float determining the time step for model evaluation.
                         Defaults to the mean of the log time range.
      final_time       - Float determining the time range of the trajectory. Defaults to 10s.
      **kwargs         - Additional parameters to pass to plot_droplet_size_temperatures.

    Returns 2 Values:

      fig_h - Figure generated for the plots. Equals none if ax_h is provided
              since no new plot is generated.
      ax_h  - Array of axes that were used for plotting.

    """

    # Use a default integration time in the middle of our temporal range if one
    # wasn't provided.
    if dt is None:
        dt = 10**(np.mean( get_parameter_ranges()["time"] ))

    NUMBER_TIME_POINTS = int( final_time // dt )
    t_eval             = dt * np.arange( 0, NUMBER_TIME_POINTS )

    print( "Inputs: {}\n"
           "dt: {}".format(
               input_parameters,
               dt ) )

    # Get truth from the ODEs.
    y0         = (input_parameters[0], input_parameters[1])
    bdf_output = timed_solve_ivp( dydt,
                                  [0, final_time],
                                  y0,
                                  method="BDF",
                                  t_eval=t_eval,
                                  args=(input_parameters[2:],),
                                  atol=BDF_TOLERANCE_ABSOLUTE,
                                  rtol=BDF_TOLERANCE_RELATIVE )

    # Package the BDF outputs for visualization.
    size_temperatures = { "bdf": bdf_output }

    # Plot the model's estimate if a model was provided.
    if model is not None:
        model_output = do_iterative_inference( np.tile( input_parameters,
                                                        (NUMBER_TIME_POINTS, 1) ),
                                               t_eval,
                                               model,
                                               "cpu" )
        model_nrmse  = calculate_nrmse( bdf_output, model_output )

        # Add the MLP output to the visualization.
        size_temperatures["model"] = model_output

    # Add a default title to the figure if the user did not provide one.
    if "title_string" not in kwargs:
        kwargs["title_string"] = ("Droplet Size and Temperature\n" +
                                  "Radius={:g}, Temperature={:g}, m_s={:g}, Air Temp={:g}, RH={:g}, rhoa={:g}").format(
                                      input_parameters[0],
                                      input_parameters[1],
                                      input_parameters[2],
                                      input_parameters[3],
                                      input_parameters[4],
                                      input_parameters[5])
        if model is not None:
            kwargs["title_string"] += ("\n" +
                                       "Model NRMSE: {:.3e}").format( model_nrmse )

    fig_h, ax_h = plot_droplet_size_temperatures( t_eval, size_temperatures, **kwargs )

    return fig_h, ax_h

def plot_droplet_size_temperatures_scoring( particle_dataframe, score_report, **kwargs ):
    """
    Wrapper for plot_droplet_size_temperature_dataframe. Plots the deviations recorded
    in score_report alongside the corresponding radius/temperature data.

    Takes 3 Arguments:

      particle_dataframe - Pandas DataFrame containing the row for the particle to plot.
      score_report       - ScoringReport object to plot deviations from.
      **kwargs           - Additional parameters to pass to plot_droplet_size_temperatures.

    Returns 2 Values:

      fig_h - Figure generated for the plots. Equals none if ax_h is provided
              since no new plot is generated.
      ax_h  - Array of axes that were used for plotting.

    """

    # The score report has the evaluation tags for our DataFrame.
    reference_evaluation_tag  = next( iter( score_report.reference_evaluation ) )
    comparison_evaluation_tag = next( iter( score_report.comparison_evaluation ) )

    # Add a default title to the figure if the user did not provide one.
    if "title_string" not in kwargs:

        particle_id    = particle_dataframe.name
        particle_nrmse = score_report.per_particle_nrmse[particle_id]

        # Provide the comparison's details and results.
        kwargs["title_string"] = ("Scoring Particle {:d} for {:s} vs. {:s}\n" +
                                  "NRMSE: {:.3e}\n" +
                                  "CUSUM Tolerance: {}, CUSUM Threshold: {}").format(
                                      particle_id,
                                      reference_evaluation_tag,
                                      comparison_evaluation_tag,
                                      particle_nrmse,
                                      score_report.cusum_error_tolerance,
                                      score_report.cusum_error_threshold )

    fig_h, ax_h = plot_droplet_size_temperatures_dataframe( particle_dataframe,
                                                            [reference_evaluation_tag,
                                                             comparison_evaluation_tag],
                                                            **kwargs )

    # We use a qualitative colormap with 9 distinct colors so we have enough for
    # datasets with a large number of deviation clusters.
    cmap = plt.get_cmap( "Set1" )

    # Plot each of the deviations found.
    for deviation_index in np.where( score_report.deviation_particle_ids == particle_dataframe.name )[0]:
        deviation_parameter = score_report.deviation_parameters[deviation_index]
        deviation_time      = score_report.deviation_times[deviation_index]
        deviation_cluster   = score_report.deviation_clusters[deviation_index]
        deviation_label     = "{:s} deviation, cluster {:d}".format( deviation_parameter.name.lower(),
                                                                     deviation_cluster )

        for ax in ax_h.flat:
            ax.axvline( x=deviation_time, linewidth=1, linestyle="--", label=deviation_label,
                        color=cmap( deviation_cluster ) )

    # Update the legends, if any, so as to include the deviations' labels.
    for ax in ax_h[0]:
        ax.legend()

    return fig_h, ax_h

def plot_particles( particles_df, force_flag=False, time_range=[-np.inf, np.inf], evaluation_tags=[] ):
    """
    Plots one or more particles' characteristics along their lifetimes.  Creates
    a figure with 7 sub-plots, arranged in a vertical configuration, visualizing
    the following characteristics:

      1. Integration time, in seconds
      2. Output BE radius, in meters
      3. Output BE particle temperature, in Kelvin
      4. Background air temperature, in Kelvin
      5. Background relative humidity, as a percentage
      6. Salt mass, in kilograms
      7. Air density, in kg/m^3

    Can plot multiple evaluations for each particle to analyze model performance.

    Takes 4 arguments:

      particles_df    - Particles DataFrame whose particles should be visualized.
      force_flag      - Optional flag indicating that large numbers of particles
                        should be visualized.  If omitted, defaults to False and
                        ValueError is raised when too many particles are present
                        in particles_df.
      time_range      - Optional tuple containing the time bounds, lower and upper,
                        to plot particles.  If omitted, defaults to [-np.inf, np.inf]
                        and all particles are plotted.  If a particle in particles_df
                        does not have any observations in time_range then it is
                        not plotted.
      evaluation_tags - Optional list of evaluation tags specifying evaluations
                        to plot for each particle.  If omitted, defaults to an
                        empty list which visualizes the backward Euler (BE)
                        observations.  When provided, *only* the tags specified
                        are visualized which will not include the BE observations
                        unless included in evaluation_tags.

    Returns 2 values:

      fig_h - The Matplotlib figure handle created.
      ax_h  - The Matplotlib axes handle containing the line plots.

    """

    def _get_particle_label( particle_df ):
        """
        Constructs a string label for the supplied particle DataFrame from its
        particle identifier.

        Takes 1 argument:

          particle_df - Particle DataFrame to label.

        Returns 1 value:

          particle_label - Label associated with the particle supplied.

        """

        return str( particle_df.name )

    def _get_evaluation_label( particle_label, evaluation_tag ):
        """
        Constructs a label for a specific evaluation tag given the provided
        particle label.

        Takes 2 arguments:

          particle_label - Particle label string.
          evaluation_tag - Evaluation tag string.

        Returns 1 value:

          evaluation_label - Label associated with the particle's evaluation.

        """

        return "{:s} - {:s}".format(particle_label, evaluation_tag )

    # Limit ourselves to the number of unique colors Matplotlib has with the
    # default color scheme by default.  More particles than this will reuse
    # colors making it difficult to distinguish what is what.
    # Get the default Matplotlib color sequence.  We cycle
    colors_list                 = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    MAXIMUM_NUMBER_OF_PARTICLES = len( colors_list )

    if isinstance( particles_df, pd.Series ):
        particles_df = particles_df.to_frame().T

    # Visualize BE for each particle if no evaluations were requested, otherwise
    # make it the first evaluation plotted for consistency.  We build a map
    # from the requested evaluation tags to the DataFrame column names that
    # provide the data to plot.
    #
    # NOTE: Care is taken to add the tags to the map so that BE is first when
    #       requested, so as to match the plot colors across the different axes.
    #       This requires Python 3.6 or newer.
    #
    evaluations_map = {}
    if len( evaluation_tags ) == 0 or BE_TAG_NAME in evaluation_tags:
        evaluations_map[BE_TAG_NAME] = get_evaluation_column_names( BE_TAG_NAME )

    for evaluation_tag in evaluation_tags:
        evaluations_map[evaluation_tag] = get_evaluation_column_names( evaluation_tag )

    # Blow up if we were provided too many particles unless we're forced to.
    # Things get *really* busy really with more than a handful of particles so
    # we don't want to plot them if it's not imperative.  Additionally, for very
    # large numbers of particles (hundreds or thousands) this can take quite a
    # long time to complete and we want to protect against accidentally plotting
    # an entire DataFrame and not just a subset of rows.
    if len( particles_df ) * len( evaluations_map ) > MAXIMUM_NUMBER_OF_PARTICLES and not force_flag:
        raise ValueError( "Too many particles to plot ({:d} particle{:s} each with {:d} evaluation{:s})!".format(
            len( particles_df ),
            "" if len( particles_df ) == 1 else "s",
            len( evaluations_map ),
            "" if len( evaluations_map ) == 1 else "s" ) )

    fig_h, ax_h = plt.subplots( 7, 1, figsize=(8, 15), sharex=True )

    fig_h.suptitle( "{:d} Particle{:s}\n".format(
        len( particles_df ),
        "" if len( particles_df ) == 1 else "s"
    ) )

    # Initialize our temperature range to an infinite range that will get
    # necked down with each particle processed.
    temperature_min =  np.inf
    temperature_max = -np.inf

    # We need our line plots' colors to match so that the environmental plots'
    # (air temperature, relative humidity, etc) colors are synchronized with the
    # BE evaluations colors (or the first evaluation when BE isn't present).
    # This is trivial when BE is the only evaluation though gets complicated
    # when there are additional.
    #
    # We build a mapping plot label names to colors to avoid complicated modular
    # arithmetic while keeping the logic to handle the environmental/evaluation
    # asymmetry as minimal as possible.
    colors_map = {}
    for particle_index in range( len( particles_df ) ):
        # Get the particle's label and its (equivalent) first evaluation label.
        particle_label         = _get_particle_label( particles_df.iloc[particle_index] )
        first_evaluation_label = _get_evaluation_label( particle_label,
                                                        list( evaluations_map.keys() )[0] )

        # Figure out where this particle's colors start.
        environmental_color_index = (particle_index * len( evaluations_map )) % MAXIMUM_NUMBER_OF_PARTICLES

        # The particle's label, used for environmental plots, and the first
        # evaluation's label share the same color so users can understand which
        # particles traversed what environment.
        colors_map[particle_label]         = colors_list[environmental_color_index]
        colors_map[first_evaluation_label] = colors_list[environmental_color_index]

        # The remaining evaluation labels use the colors available before the
        # next particle's color.
        for evaluation_tag_index in range( 1, len( evaluations_map ) ):
            evaluation_label = _get_evaluation_label( particle_label,
                                                      list( evaluations_map.keys() )[evaluation_tag_index] )

            evaluation_color_index       = (environmental_color_index + evaluation_tag_index) % MAXIMUM_NUMBER_OF_PARTICLES
            colors_map[evaluation_label] = colors_list[evaluation_color_index]

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
        for evaluation_tag in evaluations_map:
            evaluation_temperatures_name = evaluations_map[evaluation_tag][1]

            temperature_min = min( temperature_min,
                                   particle[evaluation_temperatures_name][times_mask].min() )
            temperature_max = max( temperature_max,
                                   particle[evaluation_temperatures_name][times_mask].max() )

        # Create label for this particle.
        particle_label = _get_particle_label( particle )

        # Plot our quantities with labels
        ax_h[0].plot( times[times_mask],
                      particle["integration times"][times_mask],
                      color=colors_map[particle_label],
                      label=particle_label )
        ax_h[0].set_title( "Integration Time" )
        ax_h[0].set_ylabel( "dt (s)" )

        # Plot the radii and temperatures for each of the evaluations.
        for evaluation_tag in evaluations_map:
            (evaluation_radii_name,
             evaluation_temperatures_name) = evaluations_map[evaluation_tag]

            evaluation_label = _get_evaluation_label( particle_label,
                                                      evaluation_tag )

            #
            # NOTE: These are labeled with the evaluation tag as a suffix.
            #
            ax_h[1].plot( times[times_mask],
                          particle[evaluation_radii_name][times_mask],
                          color=colors_map[evaluation_label],
                          label=evaluation_label )
            ax_h[2].plot( times[times_mask],
                          particle[evaluation_temperatures_name][times_mask],
                          color=colors_map[evaluation_label],
                          label=evaluation_label )

        ax_h[1].set_title( "Output Radius" )
        ax_h[1].set_ylabel( "Size (m)" )

        ax_h[2].set_title( "Particle Temperature" )
        ax_h[2].set_ylabel( "Kelvin" )

        # Now plot the environments for the particle.
        ax_h[3].plot( times[times_mask],
                      particle["air temperatures"][times_mask],
                      color=colors_map[particle_label],
                      label=particle_label )
        ax_h[3].set_title( "Air Temperature" )
        ax_h[3].set_ylabel( "Kelvin" )

        ax_h[4].plot( times[times_mask],
                      particle["relative humidities"][times_mask] * 100,
                      color=colors_map[particle_label],
                      label=particle_label )
        ax_h[4].set_title( "Relative Humidity" )
        ax_h[4].set_ylabel( "Percentage (%)" )

        ax_h[5].plot( times[times_mask],
                      particle["salt masses"][times_mask],
                      color=colors_map[particle_label],
                      label=particle_label )
        ax_h[5].set_title( "Salt Mass" )
        ax_h[5].set_ylabel( "kg" )

        ax_h[6].plot( times[times_mask],
                      particle["air densities"][times_mask],
                      color=colors_map[particle_label],
                      label=particle_label )
        ax_h[6].set_title( "Air Density" )
        ax_h[6].set_ylabel( "kg/m^3" )
        ax_h[6].set_xlabel( "Time (s)" )

    # Adjust our temperature bounds *very* slightly so our line plots don't
    # touch the edge of the axes.
    #
    # NOTE: We choose 1/20th of a percent of 300K is 0.15K.  Enough to give a
    #       slight amount of padding but not enough to change the zoom on the
    #       temperature plots so as to preserve their dynamics.
    #
    temperature_min *= 0.9995
    temperature_max *= 1.0005

    # Set the particle and air temperatures to the same scale.
    ax_h[2].set_ylim( [temperature_min, temperature_max] )
    ax_h[3].set_ylim( [temperature_min, temperature_max] )

    # Turn on minor ticks on both the X- and Y-axis so its easier to interpret
    # the data.
    for axes_index in range( len( ax_h ) ):
        ax_h[axes_index].minorticks_on()

    # Create a single legend positioned outside the plot area on the right,
    # centered on the middle axis.
    ax_h[2].legend( bbox_to_anchor=(1.01, 0.5),
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
