import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .data import BE_TAG_NAME, \
                  get_evaluation_column_names, \
                  insert_timeseries_gaps
from .models import do_iterative_inference
from .physics import BDF_TOLERANCE_ABSOLUTE, \
                     BDF_TOLERANCE_RELATIVE, \
                     dydt, \
                     get_parameter_ranges, \
                     timed_solve_ivp
from .scoring import calculate_nrmse

from itertools import islice

def average_particles_data( particles_df, evaluation_tags, simulation_times, background_columns=[] ):
    """
    Gives radius/temperatures averages for every each time step across all particles
    for given list of evaluation tags.

    Takes 4 Arguments:
      particles_df       - Pandas DataFrame of particle data to average
      evaluation_tags    - String list of evaluation tags to average
      simulation_times   - Numpy Array of the time steps in the simulation to average.
      background_columns - Optional String List containing list of other columns
                           to average

    Returns 2 Values:
      rt_averages         - Dictionary of evaluation_tag: averages, where averages is an array with:
                            at index 0: an array of radius averages for the given evaluation tag at
                                each time step.
                            at index 1: an array of temperature averages for the given evaluation tag
                                at each time step.
      background_averages - Dictionary of background column name: averages, where averages is an array
                            of the average value of that variable at each time step.
    """

    rt_averages         = {
        evaluation_tag: np.zeros( shape=(simulation_times.shape[0], 2) )
        for evaluation_tag in evaluation_tags
    }
    background_averages = {
        column: np.zeros( shape=( simulation_times.shape[0] ) )
        for column in background_columns
    }
    particle_counts = np.zeros( shape=simulation_times.shape[0], dtype=np.int32 )

    # Adds up the values for each evaluation tag
    # and each background column for a particle
    # at each time step.
    def _sum_particle( particle_df ):
        EPSILON = 0.0005

        # Identify the time indexes in the particle's
        # lifetime that are being averaged
        time_indexes = np.searchsorted( simulation_times, particle_df["times"] - EPSILON )

        mask = np.abs( simulation_times[[time_indexes]][0] - particle_df["times"] ) < EPSILON
        time_indexes = time_indexes[mask]

        # Sum radius/temperature and background columns for
        # relavent times
        particle_counts[time_indexes] += 1
        for evaluation_tag in evaluation_tags:
                rt_averages[evaluation_tag][time_indexes, 0] += (
                    particle_df["output {:s} radii".format( evaluation_tag )][mask])
                rt_averages[evaluation_tag][time_indexes, 1] += (
                    particle_df["output {:s} temperatures".format( evaluation_tag )][mask])

        for column in background_columns:
            background_averages[column][time_indexes] += particle_df[column][mask]

    particles_df.apply( _sum_particle, axis=1 )

    # Average the results
    for evaluation_tag in evaluation_tags:
        rt_averages[evaluation_tag][:, 0] /= particle_counts
        rt_averages[evaluation_tag][:, 1] /= particle_counts

    for _, average in background_averages.items():
        average /= particle_counts

    return rt_averages, background_averages

def bin_particles_data( particles_df, evaluation_tags, histogram_times, radbins, tempbins ):
    """
    Generates histogram counts for radius/temperature at the specified times.

    Takes 5 arguments:
      particles_df    - DataFrame row corresponding to a particular particle.
      evaluation_tags - List of evaluation tags to average.
      histogram_times - Array of times to generate histograms.
      radbins         - Array of radius histogram bin edges.
      tempbins        - Array of temperature histogram bin edges.

    Returns 1 value:
      histograms - Dictionary. Keys are evaluation tags. Values are each a list containing:
                   At index 0: Radius histograms. Array sized histogram_times.shape[0] x radbins.shape[0].
                               Contains the counts for each bin at each histogram time.
                   At index 1: Temperature histograms. Array sized histogram_times.shape[0] x tempbins.shape[0].
                               Contains the counts for each bin at each histogram time.
    """
    histograms = {
            evaluation_tag: [np.zeros( shape=(histogram_times.shape[0], radbins.shape[0]) ),
                             np.zeros( shape=(histogram_times.shape[0], tempbins.shape[0]) )]
            for evaluation_tag in evaluation_tags
    }

    def _bin_particle( particle_df ):
        """
        Bins a particular particle's radius/temperature for each histogram time.
        """

        # Skip particles that don't have observations.
        if particle_df["number evaluations"] == 0:
            return

        # Find the histogram times that fall in a particle's lifetime
        # with some rounding to account for any rounding errors.
        EPSILON = 0.0005
        timewindow_start         = np.searchsorted( histogram_times,
                                                    particle_df["times"][0]  + EPSILON ) - 1
        timewindow_end           = np.searchsorted( histogram_times,
                                                    particle_df["times"][-1] + EPSILON ) - 1
        relavent_histogram_times = histogram_times[timewindow_start:timewindow_end]

        # Find the indexes in the particles lifetime that correspond to each histogram
        time_indexes      = np.searchsorted( particle_df["times"], relavent_histogram_times)
        histogram_indexes = np.arange( timewindow_start, timewindow_end )

        # BE failures can cause misalignments in the times array. Remove any "matched" times that are
        # more than EPISLON different from the stipulated histogram time
        mask                     = np.where( np.abs( relavent_histogram_times
                                                   - particle_df["times"][[time_indexes]][0] ) < EPSILON )
        time_indexes             = time_indexes[mask]
        relavent_histogram_times = relavent_histogram_times[mask]
        histogram_indexes        = histogram_indexes[mask]

        if relavent_histogram_times.shape[0] == 0:
            return

        for evaluation_tag, hist in histograms.items():
            # Record the radius/temperature for each histogram time
            hist_radii        = particle_df["output {:s} radii".format(
                                            evaluation_tag )][[time_indexes]][0]
            hist_temperatures = particle_df["output {:s} temperatures".format(
                                            evaluation_tag )][[time_indexes]][0]

            for data_index, histogram_index in enumerate( histogram_indexes ):
                # Increment the bin corresponding to the radius/temperature
                hist[0][histogram_index, np.searchsorted( radbins, hist_radii[data_index] )]         += 1
                hist[1][histogram_index, np.searchsorted( tempbins, hist_temperatures[data_index] )] += 1

    particles_df.apply( _bin_particle, axis=1 )

    return histograms

def get_particles_data_simulation_times( particles_df ):
    """
    Generates one "simulation times" numpy array that includes
    all the time steps for all of the particles supplied.

    Takes 1 Argument:
      particles_df - DataFrame containing all of the particles from
                     which to construct a timeline.

    Returns 1 Argument:
      simulation_times - Numpy Array containing all time steps present
                         in the DataFrame.
    """
    simulation_times = particles_df.iloc[0]["times"]

    def _extend_timeline( particle_df ):
        # This used to be fancier - it would just extend the array if there were
        # missing steps at the beginning or end. However, since there are occasionally
        # gaps in our timelines due to bugs in Spray, we have to do the quick and
        # dirty method of just intersecting all of the arrays.
        nonlocal simulation_times
        simulation_times = np.union1d( simulation_times, particle_df["times"] )

    particles_df[1:].apply( _extend_timeline, axis=1 )

    return np.array( simulation_times )

def plot_droplet_size_temperatures( times, size_temperatures, background_parameters={},
                                    compare_flag=None, ax_h=None, title_string=None ):
    """
    Plots a single particle's radius and temperature data and, optionally, their
    associated background parameters.  Also plots the absolute and relative
    differences (errors) relative to the first radius and temperature time series
    when multiple are provided, one per additional time series.

    Data are plotted as a grid of N x 2, rows by columns, with the following
    layout:

                         Columns
     Row        1                         2                         Notes
     ---   ------------            -------------------       ----------------------
      1    Radii series            Temperature series
      2    Radii differences       Temperature differences   Not present when compare_flag == False,
      i    Background series 1     Background series 2
     i+1   Background series 3     Background series 4       Only when len( background_parameters ) > 2
     ...
     i+j   Background series J-1   Background series J       j = (J+1) // 2 - 1

    More succinctly, the total number of rows is defined:

        N = 1 + int( compare_flag ) + ((J + 1) // 2)

    Where there are J background time series.

    Raises ValueError if only a single radius and temperature time series are
    provided but comparison was requested.

    Takes 6 arguments:

      times                 - 1D array, containing the times corresponding to the
                              provided time series data.
      size_temperatures     - Dictionary contains timeseries labels and their
                              associated data to graph.  Each timeseries is an array
                              shaped N x 2, where N is the length of times.
      background_parameters - Optional dictionary, contains labels and time series
                              data for additional background parameters to plot.
                              If omitted, only time series in size_temperatures are
                              visualized.
      compare_flag          - Optional Boolean, determines whether to generate
                              absolute/relative difference plots between
                              radius/temperature data provided.  If omitted,
                              defaults to False when one radius/temperature time
                              series is provided and True if two or more are
                              provided.
      ax_h                  - Optional axes array to use.  If omitted, defaults to None
                              and a new figure with N * 2 axes is created.
      title_string          - Optional string to title the plot with.  If omitted,
                              defaults to "Droplet Size and Temperature."

    Returns 2 values:

      fig_h - Figure generated for the plots.  None if ax_h is provided since a new
              figure isn't created.
      ax_h  - List of list of axes that were used for plotting.  The outer list
              (first index) specifies the row while the inner list (second
              index) specifies the column.

    """

    if compare_flag is None:
        compare_flag = len( size_temperatures ) > 1
    if title_string is None:
        title_string = "Droplet Size and Temperature"

    time_series_count = len( size_temperatures )

    # Determine programmatically how many rows are needed for the plot.
    subplot_height = 2 if compare_flag else 1
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

    if compare_flag:
        if time_series_count == 1:
            raise( ValueError( "Error: compare flag true but no other time series to compare!" ) )

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
                             label="{:s} relative".format( label ) )
            ax_h[1][1].plot( times,
                             (np.abs( reference_data[:, 1] - comparison_data[:, 1] ) /
                              reference_data[:, 1] * 100),
                             color=color,
                             label="{:s} error".format( label ) )
            ax_h_twin_radius.plot( times,
                                   np.abs( reference_data[:, 0] - comparison_data[:, 0] ),
                                   c=color,
                                   label="{:s} absolute".format( label ),
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

        # Collect the lines and labels from both radius comparison axes so we
        # can create a single legend containing both.
        radius_relative_lines, radius_relative_labels = ax_h[1][0].get_legend_handles_labels()
        radius_absolute_lines, radius_absolute_labels = ax_h_twin_radius.get_legend_handles_labels()

        radius_lines  = radius_relative_lines + radius_absolute_lines
        radius_labels = radius_relative_labels + radius_absolute_labels

        ax_h[1][0].legend( radius_lines, radius_labels )
        ax_h[1][1].legend()

        # Turn on minor ticks for the twins to match what we do for the other
        # axes below.
        ax_h_twin_radius.minorticks_on()
        ax_h_twin_temperature.minorticks_on()

    # Plot background parameters in the remaining subplots.
    starting_index = 4 if compare_flag else 2
    for index, (label, time_series) in enumerate( background_parameters.items() ):
        axis_row_index    = (starting_index + index) // 2
        axis_column_index = (starting_index + index) % 2
        current_axis      = ax_h[axis_row_index][axis_column_index]

        current_axis.set_title( label )
        #
        # NOTE: We set the line's label though we do not show it as it is the
        #       only thing on the plot.  We do this as a convenience for the
        #       caller should they want to add additional data and turn on the
        #       legend after we return.
        #
        current_axis.plot( times, time_series, label=label )
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
    Wrapper for plot_droplet_size_temperatures() that plots a single particle's
    radius and temperatures provided in a DataFrame.  Time series columns are
    specified by one or more evaluation tags.

    Takes 2 arguments:

      particle_dataframe - Pandas DataFrame-like containing the evaluation timeseries
                           to plot.
      evaluation_tags    - List of evaluation tag strings to visualize, with the
                           first being the reference for comparison if two or
                           more tags are supplied.  As a convenience, may be
                           specified as a scalar string which is treated as if
                           it was a list with one element.
      **kwargs           - Optional keyword arguments to pass to plot_droplet_size_temperatures().

    Returns 2 values:

      fig_h - Figure generated for the plots.  See
              plot_droplet_size_temperatures() for details.
      ax_h  - Array of axes that were used for plotting.  See
              plot_droplet_size_temperatures() for details.

    """

    # If just one tag was provided, wrap it into an array
    if isinstance( evaluation_tags, str ):
        evaluation_tags = [evaluation_tags]

    # Add a default title to the figure if the user did not provide one.
    if "title_string" not in kwargs:
        evaluation_string = ", ".join( [evaluation_tag for evaluation_tag in evaluation_tags] )
        particle_id       = particle_dataframe.name
        kwargs["title_string"] = ("Droplet size/temperatures for {:s}\n" +
                                  "Particle {:d}").format( evaluation_string, particle_id )

    # Get the data to plot as series.
    times = insert_timeseries_gaps( particle_dataframe["times"],
                                    particle_dataframe["gap indices"] )

    size_temperatures = {}
    for evaluation_tag in evaluation_tags:
        # Get this evaluation's radii and temperatures column names.
        column_names = get_evaluation_column_names( evaluation_tag )

        # Convert the columns into a two column array.
        #
        # NOTE: We have to stack columns to get a NumPy array as the DataFrame's
        #       .values will return a list of 1D NumPy arrays.
        #
        size_temperatures[evaluation_tag] = insert_timeseries_gaps( np.column_stack( particle_dataframe[column_names].values ),
                                                                    particle_dataframe["gap indices"] )

    # Add gaps into the background parameters if they were provided.  We take
    # care to work with a shallow copy of them so we don't modify the caller's
    # version.
    if "background_parameters" in kwargs:
        background_parameters = kwargs["background_parameters"].copy()

        for parameter_name, parameter_values in background_parameters.items():
            background_parameters[parameter_name] = insert_timeseries_gaps( parameter_values,
                                                                            particle_dataframe["gap indices"] )

        kwargs["background_parameters"] = background_parameters

    fig_h, ax_h = plot_droplet_size_temperatures( times, size_temperatures, **kwargs )

    return fig_h, ax_h

def plot_droplet_size_temperatures_domain( input_parameters, model=None, dt=None, final_time=10.0, **kwargs ):
    """
    Wrapper for plot_droplet_size_temperatures() that evaluates the BDF for the
    supplied radii and temperatures time series with the given background
    conditions.  Solutions are evaluated from zero until a fixed end time.  When
    a MLP model is provided then its evaluation for the same radii,
    temperatures, and backgrounds is plotted for comparison.

    Takes 5 arguments:

      input_parameters - 1D array, of length 6, containing the initial droplet parameters.
      model            - Optional PyTorch model to compare with BDF.  If omitted,
                         defaults to None, the MLP evaluation is skipped, and only
                         the BDF evaluation is visualized.
      dt               - Optional float determining the time step for model
                         evaluation.  If omitted, defaults to the 10 to the mean
                         of the log time range's exponent.  Must be positive.
      final_time       - Optional positive float specifying the end of time
                         range to evaluate.  If omitted, defaults to 10 seconds.
                         Must be no smaller than dt.
      **kwargs         - Optional keyword arguments to pass to plot_droplet_size_temperatures().

    Returns 2 values:

      fig_h - Figure generated for the plots.  See
              plot_droplet_size_temperatures() for details.
      ax_h  - Array of axes that were used for plotting.  See
              plot_droplet_size_temperatures() for details.

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
                                               np.full( (NUMBER_TIME_POINTS, 1),
                                                        dt,
                                                        dtype=np.float32 ),
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
    Wrapper for plot_droplet_size_temperature_dataframe() that plots the
    a single particle's differences between a reference and an evaluation, and
    overlays a scoring report's deviations.  The supplied DataFrame must contain
    columns for the evaluation tags present in the scoring report.

    Takes 3 Arguments:

      particle_dataframe - Pandas DataFrame-like containing the evaluation timeseries
                           to plot.
      score_report       - ScoringReport object to plot deviations from.  The
                           .reference_evaluation and .comparison_evaluation tags
                           must exist in particles_dataframe.
      **kwargs           - Optional keyword arguments to pass to
                           plot_droplet_size_temperatures_dataframe().

    Returns 2 Values:

      fig_h - Figure generated for the plots.  See
              plot_droplet_size_temperatures() for details.
      ax_h  - Array of axes that were used for plotting.  See
              plot_droplet_size_temperatures() for details.

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
      6. Salt solute, in kilograms
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

        # Identify the observations of interest.
        times_mask = (particle["times"] >= time_range[0]) & (particle["times"] <= time_range[1])

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

        # Plot our quantities with labels.  We insert NaNs where gaps are
        # located so quantities plotted are segmented and do not span gaps.
        #
        # NOTE: We convert a single quantity at a time to reduce memory usage.
        #
        particle_times    = insert_timeseries_gaps( particle["times"][times_mask],
                                                    particle["gap indices"] )
        particle_quantity = insert_timeseries_gaps( particle["integration times"][times_mask],
                                                    particle["gap indices"] )

        ax_h[0].plot( particle_times,
                      particle_quantity,
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
            particle_quantity = insert_timeseries_gaps( particle[evaluation_radii_name][times_mask],
                                                        particle["gap indices"] )
            ax_h[1].plot( particle_times,
                          particle_quantity,
                          color=colors_map[evaluation_label],
                          label=evaluation_label )

            particle_quantity = insert_timeseries_gaps( particle[evaluation_temperatures_name][times_mask],
                                                        particle["gap indices"] )
            ax_h[2].plot( particle_times,
                          particle_quantity,
                          color=colors_map[evaluation_label],
                          label=evaluation_label )

        ax_h[1].set_title( "Output Radius" )
        ax_h[1].set_ylabel( "Size (m)" )

        ax_h[2].set_title( "Particle Temperature" )
        ax_h[2].set_ylabel( "Kelvin" )

        # Now plot the environments for the particle.
        particle_quantity = insert_timeseries_gaps( particle["air temperatures"][times_mask],
                                                    particle["gap indices"] )
        ax_h[3].plot( particle_times,
                      particle_quantity,
                      color=colors_map[particle_label],
                      label=particle_label )
        ax_h[3].set_title( "Air Temperature" )
        ax_h[3].set_ylabel( "Kelvin" )

        particle_quantity = insert_timeseries_gaps( particle["relative humidities"][times_mask],
                                                    particle["gap indices"] )
        ax_h[4].plot( particle_times,
                      particle_quantity * 100,
                      color=colors_map[particle_label],
                      label=particle_label )
        ax_h[4].set_title( "Relative Humidity" )
        ax_h[4].set_ylabel( "Percentage (%)" )

        particle_quantity = insert_timeseries_gaps( particle["salt solutes"][times_mask],
                                                    particle["gap indices"] )
        ax_h[5].plot( particle_times,
                      particle_quantity,
                      color=colors_map[particle_label],
                      label=particle_label )
        ax_h[5].set_title( "Salt Solute" )
        ax_h[5].set_ylabel( "kg" )

        particle_quantity = insert_timeseries_gaps( particle["air densities"][times_mask],
                                                    particle["gap indices"] )
        ax_h[6].plot( particle_times,
                      particle_quantity,
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
