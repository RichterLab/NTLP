from dataclasses import dataclass
from enum import Enum

import matplotlib.pyplot as plt
from matplotlib import colors

import numpy as np

from sklearn.preprocessing import StandardScaler
from sklearn.mixture import BayesianGaussianMixture

from torch import multiprocessing

from .data import be_success_mask, read_particles_data
from .models import do_inference, do_iterative_inference
from .physics import set_parameter_ranges

def calculate_cusum( differences, tolerance ):
    """
    Calculates the cumulative sum (CUSUM) of an array of differences.
    Can be 1D for one variable or 2D for multiple variables at once.
    Usually used on an arrays of radius/temperature differences in
    MLP and BE output.

    Takes 2 arguments:
      differences - NumPy array, 1D for one variable or 2D for simultaneous
                    evaluation, sized number_variables x number_observations,
                    contains arrays of the differences to analyze.
      tolerance   - NumPy array or Float, contains the error tolerances for each
                    variable and corresponds to `k` in CUSUM formula.

    Returns 1 value:
      cusum - NumPy array, sized number_variables x 2 x number_observations,
              contains:
                  an array of the positive cumulative sums for each
                      variable at each time step at index [..., 0, :]
                  an array of the negative cumulative sums for each
                      variable at each time step at index [..., 1, :]

    """

    cusum            = np.zeros( shape=(*differences.shape[:-1], 2, differences.shape[-1]) )
    cusum[..., 0, 0] = np.maximum( differences[..., 0] - tolerance, 0 )
    cusum[..., 1, 0] = np.minimum( differences[..., 0] + tolerance, 0 )

    for time_step in range( 1, differences.shape[-1] ):
        # Positive CUSUM
        cusum[..., 0, time_step] = np.maximum( cusum[..., 0, time_step-1] +
                                               differences[..., time_step] -
                                               tolerance,
                                               0 )
        # Negative CUSUM
        cusum[..., 1, time_step] = np.minimum( cusum[..., 1, time_step-1] +
                                               differences[..., time_step] +
                                               tolerance,
                                               0 )
    return cusum

def calculate_nrmse( truth_output, model_output ):
    """
    Calculates the normalized root-mean squared error (NRMSE) between the
    provided truth and model outputs. Returns the NRMSE across all observations
    (rows). NRMSE is normalized on the average of the truth value.

    Takes 2 arguments:

      truth_output - NumPy array, sized number_observations x 2, containing
                     the truth values for droplet radii and temperatures.
      model_output - NumPy array, sized number_observations x 2, containing
                     the model values for droplet radii and temperatures.

    Returns 1 value:

      NRMSE - the normalized root mean squared error of model_output against
              truth_output.
    """

    RMSE           = np.sqrt( np.mean( (truth_output - model_output)**2, axis=0 ) )
    truth_averages = np.abs( np.mean( truth_output, axis=0 ) )

    return (RMSE / truth_averages).sum()

def detect_cusum_deviations( cusum_data, threshold ):
    """
    Creates a mask the shape of the incoming CUSUM data that identifies
    the indexes where the CUSUM begins to exceed the error threshold (h).

    Takes 2 arguments:
      cusum_data - Array, sized number_variables x 2 x number_observations.
                   Usually, this is output data from `calculate_cusum`.
                   Optionally can provide a 2D array of 2 x number_observations
                   to detect deviations for just one variable. In either case,
                   the second dimension of the array corresponds to
                   positive/negative CUSUMs.
      threshold  - Array, sized number_variables, or a float if evaluating
                   on one variable, containing the threshold to apply to the
                   CUSUM array for each variable respectively. Corresponds to
                   `h` in the CUSUM formula.

    Returns 1 value:
      cusum_edge_mask - Array, sized number_variables x 2 x number_observations,
                        containing an array of masks for each variable and
                        positive/negative CUSUM. False everywhere except
                        at the index where the CUSUM starts exceeding `threshold`.
    """

    if cusum_data.ndim == 3:
        threshold = threshold[:, None, None]

    # If the threshold is exceeded in the first observation, record that
    # as a deviation (True). Then, record all observations# where CUSUM exceeds
    # the threshold but was within the threshold
    # in the previous observation.
    cusum_edge_mask          = np.zeros( cusum_data.shape, dtype=bool )
    cusum_threshold          = np.abs( cusum_data ) > threshold
    cusum_edge_mask[..., 0]  = cusum_threshold[..., 0]
    cusum_edge_mask[..., 1:] = cusum_threshold[..., 1:] & (~cusum_threshold[..., :-1])

    return cusum_edge_mask

class DeviationDirection( Enum ):
    """
    This enum records the direction of a deviation.
    """

    NEGATIVE = 0
    POSITIVE = 1

class DeviationParameter( Enum ):
    """
    This enum records the parameter associated with a deviation.
    """

    RADIUS      = 0
    TEMPERATURE = 1

def identity_norm( rt_data ):
    """
    A blank norm to use as a placeholder in analysis functions. Does
    nothing to `rt_data`.

    Takes 1 argument:
      rt_data - NumPy array, time steps x 2 containing
                    radius data in natural ranges at 0
                    temperature data in natural ranges at 1

    Returns 1 value:
      rt_data - NumPy array, time steps x 2 containing
                    radius data in normed ranges at 0
                    temperature data in normed ranges at 1

    """

    return rt_data


# This may be rewritten in the future as `ParticleScores` in
# order to be more data-oriented. Instead of containing values
# for an individual particle, each object would contain an array
# of values. `particle_pipeline` would return one `ParticleScores`
# object, and `Scoring_Report` `__init__` would concatenate all of
# these arrays together.
@dataclass
class ParticleScore:
    """
    This is a data class to structure output from `particle_pipeline`.

    Contains 10 attributes:

      particle_id            - Integer, id of the particle processed.
      particle_nrmse         - Float, NRMSE of the particle against BE outputs.
      square_error           - NumPy Array, sized 2, contains the square of the
                               difference between BE and MLP output for radius
                               and temperature respectively. Used to calculate
                               overall NRMSE.
      truth_sum              - NumPy Array, sized 2, contains the sum of the truth
                               values for radius and temperature respectively. Used
                               to calculate overall NRMSE.
      observation_count      - Integer, number of observations in this particle's
                               trajectory. Used to calculate overall NRMSE.
      deviations             - NumPy array, sized number_deviations x 7,
                               contains the droplet parameters each deviation,
                               including integration time.
      deviation_directions   - NumPy array of DeviationDirections, sized
                               number_deviations, containing whether each
                               deviation was positive or negative.
      deviation_parameters   - NumPy array of DeviationParameters, sized
                               number_deviations, containing whether each
                               deviation occurred in the radius or temperature
                               output.
      deviation_times        - NumPy array of floats, sized number_deviations,
                               containing the time at which each deviation occurred.
      deviation_particle_ids - NumPy array of integers containing the particle id
                               of each deviation. Technically, this is redundant to
                               particle_id, but it is used by ScoringReport to
                               cleanly generate the overall deviation_particle_ids array.

    """

    particle_id:            int
    particle_nrmse:         np.float64
    square_error:           np.ndarray
    truth_sum:              np.ndarray
    observation_count:      int
    deviations:             np.ndarray
    deviation_directions:   np.ndarray
    deviation_parameters:   np.ndarray
    deviation_times:        np.ndarray
    deviation_particle_ids: np.ndarray

def particle_scoring_pipeline( particles_root, particle_ids, dirs_per_level, model, device, cusum_error_tolerance,
                               cusum_error_threshold, norm, iterative, parameter_ranges ):
    """
    Calculates NRMSE and deviations for a set of particles on a given model. Loads
    particle data from per particle data frames. Also provides statistics needed
    to calcualte overall dataset NRMSE.

    Takes 10 arguments:

      particles_root        - Path to the top-level directory to read raw
                              particle files from.
      particle_ids          - Sequence of particle identifiers to process.
      dirs_per_level        - Number of subdirectories per level in the
                              particle files directory hierarchy.
      model                 - PyTorch model to use
      device                - The device to evaluate on.
      cusum_error_tolerance - NumPy Array, sized 2, containing the
                              radius/temperature tolerances to pass to
                              `calculate_cusum`.
      cusum_error_threshold - NumPy Array, sized 2, containing the
                              radius/temperature thresholds to pass to
                              `detect_cusum_deviations`.
      norm                  - User provided norm to apply to all data before
                              scoring and deviation detection. Accepts an
                              array sized number_observations x 2 and returns
                              a normed array of the same length.
      iterative             - Boolean, indicates whether to evaluate the model
                              directly on every time step, or to iterate radius
                              and temperature from the first time step.
      parameter_ranges      - Dictionary of any changes to the droplet parameter
                              ranges required for the given model.

    Returns 1 Value:

      ParticleScores - List of ParticleScore objects.

    """

    set_parameter_ranges( parameter_ranges )

    particle_scores = []

    for particle_ids_batch in particle_ids:

        particles_df = read_particles_data( particles_root, particle_ids_batch, dirs_per_level )

        for particle_index, particle_id in enumerate( particle_ids_batch ):
            # Load parameters and mask out BE failures
            p_df                = particles_df.iloc[particle_index]
            particle_parameters = np.stack( p_df[[
                                      "input radii",
                                      "input temperatures",
                                      "salt masses",
                                      "air temperatures",
                                      "relative humidities",
                                      "air densities",
                                      "integration times"
                                  ]].to_numpy(), axis=-1 )
            # Calculate actual particle times in simulation time and particle times after
            # removing BE failures. The particle times in normal simulation time are used
            # for reporting where the deviations are. The part times after removing BE failures
            # are used for iterative model inference.
            be_mask               = be_success_mask( particle_parameters[:, 0] )
            actual_particle_times = np.delete( np.cumsum( np.insert( particle_parameters[:, -1],
                                                                     0,
                                                                     0.0 )[:-1] ),
                                               ~be_mask ) + p_df["birth time"]
            particle_parameters   = particle_parameters[be_mask, :]
            particle_times        = np.cumsum( np.insert( particle_parameters[:, -1],
                                                          0,
                                                          0.0 )[:-1] ) + p_df["birth time"]

            # If there is only one trial for a given particle,
            # there will be no MLP outputs to analyze, so we should skip
            if particle_parameters.shape[0] == 1:
                continue

            # Ignore the first trial because it is not an output
            normed_be_output  = norm( particle_parameters[1:, :2] )
            normed_mlp_output = np.empty( shape=(particle_parameters.shape[0] - 1, 2) )

            # Calculate MLP results

            # We remove the last trial of `do_inference` since we don't have a reference
            # radius/temperature for the time step after the last one.
            # We also remove the first trial of `do_iterative_inference` since it just returns
            # the first radius/temperature; i.e. it is not a true MLP output.
            # This also aligns the indices of the results.
            if iterative:
                normed_mlp_output = norm( do_iterative_inference(
                                                            particle_parameters[:, :-1],
                                                            particle_times,
                                                            model,
                                                            device)[1:] )
            else:
                normed_mlp_output = norm( do_inference( particle_parameters[:, :-1],
                                                        particle_parameters[:, -1],
                                                        model,
                                                        device)[:-1] )

            # Calculating the overall NRMSE directly would require copying all of the particle
            # data frames together. Instead, we just copy the square error and the sum of the truth
            # parameters for each particle so that net NRMSE can be calculated later.
            particle_nrmse      = calculate_nrmse( normed_be_output, normed_mlp_output )
            square_error        = np.sum( (normed_be_output - normed_mlp_output)**2, axis=0 )
            truth_sum           = np.abs( np.sum( normed_be_output, axis=0 ) )
            number_observations = particle_parameters.shape[0]

            # Calculate CUSUM. There structure of the result will be
            # [ [ Positive Radius CUSUM array, Negative Radius CUSUM array ],
            #   [ Positive Temperature CUSUM array, Negative Temperature CUSUM array ] ]
            normed_output_differences = (normed_mlp_output - normed_be_output).T
            cusum                     = calculate_cusum( normed_output_differences, cusum_error_tolerance )

            # Create a mask for when positive/negative radius/temperature CUSUM exceeds
            # the threshold. Flatten so that we can zip it into results with the enum arrays.
            # The structure of this array will be
            # [ Positive Radius Deviation Mask, Negative Radius Deviation Mask,
            #   Positive Temperature Deviation Mask, Negative Temperature Deviation Mask ]
            # Also adds a "False" to the start since the first time step isn't an output and so
            # it is not in the CUSUM, but it also never a deviation.
            deviation_masks  = detect_cusum_deviations( cusum, cusum_error_threshold ).reshape( (4, -1) )
            deviation_masks  = np.insert( deviation_masks, 0, False, axis=-1 )
            deviation_counts = deviation_masks.sum( axis=1 )

            # These arrays allow the deviations identified with each
            # CUSUM array to be zipped with their corresponding enums.

            deviation_direction_vector = np.array( [ DeviationDirection.POSITIVE,
                                                     DeviationDirection.NEGATIVE,
                                                     DeviationDirection.POSITIVE,
                                                     DeviationDirection.NEGATIVE ] )
            deviation_parameter_vector = np.array( [ DeviationParameter.RADIUS,
                                                     DeviationParameter.RADIUS,
                                                     DeviationParameter.TEMPERATURE,
                                                     DeviationParameter.TEMPERATURE ] )

            # Record data about the deviations
            deviations             = np.vstack( [ particle_parameters[mask] for mask in deviation_masks ] )
            deviation_directions   = np.hstack( [ np.full( count, direction ) for count, direction
                                                  in zip( deviation_counts, deviation_direction_vector ) ] )
            deviation_parameters   = np.hstack( [ np.full( count, direction ) for count, direction
                                                  in zip( deviation_counts, deviation_parameter_vector ) ] )
            deviation_times        = np.hstack( [ actual_particle_times[mask] for mask in deviation_masks ] )
            deviation_particle_ids = np.full( deviation_counts.sum(), particle_id )

            # In the future, if we want access to CUSUM or MLP radius data, perhaps write it to a new
            # per particle analysis dataframe here

            particle_score = ParticleScore( particle_id, particle_nrmse, square_error, truth_sum,
                                            number_observations, deviations, deviation_directions,
                                            deviation_parameters, deviation_times, deviation_particle_ids )

            particle_scores.append( particle_score )

    return particle_scores

class ScoringReport():
    """
    A class to handle the information when scoring a model.

    Contains 16 attributes:

      cluster_centers        - NumPy array, sized number_clusters
                               times 7, containing the droplet
                               parameters for the centroid of each
                               cluster.
      cluster_model          - The Gaussian mixture model that was
                               used to cluster the data.
      cusum_error_threshold  - NumPy array, sized 2, containing the
                               threshold used to calculate the CUSUM
                               for radius and temperature respectively.
      cusum_error_tolerance  - NumPy array, sized 2, containing the
                               tolerance used to calculate the CUSUM
                               for radius and temperature respectively.
      deviations             - NumPy array, sized number_deviations x 7
                               containing the droplet parameters for
                               each deviation.
      deviation_clusters     - NumPy array of integers, size number_deviations,
                               containing the index of the cluster that the
                               deviation belongs to.
      deviation_directions   - NumPy array of DeviationDirections, sized
                               number_deviations containing whether each
                               deviation was positive or negative.
      deviation_particle_ids - NumPy array of integers containing the particle
                               id of each deviation.
      deviation_parameters   - NumPy array of DeviationParameters, sized
                               number_deviations containing whether each
                               deviation occurred in the radius or temperature
                               output.
      deviation_times        - NumPy array of floats, sized number_deviations,
                               containing the time at which each deviation
                               occurred.
      iterative              - Boolean, False if the scoring was done
                               with direct model outputs. True if the
                               scoring was done with iterative outputs.
      model_name             - String, the name of the model evaluated.
      net_nrmse              - Float, NRMSE of the entire dataset.
      per_particle_nrmse     - Dictionary, keys are particle IDs, values are floats
                               corresponding to the model's NRMSE on the
                               given particle.
      z_score_model          - Scikit-learn StandardScaler, z-scores deviation
                               parameters.
      z_score_mean           - The mean of each droplet parameter pre-zscore.
                               Used to convert between z-scored coordinates and
                               natural ranges.

    """

    def __init__( self, particles_root, particle_ids, dirs_per_level, model, model_name, device,
                  cusum_error_tolerance, cusum_error_threshold, norm=None, iterative=False,
                  max_clusters=7, number_processes=0, number_batches=1, parameter_ranges=None ):
        """
        Takes 12 arguments:
          particles_root        - Path to the top-level directory to read raw
                                  particle files from.
          particle_ids          - Sequence of particle identifiers to process.
          dirs_per_level        - Number of subdirectories per level in the
                                  particle files directory hierarchy.
          model                 - PyTorch model to use
          model_name            - The name of the model being evaluated.
          device                - The device to evaluate on.
          cusum_error_tolerance - NumPy Array, sized 2, containing the
                                  radius/temperature tolerances to pass to
                                  `calculate_cusum`.
          cusum_error_threshold - NumPy Array, sized 2, containing the
                                  radius/temperature thresholds to pass
                                  to `detect_cusum_deviations`.
          norm                  - Optional user provided norm to apply to all
                                  data before scoring and deviation detection.
                                  Accepts an array sized number_observations x 2
                                  and returns a normed array of the same length.
                                  Defaults to `identity_norm`.
          iterative             - Optional Boolean, determines whether the back end
                                  used to calculate MLP output is `do_inference`
                                  (if False) or `do_iterative_inference` (if True).
                                  Defaults to False.
          max_clusters          - Optional Integer, identifies the max number of
                                  clusters for the Bayesian Gaussian mixture model
                                  to attempt. Defaults to 7.
          number_processes      - Optional Integer, the number of workers to
                                  parallelize error analysis across. Defaults
                                  to `multiprocessing.cpu_count()`.
          number_batches        - Optional Integer, determines the number of
                                  batches for each core to break up the particles
                                  into. Higher batches results in fewer particles
                                  loaded at one time. Defaults to 1.
          parameter_ranges      - Optional dictionary of any changes to the droplet
                                  parameter ranges required for the given model.
                                  Defaults to None.
        """

        # The processes must be created with spawn; it seems the pipeline
        # will not parallelize properly with fork.
        multiprocessing.set_start_method( 'spawn', force=True )

        if norm is None:
            norm = identity_norm
        if number_processes == 0:
            number_processes = multiprocessing.cpu_count()

        self.cusum_error_tolerance = cusum_error_tolerance
        self.cusum_error_threshold = cusum_error_threshold
        self.iterative             = iterative
        self.model_name            = model_name

        # Splits all of the requested particle ids across processors
        # and into batches to avoid loading all of the particles
        # at once
        particle_id_chunks  = [np.array_split( process_particle_ids, number_batches )
                               for process_particle_ids in np.array_split( particle_ids,
                                                                           number_processes )]
        pipeline_parameters = [[particles_root, process_particle_ids, dirs_per_level,
                                model, device, cusum_error_tolerance, cusum_error_threshold,
                                norm, iterative, parameter_ranges]
                               for process_particle_ids in particle_id_chunks]

        if number_processes > 1:
            with multiprocessing.Pool( processes=number_processes ) as pool:
                particle_scores_list = pool.starmap( particle_scoring_pipeline, pipeline_parameters )
        else:
            particle_scores_list = [particle_scoring_pipeline( *pipeline_parameters[0] )]

        # Collect results
        self.per_particle_nrmse = {}

        overall_square_error      = np.zeros( 2 )
        overall_truth_sum         = np.zeros( 2 )
        overall_observation_count = 0

        deviations             = []
        deviation_directions   = []
        deviation_parameters   = []
        deviation_times        = []
        deviation_particle_ids = []

        for particle_scores in particle_scores_list:
            for particle_score in particle_scores:
                self.per_particle_nrmse[particle_score.particle_id] = particle_score.particle_nrmse

                overall_square_error      += particle_score.square_error
                overall_truth_sum         += particle_score.truth_sum
                overall_observation_count += particle_score.observation_count

                deviations.append(             particle_score.deviations )
                deviation_directions.append(   particle_score.deviation_directions )
                deviation_parameters.append(   particle_score.deviation_parameters )
                deviation_times.append(        particle_score.deviation_times )
                deviation_particle_ids.append( particle_score.deviation_particle_ids )

        self.net_nrmse = np.sum( np.sqrt( overall_square_error / overall_observation_count ) /
                                 (overall_truth_sum / overall_observation_count) )

        self.deviations             = np.vstack( deviations )
        self.deviation_directions   = np.hstack( deviation_directions )
        self.deviation_parameters   = np.hstack( deviation_parameters )
        self.deviation_times        = np.hstack( deviation_times )
        self.deviation_particle_ids = np.hstack( deviation_particle_ids )

        # Clustering
        self.z_score_model           = StandardScaler()
        z_scored_deviations          = self.z_score_model.fit_transform( self.deviations )
        filtered_z_scored_deviations = z_scored_deviations[np.all( np.abs( z_scored_deviations < 3 ),
                                                                   axis=1)]
        self.z_score_mean            = np.mean( z_scored_deviations, axis=0 )

        self.cluster_model = (BayesianGaussianMixture( init_params="k-means++",
                                                       n_init=5,
                                                       n_components=max_clusters )
                              .fit( filtered_z_scored_deviations ))

        self.deviation_clusters = self.cluster_model.predict( z_scored_deviations )
        self.cluster_centers    = (self.z_score_model.inverse_transform( self.cluster_model.means_ ) +
                                   self.z_score_mean)

    def plot_deviations( self, label_centers=True, thinning_ratio=1 ):
        """
        Plots the deviations in their clusters in 3d coordinates:

          x: log10 radius
          y: particle temperature - air temperature
          z: Relative Humidity

        Likewise labels the coordinates of cluster centers if
        `label_centers=True`.

        Takes 2 argument:

          label_centers  - Optional Boolean, determines whether to label
                           the center of each deviation cluster. Defaults
                           to True.
          thinning_ratio - Optional Integer, determines what fraction of
                           droplets to plot. Defaults to 1.

        Returns 2 values:

          fig - the matplotlib figure of the graph
          ax  - the axis of the graph

        """

        colormap = colors.ListedColormap( ["red",
                                           "blue",
                                           "green",
                                           "orange",
                                           "black",
                                           "yellow",
                                           "teal",
                                           "gold",
                                           "magenta",
                                           "lightgreen",
                                           "navy",
                                           "dimgray",
                                           "lightcoral"] )

        # Gather all of the data in the desired coordinates. Thin the data out based on thinning_ratio
        plot_data             = np.array( [np.log10( self.deviations[::thinning_ratio, 0] ),
                                           self.deviations[::thinning_ratio, 1] - self.deviations[::thinning_ratio, 3],
                                           self.deviations[::thinning_ratio, 4]] )
        cluster_coordinates   = np.array( [np.log10( self.cluster_centers[:, 0] ),
                                           self.cluster_centers[:, 1] - self.cluster_centers[:, 3],
                                           self.cluster_centers[:, 4]] )
        cluster_center_labels = [f"Cluster { cluster_index }: { cluster_coordinates.T[cluster_index] }"
                                 for cluster_index in range( self.cluster_model.n_components )]

        fig = plt.figure()
        ax  = plt.axes( projection="3d" )

        for cluster_id in np.unique( self.deviation_clusters ):
            cluster_plot_data = plot_data[:, self.deviation_clusters[::thinning_ratio] == cluster_id]
            cluster_colors    = np.full( (cluster_plot_data.shape[1], 4), colormap( cluster_id ) )
            ax.scatter( cluster_plot_data[0],
                        cluster_plot_data[1],
                        cluster_plot_data[2],
                        color=cluster_colors,
                        label=f"Cluster {cluster_id}" )

        # Label the cluster centers
        ax.scatter( cluster_coordinates[0],
                    cluster_coordinates[1],
                    cluster_coordinates[2],
                    c=np.arange( self.cluster_model.n_components ),
                    s=20 )

        if label_centers:
            for cluster_index in range( self.cluster_model.n_components ):
                ax.text( cluster_coordinates[0][cluster_index], cluster_coordinates[1][cluster_index],
                         cluster_coordinates[2][cluster_index], cluster_center_labels[cluster_index],
                         fontweight="bold",
                         color="r" )

        ax.set_xlabel( "Log Radius (m)" )
        ax.set_ylabel( "Temperature Difference (K)" )
        ax.set_zlabel( "Relative Humidity (%)" )

        # Do not pad these curly braces with spaces.
        # ':3f ' is not the same as ':3f'. The former
        # will throw an error.
        ax.set_title( f"{'Iterative' if self.iterative else ''} Deviation Clusters " +
                      f"for { self.model_name }\n" +
                      f"Radius Tolerance: {100*self.cusum_error_tolerance[0]:.3f}%," +
                      f" Temperature Tolerance: {self.cusum_error_tolerance[1]:.3f} degrees\n" +
                      f"Radius Threshold: {self.cusum_error_threshold[0]:.3f}," +
                      f" Temperature Threshold: {self.cusum_error_threshold[1]:.3f}\n" )

        # Zooming out prevents clipping the z-axis label

        ax.legend()
        fig.tight_layout()

        return fig, ax

    def recluster( self, z_score_outlier_threshold = 3, n_init = 5, n_components = 7, init_params="k-means++", **kwargs ):
        """
        Reruns deviation clustering with the given arguments.

        Takes 5 arguments:
          z_score_outlier_threshold - Optional integer, determines how many standard deviations
                                      off of average deviation droplet parameters to include in clustering.
                                      Defaults to 3.
          n_init                    - Optional integer, specifies how many times to initialize clustering.
                                      Defaults to 5.
          n_components              - Optional integer, specifies how many clusters to find.
                                      Defaults to 7.
          init_params               - Optional String, specifies what method BayesianGaussianMixture
                                      uses for initialization. Defaults to 'k-means++'.
          **kwargs                  - Additional arguments for BayesianGaussianMixture.

        Returns nothing.
        """

        self.z_score_model           = StandardScaler()
        z_scored_deviations          = self.z_score_model.fit_transform( self.deviations )
        filtered_z_scored_deviations = z_scored_deviations[np.all( np.abs( z_scored_deviations < z_score_outlier_threshold ),
                                                                   axis=1)]
        self.z_score_mean            = np.mean( z_scored_deviations, axis=0 )

        self.cluster_model = ( BayesianGaussianMixture( init_params=init_params,
                                                        n_init=n_init,
                                                        n_components=n_components,
                                                        **kwargs )
                                                       .fit( filtered_z_scored_deviations ) )

        self.deviation_clusters = self.cluster_model.predict( z_scored_deviations )
        self.cluster_centers    = ( self.z_score_model.inverse_transform( self.cluster_model.means_ )
                                  + self.z_score_mean )

def standard_norm( rt_data ):
    """
    Normalizes rt_data by applying log10 to radius.

    Takes 1 argument:

      rt_data - NumPy array, number_observations x 2 containing:

                  radius data in natural ranges at index 0
                  temperature data in natural ranges at index 1

    Returns 1 value:

      normalized_rt_data - NumPy array, time steps x 2 containing:

                             log10 radius data in natural ranges at index 0
                             temperature data in natural ranges at index 1

    """

    normalized_rt_data       = np.empty_like( rt_data )
    normalized_rt_data[:, 0] = np.log10( rt_data[:, 0] )
    normalized_rt_data[:, 1] = rt_data[:, 1]

    return normalized_rt_data
