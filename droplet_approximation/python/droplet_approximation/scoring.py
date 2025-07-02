from enum import Enum

import matplotlib.pyplot as plt

import numpy as np

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import GridSearchCV
from sklearn.mixture import GaussianMixture

from .data import be_success_mask
from .models import do_inference, do_iterative_inference

def calculate_cusum( differences, tolerance ):
    """
    Calculates the cumulative sum (CUSUM) on an array. Can be 1D for
    one variable or 2D for calculation of multiple variables at once.
    Usually used on an arrays of differences in MLP and BE output.

    Takes 2 arguments:
 
      rt_differences  - Array, 1D for one variable or 2D for simultaneous
                        evaluation, sized number_variables x number_observations.
                        Contains arrays of the differences to analyze.
      tolerance       - Array or Float, contains the error tolerances for each
                        variable. Corresponds to `k` in CUSUM formula.

    Returns 1 value:
      cusum           - Array, sized 2 x number_variables x number_observations,
                        containing:
                            an array of the positive cumulative sums for each variable
                                at each time step at index 0
                            an array of the negative cumulative sums for each variable
                                at each time step at index 1

    """
    cusum = np.zeros( shape=( *differences.shape[:-1], 2, differences.shape[-1] ) )

    for i in range( 1, differences.shape[-1] ):
        cusum[..., 0, i] = np.maximum( cusum[..., 0, i-1] + differences[..., i-1] - tolerance, 0 )
        cusum[..., 1, i] = np.minimum( cusum[..., 1, i-1] + differences[..., i-1] + tolerance, 0 )

    return cusum

def calculate_nrmse( truth_output, model_output ):
    """
    Calculates the normalized root-mean squared error (NRMSE) between the
    provided truth and model outputs. Returns the NRMSE across all observations
    (rows). NRMSE is normalized on the average of the truth value.

    Takes 3 arguments:

      truth_output      - NumPy array, sized number_observations x 2, containing
                          the truth values for radii and temperatures.
      model_output      - NumPy array, sized number_observations x 2, containing
                          the model values for radii and temparatures.

    Returns 2 values:

      radii_nrmse        - NRMSE of the radii.
      temperatures_nrmse - NRMSE of the temperatures.

    """

    RMSE = np.sqrt( np.mean( ( truth_output - model_output ) ** 2, axis=0 ) )
    truth_averages = np.abs( np.mean( truth_output, axis=0 ) )

    return ( RMSE / truth_averages ).sum()

def detect_cusum_deviations( cusum_data, threshold ):
    """
    Creates a mask the size of the incoming CUSUM data that identifies
    where the CUSUM first exceeds the threshold.

    Takes 2 arguments:
      cusum_data      - Array, sized number_variables x 2 x number_observations. 
                        Optionally can provide a 2D array of 2 x number_observations
                        to detect deviations for just one variable.
      threshold       - Array, sized number_variables or just one float if evaluating
                        on one variable. Contains desired threshold for each CUSUM array 
                        in the last dimension of `cusum_data`. Corresponds to
                        `h` in the CUSUM formula.

    Returns 1 value:
      cusum_edge_mask - Array, sized ... x number_observations. Contains
                        a mask with `True` at the index where the CUSUM
                        exceeded `threshold` and `False` everywhere else.
    """

    if cusum_data.ndim == 2:
        threshold = threshold[:, None]
    if cusum_data.ndim == 3:
        threshold = threshold[:, None, None]

    cusum_threshold = np.abs( cusum_data ) > threshold
    cusum_edge_mask = ( cusum_threshold[..., 1:] == 1 ) & ( cusum_threshold[..., :-1] == 0 )
    cusum_edge_mask = np.insert(cusum_edge_mask, 0, False, axis=-1 )

    return cusum_edge_mask

class DeviationDirection( Enum ):
    """
    Encodes the direction of each deviation
    """
    NEGATIVE = 0
    POSITIVE = 1

class DeviationParameter( Enum ):
    """
    Encodes the parameter associated with each deviation
    """
    RADIUS = 0
    TEMPERATURE = 1

class ScoringReport():
    """
    A class to handle the information when scoring a model.  

    Properties:
      cluster_centers             - NumPy array, sized number_clusters
                                    times 7, containing the droplet
                                    parameters for the centroid of each
                                    cluster
      cluster_model               - The Gaussian mixture model that was
                                    used to cluster the data
      cusum_error_threshold       - NumPy array, sized 2, containing the
                                    threshold used to calculate the CUSUM
                                    for radius and temperature respectively.
      cusum_error_tolerance       - NumPy array, sized 2, containing the
                                    tolerance used to calculate the CUSUM
                                    for radius and temperature respectively.
      deviations                  - NumPy array, sized number_deviations x 7
                                    containing the droplet parameters for
                                    each deviation
      deviation_clusters          - NumPy array of integers, size number_deviations,
                                    containing the index of the cluster that the
                                    deviation belongs to.
      deviation_directions        - NumPy array of DeviationDirections, sized
                                    number_deviations containing whether each
                                    deviation was positive or negative.
      deviation_particle_ids      - NumPy array of integers containing the particle
                                    id of each deviation.
      deviation_parameters        - NumPy array of DeviationParameters, sized
                                    number_deviations containing whether each
                                    deviation occured in the radius or temperature
                                    output.
      deviation_times             - NumPy array of floats, sized number_deviations,
                                    containing the time at which each deviation
                                    occured.
      iterative                   - Boolean, False if the scoring was done
                                    with direct model outputs. True if the
                                    scoring was done with iterative outputs.
      model_name                  - String, the name of the model evaluated.
      net_nrmse                   - Float, NRMSE of the entire dataset
      pca_model                   - Scikit-learn PCA, PCA object for deviation parameters
      per_particle_nrmse          - Dictionary, keys are particle IDs, values are floats
                                    corresponding to the model's NRMSE on the
                                    given particle.
      z_scorer                    - Scikit-learn StandardScaler, z-scores deviation parameters
      zs_mean                     - The mean of each droplet parameter pre-zscore. Used to convert
                                    between PCA coordinates and natural ranges.

    """
    def __init__( self, droplet_parameter_list, particle_ids, model, model_name, device,
                 cusum_error_tolerance, cusum_error_threshold, norm=None, iterative=False, 
                 pca_threshold = 0.90, max_clusters = 7 ):
        """
        Takes 12 arguments:
          droplet_parameter_list      - List of NumPy arrays containing the droplet
                                        parameters for each particle.
          particle_ids                - List containing the particle ids to identify
                                        the elements in droplet_parameter_list.
          model_name                  - The name of the model evalutes.
          device                      - The device to evaluate on.
          cusum_error_tolerance       - NumPy Array, sized 2, containing the radius/tempearture
                                        tolerances to pass to `calculate_cusum`.
          cusum_error_threshold       - NumPy Array, sized 2, containing the radius/temperature
                                        thresholds to pass to  `detect_cusum_deviations`
          norm                        - User provided norm to apply to all data before
                                        scoring and deviation detection. Accepts an array
                                        sized number_observations x 2 and returns a normed
                                        array of the same length. Defaults to `identity_norm`.
          iterative                   - Boolean, indicates whether to evaluate the model directly
                                        on every time step, or to iterate radius and temperature
                                        from the first time step.
          pca_threshold               - Float between 0 and 1, determines the minimum data
                                        retention PCA must attain. Defaults to 90%.
          max_clusters                - Integer, identifies the max number of clusters for
                                        the Gaussian mixture model to attempt. Defaults to 7.
        """

        if norm is None:
            norm = identity_norm

        self.cusum_error_tolerance = cusum_error_tolerance
        self.cusum_error_threshold = cusum_error_threshold

        self.iterative = iterative
        self.model_name = model_name

        particle_count = len( droplet_parameter_list )

        # Clean bad BE trials from data
        for particle_index in range( particle_count ):
            droplet_parameter_list[particle_index] = droplet_parameter_list[particle_index][
                be_success_mask( droplet_parameter_list[particle_index][:, 0] ) ]

        particle_times_list = [ np.cumsum( droplet_parameter_list[particle_index][:, -1] )
                                for particle_index in range( particle_count ) ]

        # Calculate MLP results by particle
        # Extract the be outputs from the be inputs
        # Ignoring the first input in non-iterative case
        # Since it isn't an output of do_inference
        if iterative:
            normed_be_output_list = [ norm( particle_parameters[:, :2] )
                                      for particle_parameters in droplet_parameter_list ]
        else:
            normed_be_output_list = [ norm( particle_parameters[1:, :2] ) 
                                      for particle_parameters in droplet_parameter_list ]

        # Generate list of normed MLP output parameters
        normed_mlp_output_list = [ np.empty( shape=( particle_parameters.shape[0] - 1, 2) ) 
                                   for particle_parameters in droplet_parameter_list ]
        for particle_index, particle_parameters in enumerate( droplet_parameter_list ):
            # It feels a little silly to run do_inference and seperate the parameters
            # And then stack them again in `do_inference...`
            if iterative:
                normed_mlp_output_list[particle_index] = norm( do_iterative_inference( 
                                                particle_parameters[:, :-1],
                                                np.cumsum( particle_parameters[:, -1] ),
                                                model,
                                                device) )
            else:
                normed_mlp_output_list[particle_index] = norm( do_inference( 
                                                particle_parameters[:, :-1],
                                                particle_parameters[:, -1],
                                                model,
                                                device)[:-1] )

        # Calculate overall NRMSE
        merged_normed_be_outputs = np.vstack( normed_be_output_list )
        merged_normed_mlp_ouputs = np.vstack( normed_mlp_output_list )

        self.net_nrmse = calculate_nrmse( merged_normed_be_outputs, merged_normed_mlp_ouputs )

        del merged_normed_be_outputs # Delete these large arrays since they are no longe rneeded
        del merged_normed_mlp_ouputs

        # Per particle NRMSE
        self.per_particle_nrmse =   { particle_ids[i]: calculate_nrmse(
                                        normed_be_output_list[i],
                                        normed_mlp_output_list[i] )
                                    for i in range(particle_count) }

        # Per particle CUSUM detection
        deviations = []
        deviation_directions = []
        deviation_particle_ids = []
        deviation_parameters = []
        deviation_times = []

        # Write one deviations array and store seperately the metadata in some other array
        for particle_index in range( particle_count ):
            # Create an array of the radius differences and temperature differences
            normed_output_differences = ( normed_mlp_output_list[particle_index]
                                        - normed_be_output_list[particle_index] ).T

            # Calculate CUSUM. There structure of this array will be
            # [ [ Positive Radius CUSUM, Negative Radius CUSUM ],
            #   [ Positive Temperature CUSUM, Negative Temperature CUSUM ] ]
            cusum = calculate_cusum( normed_output_differences, cusum_error_tolerance )

            # Create a mask for when positive/negative radius/temperature error appear
            # The structure of this array will be
            # [ Positive Radius Deviation Mask, Negative Radius Deviation Mask,
            #   Positive Temperature Deviation Mask, Negative Temperature Deviation Mask ]
            deviation_masks =  detect_cusum_deviations( cusum, cusum_error_threshold ).reshape( (4, -1) )

            # Encode the direction/parameter for each mask
            deviation_direction_vector = np.array( [ DeviationDirection.POSITIVE,
                                                    DeviationDirection.NEGATIVE,
                                                     DeviationDirection.POSITIVE,
                                                     DeviationDirection.NEGATIVE ] )
            deviation_parameter_vector = np.array( [ DeviationParameter.RADIUS,
                                                     DeviationParameter.RADIUS,
                                                     DeviationParameter.TEMPERATURE, 
                                                     DeviationParameter.TEMPERATURE ] )

            # Correct for different output lengths of inference and standard inference
            if not iterative:
                deviation_masks = np.insert( deviation_masks, 0, False, axis=-1 )

            # Count the different types of deviations
            deviation_counts = deviation_masks.sum( axis=1 )

            # Add the droplet parameters to the deviation array
            deviations.extend( np.vstack(
                [ droplet_parameter_list[particle_index][mask] for mask in deviation_masks ] ) )

            # Add the particle ids to the deviation particle id array
            deviation_particle_ids.extend( np.full( deviation_counts.sum(),
                                                    particle_ids[particle_index] ) )

            # Add the deviation directions to the deviation direction array
            deviation_directions.extend( np.hstack( [ np.full( count, direction ) for count, direction
                                            in zip( deviation_counts, deviation_direction_vector ) ] ) )

            # Add the deviation parameters to the deviation direction array
            deviation_parameters.extend( np.hstack( [ np.full( count, direction ) for count, direction
                                            in zip( deviation_counts, deviation_parameter_vector ) ] ) )

            # Add the deviation times to the deviation time array
            deviation_times.extend( np.hstack(
                [ particle_times_list[particle_index][mask] for mask in deviation_masks ] ) )

        # Put the data into np arrays for easier manupulation
        self.deviations = np.array( deviations )
        self.deviation_directions = np.array( deviation_directions )
        self.deviation_particle_ids = np.array( deviation_particle_ids )
        self.deviation_parameters = np.array( deviation_parameters )
        self.deviation_times = np.array( deviation_times )

        # Clustering
        self.z_scaler = StandardScaler()
        z_scored_deviations = self.z_scaler.fit_transform( self.deviations )
        self.zs_mean = np.mean( z_scored_deviations, axis=0 ) # Used later to reverse PCA

        self.pca_model = PCA( n_components = pca_threshold )
        deviations_pca = self.pca_model.fit_transform( z_scored_deviations )

        # Determine optimal gaussian mixture cluster parameters
        def gmm_bic_score( estimator, X ):
            return -estimator.bic(X)

        param_grid = {
            "n_components": range(1, max_clusters)
        }

        grid_search = GridSearchCV(
            GaussianMixture( init_params="k-means++", n_init=15, covariance_type="full" ),
                param_grid=param_grid, scoring=gmm_bic_score
        )

        grid_search.fit( deviations_pca )

        # Cluster based on the best BIC score
        self.cluster_model = grid_search.best_estimator_
        self.deviation_clusters = self.cluster_model.predict( deviations_pca )

        self.cluster_centers = self.z_scaler.inverse_transform(
            np.dot( self.cluster_model.means_, self.pca_model.components_ )
          + self.zs_mean )

    def plot_deviations( self ):
        """
        Plots the deviations in their clusters in 3d coordinates:
            x: log10 radius
            y: particle temperature - air temperature
            z: Relative Humidity 

        Likewise labels the coordinates of cluster centers

        Takes no arguments.

        Returns 2 values:
            fig   - the matplotlib figure of the graph
            ax    - the axis of the graph
        """
        plot_data = np.array( [
            np.log10( self.deviations[:, 0] ),
            self.deviations[:, 1] - self.deviations[:, 3],
            self.deviations[:, 4] ] )

        cluster_coordinates = np.array( [
            np.log10( self.cluster_centers[:, 0] ),
            self.cluster_centers[:, 1] - self.cluster_centers[:, 3],
            self.cluster_centers[:, 4] ] )

        cluster_center_labels = [f"Cluster {i}: {cluster_coordinates.T[i]}"
                                 for i in range( self.cluster_model.n_components ) ]

        fig = plt.figure()
        ax = plt.axes(projection="3d")
        ax.scatter(plot_data[0], plot_data[1], plot_data[2],
            c=self.deviation_clusters )

        # Label the cluster centers
        ax.scatter(cluster_coordinates[0], cluster_coordinates[1], cluster_coordinates[2],
            c=np.arange( self.cluster_model.n_components ), s=20 )

        for i in range( self.cluster_model.n_components ):
            ax.text( cluster_coordinates[0][i], cluster_coordinates[1][i],
                cluster_coordinates[2][i], cluster_center_labels[i],
                fontweight="bold", color="r" )

        ax.set_xlabel("Log Radius (m)")
        ax.set_ylabel("Temperature Difference (K)")
        ax.set_zlabel("Relative Humidity (%)")

        ax.set_title(f"{'Iterative' if self.iterative else ''} Deviation Clusters for {self.model_name}\n"
                + f"Radius Tolerance: {100*self.cusum_error_tolerance[0]:.3f}%,"
                + f" Temperature Tolerance: {self.cusum_error_tolerance[1]:.3f} degrees\n"
                + f"Radius Threshold: {self.cusum_error_threshold[0]:.3f},"
                + f" Temperature Threshold: {self.cusum_error_threshold[1]:.3f}\n")

        ax.set_box_aspect(None, zoom=0.85)

        fig.tight_layout()

        return fig, ax

def identity_norm( rt_data ):
    """
    A blank norm to use as a placeholder in analysis functions. Does
    nothing to `rt_data`.

    Takes 1 argument:
      rt_data     - NumPy array, time steps x 2 containing
                        radius data in natural ranges at 0
                        temperature data in natural ranges at 1

    Returns 1 value:
      rt_data     - NumPy array, time steps x 2 containing
                        radius data in natural ranges at 0
                        temperature data in natural ranges at 1

    """
    return rt_data

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

def standard_norm( rt_data ):
    """
    Normalizes rt_data by applying log10 to radius.

    Takes 1 argument:
      rt_data     - NumPy array, time steps x 2 containing
                        radius data in natural ranges at 0
                        temperature data in natural ranges at 1

    Returns 1 value:
      normalized_rt_data  - NumPy array, time steps x 2 containing
                                log10 radius data in natural ranges at index 0
                                temperature data in natural ranges at index 1
    """
    normalized_rt_data = np.empty_like( rt_data )
    normalized_rt_data[:, 0] = np.log10( rt_data[:, 0] )
    normalized_rt_data[:, 1] = rt_data[:, 1]

    return normalized_rt_data
