import os
import subprocess

# Dummy class intended to stand in for the Weights and Biases (W&B) Python
# package, wandb, when it isn't install or the user elects to not use it.
# An instance of this class can be used in lieu of the wandb package and
# immediately returns from every function and method invocation.
#
# None of the methods are documented so as to save maintainer effort and cut
# down on verbosity.  The two dummy classes, NoOpWandB and NoOpRun, each have
# links to the documentation of the API they provide a subset of.
#
# NOTE: This is not a fully featured mock of the wandb package but rather
#       enough of a subset to support what the rest of this module requires.
#
class NoOpWandB:
    """
    No-op dummy class that accepts Weights and Biases (W&B) calls but does
    nothing.

    NOTE: When standing in for W&B's wandb module, this class needs to be
          instantiated instead of simply used in lieu of the wandb module!
          Without this, you'll get the following cryptic error:

          TypeError: NoOpWandB.init() missing 1 required positional argument: 'self'

    See https://docs.wandb.ai/ for more details.
    """

    class NoOpRun:
        """
        No-op dummy for the Run class.  Objects of this class accept a barebones
        subset of wandb.Run's methods and simply returns or ignores.

        See https://docs.wandb.ai/ref/python/public-api/run/ for more details.
        """

        def __init__( self, *args, **kwargs ):
            pass

        def __getattr__( self, name ):
            return self._noop

        def __setattr__( self, name, value ):
            pass

        def _noop( self, *args, **kwargs ):
            return self

        def log( self, *args, **kwargs ):
            pass

        def finish( self, *args, **kwargs ):
            pass

        def watch( self, *args, **kwargs ):
            pass

    class NoOpArtifact:
        """
        No-op dummy for the Artifact class.  Objects of this class accept a
        barebones subset of wandb.Artifact's methods and simply returns or
        ignores.

        See https://docs.wandb.ai/ref/python/experiments/artifact/ for more
        details.
        """

        def __init__( self, *args, **kwargs ):
            pass

        def add_file( self, *args, **kwargs ):
            pass

        def add_reference( self, *args, **kwargs ):
            pass

        def save( self, *args, **kwargs ):
            pass

    # Make our dummy class a class attribute so we can imitate wandb.Artifact().
    Artifact = NoOpArtifact

    def __init__( self ):
        # Expose the current Run object to the world.
        self.run = None

    def init( self, *args, **kwargs ):
        # Track the current "run".
        self.run = self.NoOpRun()

        return self.run

    def log( self, *args, **kwargs ):
        pass

    def logged_artifacts( self, *args, **kwargs ):
        pass

    def finish( self, *args, **kwargs ):
        pass

    def config( self, *args, **kwargs ):
        pass

    def __getattr__( self, name ):
        return lambda *args, **kwargs: None

def _get_repository_information( repo_path="." ):
    """
    Retrieves the Git repository information associated with a directory path.
    Returns the current commit's SHA, the branch name, and a Boolean flag
    indicating whether any tracked files have been modified.

    Takes 1 argument:

      repo_path - Optional path to the Git repository to query.  If omitted,
                  defaults to the current working directory of the process.

    Returns 3 values:

      commit_sha  - String containing the short Git commit SHA if it could be
                    determined, "unknown" otherwise.
      branch_name - String containing the Git branch name if it could be
                    determined, "unknown" otherwise.
      clean_flag  - Boolean flag indicating that none of the track files are
                    modified relative to commit_sha.  False if at least one
                    file is modified or if it could not be determined.

    """

    # We capture the output of each command as text.  Commands that exit with a
    # non-zero status code should raise subprocess.CalledProcessError.
    command_args = {
        "capture_output": True,
        "check":          True,
        "text":           True
    }

    try:
        # Check if we're in a Git repository.  We don't care about the
        # output, merely that we executed correctly.  If we don't, then
        # we're not in a repository and we bail early.
        subprocess.run(
            ["git", "-C", repo_path, "rev-parse", "--git-dir"],
            **command_args )

        # Commit SHA.
        sha_result = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            **command_args )

        # Branch name.
        branch_result = subprocess.run(
            ["git", "rev-parse", "--abbrev-ref", "HEAD"],
            **command_args )

        # Clean state.
        clean_result = subprocess.run(
            ["git", "status", "--porcelain", "--untracked-files=no"],
            **command_args )

    except (subprocess.CalledProcessError, FileNotFoundError):
        return ("unknown", "unknown", False)

    # Convert the results into strings and a Boolean flag.
    commit_sha  = sha_result.stdout.strip()
    branch_name = branch_result.stdout.strip()
    clean_flag  = len( clean_result.stdout.strip() ) == 0

    return commit_sha, branch_name, clean_flag

def prepare_wandb_run( wandb, project_name, model, criterion, optimizer, batch_size, lr_scale, file_paths, number_droplets ):
    """
    Creates a Weights and Biases (W&B) Run object to log checkpoints (artifacts)
    and training metrics (metrics) to.  This captures metadata related to the
    code, model architecture/configuration/hyperparameters, and the
    training/validation datasets.

    Takes 9 arguments:

      wandb           - W&B module handle.
      project_name    - W&B project name.  Must be non-empty.
      model           - Torch model object.
      criterion       - Loss function handle.
      optimizer       - Torch optimizer object.
      batch_size      - Integer batch size used during training.
      lr_scale        - Floating point learning rate scale factor.  The optimizer's
                        learning rate is scaled by this after each epoch.
      file_paths      - Sequence of two string paths specifying the location of
                        training and validation.  The validation data path may
                        be None when validation is not used.
      number_droplets - Sequence of two integers specifying the number of
                        training and validation samples in file_paths.  The
                        number of validation samples must be None when
                        validation is not used.

    Returns 1 value:

      run - W&B Run object.

    """

    # Unpack the training/validation data.
    (training_file,
     validation_file) = file_paths

    (number_training_samples,
     number_validation_samples) = number_droplets

    # Get version information about the code that is running.
    (commit_sha,
     branch_name,
     clean_flag) = _get_repository_information( os.path.dirname( os.path.abspath( __file__ ) ) )
    code_version_str = "Branch '{:s}' (Commit SHA {:s}{:s})".format(
        branch_name,
        commit_sha,
        "" if clean_flag else "*" )

    wandb_config = {
        # Code version.
        "code_version":                 code_version_str,

        # Model architecture.
        "activation_function":          model.activation_func(),
        "architecture":                 model.architecture(),
        "layer_sizes":                  model.layer_sizes(),
        "number_layers":                model.number_layers(),
        "number_parameters":            sum( p.numel() for p in model.parameters() ),

        # Loss function.
        "loss_function":                str( criterion ),
        "batch_size":                   batch_size,

        # Optimizer hyperparameters.
        "optimizer":                    optimizer.__class__,
        "learning_rate":                optimizer.defaults["lr"],
        "learning_rate_scale_factor":   lr_scale,
        "weight_decay":                 optimizer.defaults["weight_decay"],

        # Training and validation data.
        "samples_training":             number_training_samples,
        "samples_validation":           0 if validation_file is None else number_validation_samples,
    }

    # Create separate artifacts for the training and validation files.
    training_dataset_artifact = wandb.Artifact( name="{:s}-training-data".format( model.name() ),
                                                type="dataset" )
    validation_dataset_artifact = wandb.Artifact( name="{:s}-validation-data".format( model.name() ),
                                                  type="dataset" )

    # Record the training file as an external file reference if we have it.
    if isinstance( training_file, str ):
        wandb_config["path_training"] = os.path.abspath( training_file )
        training_dataset_artifact.add_reference( name="training-data",
                                                 uri="file://{:s}".format( wandb_config["path_training"] ) )
    else:
        wandb_config["path_training"] = "N/A"

    # Record the validation file as an external file reference if we have it.
    if isinstance( validation_file, str ):
        wandb_config["path_validation"] = os.path.abspath( validation_file )
        validation_dataset_artifact.add_reference( name="validation-data",
                                                   uri="file://{:s}".format( wandb_config["path_validation"] ) )
    else:
        wandb_config["path_validation"] = "N/A"

    # Create the Run object with all of the hyperparameters and metadata for
    # this training run.
    run = wandb.init( project=project_name,
                      name=model.name(),
                      config=wandb_config )

    # Register this run's use of our dataset artifacts.
    run.use_artifact( training_dataset_artifact )
    run.use_artifact( validation_dataset_artifact )

    # Define our timeline's steps which correspond to batch-specific and
    # epoch-specific metrics.
    run.define_metric( "batch" )
    run.define_metric( "epoch" )

    # Batch-specific metrics.
    run.define_metric( "training_loss", step_metric="batch" )

    # Epoch-specific metrics.
    run.define_metric( "learning_rate",   step_metric="epoch" )
    run.define_metric( "validation_loss", step_metric="epoch" )

    # Periodically record the model's parameters and gradients.  Aim for two
    # sets of gradients and parameters per epoch to avoid excessive logging.
    number_logs   = 2
    log_frequency = (number_training_samples // batch_size) // number_logs
    run.watch( model, log="all", log_freq=log_frequency )

    return run
