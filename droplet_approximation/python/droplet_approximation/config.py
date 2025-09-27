from configparser import ConfigParser, ExtendedInterpolation

from pathlib import Path
from functools import partial
from copy import deepcopy

import numpy as np

config      = None

CONFIG_ROOT = ""
CONFIG_PATH = ""

def _filter( value ):
    """
    Preprocessing for all custom config getter functions. Strips
    string and removes line breaks.

    """

    return value.strip().replace( "\n", "" )

def _getlist( value, dtype=None ):
    """
    Template for custom list getter functions.

    Takes 2 argument:

      value - String that contains the comma separated
              values to convert to a list. Expects the form:
                  value = "value_1,
                           ...,
                           value_n"
              Newliens are handled properly but not necessary.
      dtype - Optional datatype for the elements in the list.

    Returns 1 Value:

      list - corresponding list of the given datatype

    """

    value = _filter( value )
    if dtype is None:
        dtype = lambda x: x
    return [dtype( _filter( blob ) ) for blob in value.split(",")]

def _getdict( value, dtype=None ):
    """
    Template for custom dictionary getter functions.

    Takes 2 Arguments:

      value - String that contains the comma/colen separated
              values to convert to a dictionary. Expects the form:
                  value = key_1: value_1,
                          ...  : ...,
                          key_n: value_n
              Newlines are handled properly but not necessary.
      dtype - Optional datatype for the elements in the list.

    Returns 1 Value:

      dict - corresponding dictionary object of the given datatype.

    """

    value = _filter( value )
    if dtype is None:
        dtype = lambda x: x
    return {_filter( blob.split(":")[0] ) : dtype( _filter( blob.split(":")[1] ) )
             for blob in value.split(",")}

def _create_config():
    """
    Function to generate the custom config object. Creates a ConfigParser object with
    ExtendedInterpolation, AllowNoValues, and a number of additional getter functions.

    Takes No Arguments.

    Returns 1 Value:

      config - Empty ConfigPaser object
    """
    def _getarray_generator( dtype ):
        """
        Generates a getter function for a numpy array with dtype datatype.

        Takes 1 Argument:

          dtype - datatype cast function to apply to each string entry

        Returns 1 Value:

          getarray_lambda - function that converts a string into a numpy
                            array of the specified type
        """

        return lambda value: np.array( _getlist( value, dtype ) )
    def _getdict_generator( dtype ):
        """
        Generates a getter function for a dictionary with dtype datatype
        for each value.

        Takes 1 argument:

          dtype - datatype cast function to apply to each string entry

        Returns 1 Value:

          getdict_lambda - function that converts a string into a dictionary
                           with keys of the given type

        """

        return lambda value: _getdict( value, dtype )

    # This creates the following functions:
    #   config.getfloat32( value )       --> float32
    #   config.getfloat64( value )       --> float64
    #   config.getlist( value )          --> list
    #   config.getfloat32_array( value ) --> np.array( dtype=float32 )
    #   config.getfloat64_array( value ) --> np.array( dtype=float64 )
    #   config.getdict( value )          --> { str: str }
    #   config.getfloat32_dict( value )  --> { str: float32 }
    #   config.getfloat64_dict( value )  --> { str: float64 }

    converters = {
        "float32": np.float32,
        "float64": np.float64,
        "list":   _getlist,
        "float32_array": _getarray_generator( np.float32 ),
        "float64_array": _getarray_generator( np.float64 ),
        "dict": _getdict,
        "float32_dict": _getdict_generator( np.float32 ),
        "float64_dict": _getdict_generator( np.float64 )
    }
    return ConfigParser( interpolation=ExtendedInterpolation(),
                         allow_no_value=True,
                         converters=converters )

def display_config( section_keys=None ):
    """
    Displays specified sections of the current config. Defaults
    to displaying the entire config.

    Takes 1 argument:

      section_keys - An optional list of section keys to display.
                     Defaults to displaying all sections.

    Returns nothing.

    """

    current_config = get_config( validate=False )

    config_keys = current_config.sections()
    if section_keys is None:
        section_keys = config_keys

    # Print each section's entries if the
    # section exists. Otherwise, report
    # the section missing.
    for key in section_keys:
        if key not in config_keys:
            print( "-------- MISSING SECTION: {:s} -------- \n\n".format( key ) )
            continue

        section = dict( current_config[key] )

        # If section is empty, print the section title and continue
        if len( section.keys() ) == 0:
            print( "[{:s}]\n".format( key ))
            continue

        padding = max( [len( section_name ) + 1 for section_name in section.keys() ] )

        # We use f-strings instead of .format because
        # they support more advanced text alignment.
        section_string = "[{:s}]\n".format( key )
        for field, value in section.items():
            section_string += f"{field:<{padding}}= "
            if value == None:
                section_string = section_string[:-2] + "\n"
                continue

            # align multi-line values
            for line in value.split( "\n" ):
                section_string += f"{line}\n{' ':<{padding}}  "

            # removing trailing line
            section_string = "\n".join( section_string.split( "\n" )[:-1] ) + "\n"

        print( section_string )

def get_config( validate=True ):
    """
    Returns the current configuration if it has been loaded.

    Takes 1 Argument:

      validate - Optional boolean, determines whether to run
                 validate_config before returning config. Defaults
                 to True.

    Returns 1 value:

      config - The ConfigParser object for the loaded configuration.

    """
    global config

    if config is None:
        raise ValueError( "No config provided. Please load a config with "
                          "droplet_approximation.load_config( config_path ) "
                          "or set the config directly with "
                          "droplet_approximation.set_config( new_config )." )

    if validate:
        validate_config()

    return config

def get_config_as_dict():
    """
    Returns the config data as a dictionary.

    Takes no arguments.

    Returns 1 argument:

      config_dict - A dictionary of dictionaries
                    by section in the config.

    """
    config = get_config()
    config_dict = {}

    for section in config.sections():
        config_dict[section] = {}
        for option, value in config.items( section ):
            config_dict[section][option] = value

    return config_dict

def load_config( config_path ):
    """
    Load droplet_approximation config file from supplied path.

    Will also load any subconfig listed in the "autoloads" section.
    For example, if one had a config with NTLP's physical constants,
    one could add:
        physics: NTLP
    in the [autoloads] section of the config, and that profile would be
    loaded. If a profile is ommited, the corresponding defaults will be
    loaded instead.

    Takes 1 argument:

      config_path - path to .ini config file.

    Returns nothing.

    """
    global config, CONFIG_PATH, CONFIG_ROOT
    CONFIG_PATH = config_path
    CONFIG_ROOT = "/".join( config_path.split("/")[:-1] )

    config = _create_config()
    config.read( config_path )

    if not config.has_section( "autoloads" ):
        validate_config()
        return

    # Load in default subconfigs
    for subconfig, profile in dict( config["autoloads"] ).items():
        if profile is None:
            profile = config["default_profiles"].get( subconfig, None )
        if profile is None:
            raise ValueError( "No default specified for section \"{:s}\" autoload!".format( subconfig ) )

        load_subconfig( subconfig, profile )


def load_subconfig( subconfig, profile_name ):
    """
    Loads the .ini file corresponding to the profile/subconfig into the
    current config.

    Also checks for a an entry in the [default_profiles] section like so:
        subconfig: default_profile_name
    and loads this alongside the listed profile.

    For example, if you run load_subconfig( "simulation", "cfog" ) and the current
    config has the following entry in [default_profiles]:
        simulation: simulation_defaults
    then the config file simulation/simulation_defaults.ini will be loaded first. Then
    the config file simulation/cfog.ini will be loaded, potentially overwriting certain
    default values.

    Takes 1 argument:

      subconfig - Name of the subconfig to load.
      profile   - Name of the profile to load.

    Returns nothing.

    """
    global CONFIG_ROOT, config

    # Ensure that a config is loaded
    get_config( validate=False )

    # Store the old config to restore in case of an error
    old_config = deepcopy( config )

    # Ensure subconfig folder exists
    subconfig_folder   = Path( "{:s}/{:s}".format( CONFIG_ROOT, subconfig ) )
    if not subconfig_folder.is_dir():
        raise FileNotFoundError( "No folder for subconfig: {:s} at {}".format( subconfig, subconfig_folder ) )
    # Load defaults if they exist
    subconfig_default_profile = config["default_profiles"].get( subconfig, "" )
    if subconfig_default_profile != "":
        defaults_file = subconfig_folder / "{:s}.ini".format( config["default_profiles"][subconfig] )
        if not defaults_file.exists():
            raise FileNotFoundError( "No such default profile {:s}.ini ".format( subconfig_default_profile )
                                   + "for subconfig: {:s}".format( subconfig_default_profile, subconfig ) )

        print( "Loading defaults at: {}".format( defaults_file ) )
        config.read( defaults_file )

    # Load profile
    profile_filename = "{:s}.ini".format( profile_name )
    profile_file     = subconfig_folder / profile_filename
    if not profile_file.exists():
        raise FileNotFoundError( "No file for profile: {:s}".format( profile_filename ) )

    print( "Loading profile at: {}".format( profile_file ) )
    config.read( profile_file )

    try:
        validate_config()
    except ValueError as e:
        config = old_config
        raise( e )

def set_config_from_dict( config_dict, validate=True ):
    """
    Set the loaded config directly.

    Takes 2 arguments:

      new_config - the config to set.
      validate   - Optional boolean, determines whether to validate
                   the new config with validate_config(). Defaults to True.

    Returns nothing.

    """
    global config

    old_config = deepcopy( config )

    config     = _create_config()

    for section in config_dict:
        config.add_section( section )
        for option, value in config_dict[section].items():
            config.set( section, option, value )

    if validate:
        try:
            validate_config()
        except ValueError as e:
            config = old_config
            raise( e )

def validate_config( validation_sections=None ):
    """
    Check list of required options to make sure they all
    exist and are not None. Required options are listed
    in the [required] section of the config for the specified
    sections.

    If provided a list of sections, additionally verifies that those
    particular sections exist.

    Takes 1 Argument:

      validation_sections - Optional list of sections to check for existence.

    Returns nothing.

    """

    if validation_sections is None:
        validation_sections = []

    config = get_config( validate=False )

    # Create a list of sections/options that we are validating
    required_options = [requirement.split( "." ) for requirement in config["requires"]]
    validation_sections.extend( [section[0] for section in required_options] )

    # Create a list of missing sections and raise an error if its not empty
    missing_sections = [missing_section for missing_section in validation_sections
                        if not config.has_section( missing_section )]
    if missing_sections != []:
        raise ValueError( "Missing the following config sections: {:s}"
                              .format( ", ".join( np.unique( missing_sections ) ) )
                        + "\n\nCurrent sections:\n{:s}"
                              .format( ", ".join( config.sections() ) ) )

    # Create a lsit of missing options and raise an error if its not empty
    missing_options  = [".".join( missing_option ) for missing_option in required_options
                        if not config.has_option( missing_option[0], missing_option[-1] )
                           and len( missing_option ) == 2]
    if missing_options != []:
        raise ValueError( "Missing the following config options:\n{:s}"
                              .format( ", ".join( missing_options ) ) )

