import errno
import functools
import os
import signal

import numpy as np
from scipy.integrate import solve_ivp

# Provide the range of each of the input parameters so we can normalize
# each to [-1, 1], allowing the models to weight each equally.  If we don't
# perform this mapping then optimization process will effectively ignore
# the small parameters (e.g. radius) as they won't contribute as much as
# the big parameters (e.g. temperature).  As a result, we use log-scale
# for the parameters that have a large dynamic range.
DROPLET_RADIUS_LOG_RANGE        = np.array( (-8, -3) )
DROPLET_TEMPERATURE_RANGE       = np.array( (273, 310) )
DROPLET_SALT_MASS_LOG_RANGE     = np.array( (-22, -10) )
DROPLET_AIR_TEMPERATURE_RANGE   = np.array( (273, 310) )
DROPLET_RELATIVE_HUMIDITY_RANGE = np.array( (0.65, 1.1) )
DROPLET_RHOA_RANGE              = np.array( (0.8, 1.2) )
DROPLET_TIME_LOG_RANGE          = np.array( (-2.0, 1.0) )

# Tolerances for BDF solves.  Per solve_ivp() the effective tolerance is:
#
#   absolute_tolerance + (solution * relative_tolerance)
#
# Where the local error estimate used for convergence is kept within this value.
# We specify two absolute tolerances, one each for radius and temperature, since
# they have significantly different magnitudes and ranges.
#
# NOTE: These values are significantly tighter than solve_ivp()'s defaults
#       and are required for clean data with small particles (micron sized and
#       smaller, i.e. radius <= 1e-6m).
#
BDF_TOLERANCE_ABSOLUTE = (1e-10, 1e-4)
BDF_TOLERANCE_RELATIVE = 1e-7

def dydt( t, y, parameters ):
    """
    Differential equations governing a water droplet's radius and temperature as
    a function of time, given the specified physical parameters.

    Adapted from code provided by David Richter (droplet_integrator.py) in August 2024.
    Changes made include:

      - Accepting NumPy arrays for the parameters instead of having hard-coded values.

        NOTE: While it allows for multiple parameters to be specified via additional rows
              this does not necessarily work with scipy.solve_ivp()!

      - Parameters are promoted to 64-bit floating point values so all internal
        calculations are performed at maximum precision, resulting in 64-bit outputs.
        It is the caller's responsibility for casting the results to a different precision.

    Takes 3 arguments:

      t          - Unused argument.  Required for the use of scipy.solve_ivp().
      y          - NumPy array of droplet radii and temperatures.  May be specified
                   as either a 1D vector (of length 2), or a 2D array (sized
                   number_droplets x 2) though must be shape compatible with the
                   parameters array.
      parameters - NumPy array of droplet parameters containing salt mass, air temperature,
                   relative humidity, and rhoa.  May be specified as either a
                   1D vector (of length 2), or a 2D array (sized number_droplets x 2)
                   though must be shape compatible with the y array.

    Returns 2 values:

      dradius_dt      - The derivative of the droplet's radius with respect to time, shaped
                        number_droplets x 1.
      dtemperature_dt - The derivative of the droplet's temperature with respect to time,
                        shaped number_droplets x 1.

    """

    #
    # NOTE: We work in 64-bit precision regardless of the input
    #       so we get an as accurate as possible answer.
    #

    m_s  = parameters[0, ...].astype( "float64" )
    Tf   = parameters[1, ...].astype( "float64" )
    RH   = parameters[2, ...].astype( "float64" )
    rhoa = parameters[3, ...].astype( "float64" )

    rhow = np.float64( 1000 )
    rhos = np.float64( 2000 )
    #Cpp  = np.float64( 4190 )  #CM1
    Cpp  = np.float64( 4179 )  #NTLP
    Mw   = np.float64( 0.018015 )
    Ru   = np.float64( 8.3144 )
    Ms   = np.float64( 0.05844 )
    Gam  = np.float64( 7.28e-2 )
    #Ion  = np.float64( 2.0 )
    #Os   = np.float64( 1.093 )
    kap  = np.float64( 1.2 )
    Shp  = np.float64( 2 )
    Sc   = np.float64( 0.615 )
    Pra  = np.float64( 0.715 )
    Cpa  = np.float64( 1006.0 )
    nuf  = np.float64( 1.57e-5 )
    Lv   = np.float64( (25.0 - 0.02274*26)*10**5 )
    Nup  = np.float64( 2 )

    #einf = 611.2*np.exp(17.67*(Tf-273.15)/(Tf-29.65))  #CM1, Bolton (1980, MWR)
    einf = 610.94*np.exp((17.6257*(Tf-273.15))/(243.04+(Tf-273.15)))  #NTLP

    qinf = RH/rhoa*(einf*Mw/Ru/Tf)

    Volp = 4/3*np.pi*y[0]**3
    rhop = (m_s + Volp*rhow)/Volp
    taup = rhop*(2*y[0])**2/18/nuf/rhoa


    #qstar = einf*Mw/Ru/y[1]/rhoa*np.exp(Lv*Mw/Ru*(1/Tf - 1/y[1]) + 2*Mw*Gam/Ru/rhow/y[0]/y[1] - Ion*Os*m_s*(Mw/Ms)/(Volp*rhow))
    #exp_stuff = Lv*Mw/Ru*(1/Tf - 1/y[1]) + 2*Mw*Gam/Ru/rhow/y[0]/y[1] - (kap*m_s*rhow/rhos)/(Volp*rhow)
    #print(f"exp_stuff = {exp_stuff:.17g}")
    #qstar = einf*Mw/Ru/y[1]/rhoa*np.exp(exp_stuff)
    term1 = Lv*Mw/Ru*(1/Tf - 1/y[1])
    term2 = 2*Mw*Gam/Ru/rhow/y[0]/y[1]
    term3 = (kap*m_s*rhow/rhos)/(Volp*rhow) # !!! Shouldn't this be Mw/Ms not density??
    exp_stuff = term1 + term2  - term3
    #if exp_stuff > 10.0:
    #    print(f"exp_stuff = {exp_stuff:.17g}")
    #    print(f"term1 = {term1:.17g}")
    #    print(f"term2 = {term2:.17g}")
    #    print(f"term3 = {term3:.17g}")
    #    print(f"Volp = {Volp:.17g}")
    qstar = einf*Mw/Ru/y[1]/rhoa*np.exp(exp_stuff)


    dy1dt = 1/9*Shp/Sc*rhop/rhow*y[0]/taup*(qinf - qstar)
    dy2dt = -1/3*Nup/Pra*Cpa/Cpp*rhop/rhow/taup*(y[1] - Tf) + 3*Lv/y[0]/Cpp*dy1dt

    return [dy1dt, dy2dt]

def get_parameter_ranges():
    """
    Gets the current droplet parameter ranges as a dictionary mapping
    droplet property names to a two element 1D NumPy array containing
    the property's lower and upper bounds, respectively.

    Takes no arguments.

    Returns 1 value:

      parameter_ranges - Dictionary with the following keys: "radius",
                         "temperature", "salt_mass", "air_temperature",
                         "relative_humidity", "rhoa", "time"

    """

    parameter_ranges = {}

    parameter_ranges["radius"]            = DROPLET_RADIUS_LOG_RANGE
    parameter_ranges["temperature"]       = DROPLET_TEMPERATURE_RANGE
    parameter_ranges["salt_mass"]         = DROPLET_SALT_MASS_LOG_RANGE
    parameter_ranges["air_temperature"]   = DROPLET_AIR_TEMPERATURE_RANGE
    parameter_ranges["relative_humidity"] = DROPLET_RELATIVE_HUMIDITY_RANGE
    parameter_ranges["rhoa"]              = DROPLET_RHOA_RANGE
    parameter_ranges["time"]              = DROPLET_TIME_LOG_RANGE

    return parameter_ranges

def set_parameter_ranges( parameter_ranges ):
    """
    Sets the droplet parameter ranges based on the provided upper
    and lower bounds.  One or more parameter ranges may be specified
    depending on the keys present in the specified dictionary, only
    those present are updated leaving the unspecified properties as
    is.

    Takes 1 argument:

      parameter_ranges - Dictionary with one or more of the following keys:
                         "radius", "temperature", "salt_mass",
                         "air_temperature", "relative_humidity", "rhoa", "time"

    Returns nothing.

    """

    global DROPLET_RADIUS_LOG_RANGE, DROPLET_TEMPERATURE_RANGE, \
           DROPLET_SALT_MASS_LOG_RANGE, DROPLET_AIR_TEMPERATURE_RANGE, \
           DROPLET_RELATIVE_HUMIDITY_RANGE, DROPLET_RHOA_RANGE, \
           DROPLET_TIME_LOG_RANGE

    if "radius" in parameter_ranges:
        DROPLET_RADIUS_LOG_RANGE = parameter_ranges["radius"]
    if "temperature" in parameter_ranges:
        DROPLET_TEMPERATURE_RANGE = parameter_ranges["temperature"]
    if "salt_mass" in parameter_ranges:
        DROPLET_SALT_MASS_LOG_RANGE = parameter_ranges["salt_mass"]
    if "air_temperature" in parameter_ranges:
        DROPLET_AIR_TEMPERATURE_RANGE = parameter_ranges["air_temperature"]
    if "relative_humidity" in parameter_ranges:
        DROPLET_RELATIVE_HUMIDITY_RANGE = parameter_ranges["relative_humidity"]
    if "rhoa" in parameter_ranges:
        DROPLET_RHOA_RANGE = parameter_ranges["rhoa"]
    if "time" in parameter_ranges:
        DROPLET_TIME_LOG_RANGE = parameter_ranges["time"]

def normalize_droplet_parameters( droplet_parameters ):
    """
    Normalizes an array of droplet parameters into the range [-1, 1].  Operates
    on both input and parameters.

    Takes 1 argument:

      droplet_parameters - NumPy array of parameters to normalize, shaped
                           number_droplets x number_parameters.  May either be
                           1D, for a single droplet, or 2D for multiple droplets.
                           number_parameters must be either 2 (output parameters),
                           6 (input parameters without t_final), or 7 (input
                           parameters with t_final).

    Returns 1 value:

      normalized_droplet_parameters - NumPy array of normalized parameters, shaped
                                      number_droplets x number_parameters.  Has the
                                      same shape as droplet_parameters.

    """

    #
    # NOTE: We take care to handle arbitrary rank arrays!  While we expect either
    #       rank-1 or rank-2, the code is written to handle larger ranks naturally
    #       with the inner dimension representing a single droplet's parameters.
    #

    number_parameters = droplet_parameters.shape[-1]

    # Bail if we didn't get 1) output parameters, 2) input parameters without a t_final,
    # or 3) input parameters with a t_final.
    if number_parameters != 2 and number_parameters != 6 and number_parameters != 7:
        raise ValueError( "Unknown number of parameters to normalize ({:d})!".format( number_parameters ) )

    normalized_droplet_parameters = np.empty_like( droplet_parameters )

    # We always have radius and temperature.
    normalized_droplet_parameters[..., 0] = (np.log10( droplet_parameters[..., 0] ) - np.mean( DROPLET_RADIUS_LOG_RANGE )) / (np.diff( DROPLET_RADIUS_LOG_RANGE ) / 2)
    normalized_droplet_parameters[..., 1] = (droplet_parameters[..., 1] - np.mean( DROPLET_TEMPERATURE_RANGE )) / (np.diff( DROPLET_TEMPERATURE_RANGE ) / 2)

    # Sometimes we have the remaining parameters.
    if number_parameters > 2:
        normalized_droplet_parameters[..., 2] = (np.log10( droplet_parameters[..., 2] ) - np.mean( DROPLET_SALT_MASS_LOG_RANGE )) / (np.diff( DROPLET_SALT_MASS_LOG_RANGE ) / 2)
        normalized_droplet_parameters[..., 3] = (droplet_parameters[..., 3] - np.mean( DROPLET_AIR_TEMPERATURE_RANGE )) / (np.diff( DROPLET_AIR_TEMPERATURE_RANGE ) / 2)
        normalized_droplet_parameters[..., 4] = (droplet_parameters[..., 4] - np.mean( DROPLET_RELATIVE_HUMIDITY_RANGE )) / (np.diff( DROPLET_RELATIVE_HUMIDITY_RANGE ) / 2)
        normalized_droplet_parameters[..., 5] = (droplet_parameters[..., 5] - np.mean( DROPLET_RHOA_RANGE )) / (np.diff( DROPLET_RHOA_RANGE ) / 2)

        # Copy the evaluation times when they're present.
        if number_parameters > 6:
            normalized_droplet_parameters[..., 6] = droplet_parameters[..., 6]

    return normalized_droplet_parameters

def scale_droplet_parameters( droplet_parameters ):
    """
    Scales an array of droplet parameters from the range [-1, 1] to their
    expected ranges of physical values.  Operates on both input and parameters.

    Takes 1 argument:

      droplet_parameters - NumPy array of parameters to scale, shaped
                           number_droplets x number_parameters.  May either be
                           1D, for a single droplet, or 2D for multiple droplets.
                           number_parameters must be either 2 (output parameters),
                           6 (input parameters without t_final), or 7 (input
                           parameters with t_final).

    Returns 1 value:

      scaled_droplet_parameters - NumPy array of scaled parameters, shaped
                                  number_droplets x number_parameters.  Has the
                                  same shape as droplet_parameters.


    """

    #
    # NOTE: We take care to handle arbitrary rank arrays!  While we expect either
    #       rank-1 or rank-2, the code is written to handle larger ranks naturally
    #       with the inner dimension representing a single droplet's parameters.
    #

    number_parameters = droplet_parameters.shape[-1]

    # Bail if we didn't get 1) output parameters, 2) input parameters without a t_final,
    # or 3) input parameters with a t_final.
    if number_parameters != 2 and number_parameters != 6 and number_parameters != 7:
        raise ValueError( "Unknown number of parameters to scale ({:d})!".format( number_parameters ) )

    scaled_droplet_parameters = np.empty_like( droplet_parameters )

    # We always have radius and temperature.
    scaled_droplet_parameters[..., 0] = 10.0 ** (droplet_parameters[..., 0] * (np.diff( DROPLET_RADIUS_LOG_RANGE ) / 2) + np.mean( DROPLET_RADIUS_LOG_RANGE ))
    scaled_droplet_parameters[..., 1] = droplet_parameters[..., 1] * (np.diff( DROPLET_TEMPERATURE_RANGE ) / 2) + np.mean( DROPLET_TEMPERATURE_RANGE )

    # Sometimes we have the remaining parameters.
    if number_parameters > 2:
        scaled_droplet_parameters[..., 2] = 10.0 ** (droplet_parameters[..., 2] * (np.diff( DROPLET_SALT_MASS_LOG_RANGE ) / 2) + np.mean( DROPLET_SALT_MASS_LOG_RANGE ))
        scaled_droplet_parameters[..., 3] = droplet_parameters[..., 3] * (np.diff( DROPLET_AIR_TEMPERATURE_RANGE ) / 2) + np.mean( DROPLET_AIR_TEMPERATURE_RANGE )
        scaled_droplet_parameters[..., 4] = droplet_parameters[..., 4] * (np.diff( DROPLET_RELATIVE_HUMIDITY_RANGE ) / 2) + np.mean( DROPLET_RELATIVE_HUMIDITY_RANGE )
        scaled_droplet_parameters[..., 5] = droplet_parameters[..., 5] * (np.diff( DROPLET_RHOA_RANGE ) / 2) + np.mean( DROPLET_RHOA_RANGE )

        # Copy the evaluation times when they're present.
        if number_parameters > 6:
            scaled_droplet_parameters[..., 6] = droplet_parameters[..., 6]

    return scaled_droplet_parameters

def solve_ivp_float32_outputs( dydt, t_span, y0, **kwargs ):
    """
    Solves an initial value problem and returns the solutions in 32-bit precision.

    NOTE: This is a wrapper around scipy.integrate.solve_ivp(), so it takes the
          same arguments and returns the same values.  See that function's
          help for a (way more) detailed explanation of each argument and value.

    Takes 4 arguments:

      dydt   - Right-hand side of the system to solve.  The calling signature
               is 'dydt( t, y, parameters )'.
      t_span - 2-member sequence specifying the interval of integration.  dydt
               is integrated from t_span[0] until t_span[1].
      y0     - Initial state of the system to solve.  Must be compatible with
               dydt.
      kwargs - Optional keyword arguments to supply to solve_ivp().

    Returns 1 value:

      solution - Object containing fields related to the integration process,
                 including:

                   t - Vector of times, of length number_times, where dydt was
                       evaluated.
                   y - Array of solutions, shaped number_variables x numbe_times,
                       where dydt was evaluated.

                 Solutions will be returned according to the evaluation window,
                 t_eval, either supplied via kwargs or selected as a default from
                 t_span.

    """

    # Solve the ODE in the precision supplied by the caller.
    solution = solve_ivp( dydt, t_span, y0, **kwargs )

    # Return the outputs as the requested precision when a solution was found.
    if solution.success:
        solution.t = solution.t.astype( "float32" )
        solution.y = solution.y.astype( "float32" )

    return solution

class TimeoutError( Exception ):
    pass

def timeout( seconds=10, error_message=os.strerror( errno.ETIME ) ):

    def decorator( func ):
        def _handle_timeout( signum, frame ):
            raise TimeoutError( error_message )

        @functools.wraps( func )
        def wrapper( *args, **kwargs ):
            signal.signal( signal.SIGALRM, _handle_timeout )
            signal.alarm( seconds )

            try:
                result = func( *args, **kwargs )
            finally:
                signal.alarm( 0 )

            return result

        return wrapper

    return decorator

@timeout( 5 )
def timed_solve_ivp( *args, **kwargs ):
    """
    Solves an initial value problem and returns the solutions in 32-bit precision,
    though aborts execution if the solution is not available within 5 seconds.
    This works around the fact that dydt() sometimes converges *VERY* slowly as
    some input parameter combinations result in incredibly small steps that
    results in the solution appearing to hang.

    NOTE: This is a wrapper around scipy.integrate.solve_ivp(), so it takes the
          same arguments and returns the same values.  See that function's
          help for a (way more) detailed explanation of each argument and value.

    Takes 4 arguments:

      dydt   - Right-hand side of the system to solve.  The calling signature
               is 'dydt( t, y, parameters )'.
      t_span - 2-member sequence specifying the interval of integration.  dydt
               is integrated from t_span[0] until t_span[1].
      y0     - Initial state of the system to solve.  Must be compatible with
               dydt.
      kwargs - Optional keyword arguments to supply to solve_ivp().

    Returns 1 value:

      solution - Object containing fields related to the integration process,
                 including:

                   t - Vector of times, of length number_times, where dydt was
                       evaluated.
                   y - Array of solutions, shaped number_variables x numbe_times,
                       where dydt was evaluated.

                 Solutions will be returned according to the evaluation window,
                 t_eval, either supplied via kwargs or selected as a default from
                 t_span.
    """

    return solve_ivp_float32_outputs( *args, **kwargs )
