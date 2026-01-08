from enum import Enum
import errno
import functools
import os
import signal
import warnings

import numpy as np
import torch
from scipy.integrate import solve_ivp

# Temperature delta, in Kelvin, relative to the surrounding air to sample
# droplet temperatures from.  This constrains generated droplets to physically
# realizable scenario by respecting the implicit coupling of droplet and the
# space it occupies instead of leaving them fully independent (and creating
# exciting physics as a result).
#
# NOTE: 3 K was somewhat arbitrarily chosen based on reviewing NTLP simulation
#       outputs, but not thoroughly investigated.
#
DROPLET_DELTA_AIR_TEMPERATURE = 3.0

# Provide the range of each of the input parameters so we can normalize
# each to [-1, 1], allowing the models to weight each equally.  If we don't
# perform this mapping then optimization process will effectively ignore
# the small parameters (e.g. radius) as they won't contribute as much as
# the big parameters (e.g. temperature).  As a result, we use log-scale
# for the parameters that have a large dynamic range.
#
# NOTE: We force the data type of each array to avoid NumPy's "help" should
#       someone specify both ranges as integers.
#
DROPLET_RADIUS_LOG_RANGE        = np.array( (-8.0, -3.0),   dtype=np.float64 )
DROPLET_TEMPERATURE_RANGE       = np.array( (273.0, 310.0), dtype=np.float64 )
DROPLET_SALT_SOLUTE_LOG_RANGE   = np.array( (-22.0, -8.0),  dtype=np.float64 )
DROPLET_AIR_TEMPERATURE_RANGE   = np.array( (273.0, 310.0), dtype=np.float64 )
DROPLET_RELATIVE_HUMIDITY_RANGE = np.array( (0.55, 1.1),    dtype=np.float64 )
DROPLET_RHOA_RANGE              = np.array( (0.8, 1.3),     dtype=np.float64 )
DROPLET_TIME_LOG_RANGE          = np.array( (-2.0, 1.0),    dtype=np.float64 )

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

# Density of fresh water (kg/m^3).
RHOW = 1000.0

class DisplayType( Enum ):
    """
    Enumeration class specifying how particle parameter ranges are displayed.
    """

    # Show the ranges in human-readable form with units.
    HUMAN   = 0

    # Show the ranges as code that can be copy and pasted.
    CODE    = 1

    # Show the ranges in a form that can be parsed on the command line.
    COMMAND = 2

def display_parameter_ranges( parameter_ranges, display_type=DisplayType.HUMAN, indentation_str="" ):
    """
    Displays the supplied parameter ranges in a variety of formats depending on
    how they'll be used.  Human-readable (with units and details), code (that
    can be copied verbatim), and command line (that can be parsed) formats are
    supported.

    Takes 3 arguments:

      parameter_ranges - Dictionary of parameter ranges from get_parameter_ranges().
                         This must have all of the parameter keys!
      display_type     - Optional DisplayType enumeration specifying how to
                         display the ranges.

                            .CODE:    Ranges are printed as code that can be copied
                                      and pasted as is.
                            .COMMAND: Ranges are printed compactly in a format that
                                      can be parsed on the command line.
                            .HUMAN:   Ranges are printed for human consumption with
                                      units and additional details for interpretation.

                         If omitted, defaults to DisplayType.HUMAN.

      indentation_str  - Optional indentation string that is applied before each
                         line of output.  If omitted, defaults to an empty
                         string resulting in no indentation.

    Returns nothing.

    """

    TEMPLATE_HUMAN = \
        "{indent:s}Data ranges:\n" + \
        "\n" + \
        "{indent:s}  log10(Radius)       [{{:.2f}}, {{:.2f}}] m\n" + \
        "{indent:s}  Temperature:        [{{:.2f}}, {{:.2f}}] K\n" + \
        "{indent:s}  log10(Salt solute): [{{:.2f}}, {{:.2f}}] kg\n" + \
        "{indent:s}  Air temperature:    [{{:.2f}}, {{:.2f}}] K\n" + \
        "{indent:s}  Relative humidity:  [{{:.2f}}, {{:.2f}}] %\n" + \
        "{indent:s}  Air density:        [{{:.2f}}, {{:.2f}}] kg/m^3\n" + \
        "{indent:s}  log10(dt):          [{{:.2f}}, {{:.2f}}] s"
    TEMPLATE_CODE = \
        "{indent:s}parameter_ranges = {{{{\n" + \
        "{indent:s}    \"radius\":             ({{:f}}, {{:f}}),\n" + \
        "{indent:s}    \"temperature\":        ({{:f}}, {{:f}}),\n" + \
        "{indent:s}    \"salt_solute\":        ({{:f}}, {{:f}}),\n" + \
        "{indent:s}    \"air_temperature\":    ({{:f}}, {{:f}}),\n" + \
        "{indent:s}    \"relative_humidity\":  ({{:f}}, {{:f}}),\n" + \
        "{indent:s}    \"rhoa\":               ({{:f}}, {{:f}}),\n" + \
        "{indent:s}    \"time\":               ({{:f}}, {{:f}})\n" + \
        "{indent:s}}}}}"
    TEMPLATE_COMMAND = \
        "{indent:s}{{:f}}:{{:f}},{{:f}}:{{:f}},{{:f}}:{{:f}},{{:f}}:{{:f}},{{:f}}:{{:f}},{{:f}}:{{:f}},{{:f}}:{{:f}}"

    # Instantiate an indented template.  Blow up if we don't know which
    # template to use.
    if display_type == DisplayType.CODE:
        template = TEMPLATE_CODE.format( indent=indentation_str )
    elif display_type == DisplayType.COMMAND:
        template = TEMPLATE_COMMAND.format( indent=indentation_str )
    elif display_type == DisplayType.HUMAN:
        template = TEMPLATE_HUMAN.format( indent=indentation_str )
    else:
        raise ValueError( "Unknown DisplayType provided ({})!  Must be one of: '{:s}', '{:s}', '{:s}'.".format(
            display_type,
            "human",
            "code",
            "command_line" ) )

    # Render the ranges according to the template.
    display_str = template.format(
        parameter_ranges["radius"][0],            parameter_ranges["radius"][1],
        parameter_ranges["temperature"][0],       parameter_ranges["temperature"][1],
        parameter_ranges["salt_solute"][0],       parameter_ranges["salt_solute"][1],
        parameter_ranges["air_temperature"][0],   parameter_ranges["air_temperature"][1],
        parameter_ranges["relative_humidity"][0], parameter_ranges["relative_humidity"][1],
        parameter_ranges["rhoa"][0],              parameter_ranges["rhoa"][1],
        parameter_ranges["time"][0],              parameter_ranges["time"][1] )

    # Display the range.
    print( display_str )

def droplet_equilibrium( droplet_parameters ):
    """
    Uses ported BE code to calculate the droplet equilibrium.

    Takes 1 Argument:
      droplet_parameters - NumPy array of parameters to calculate the equilibria of.
                           Sized number_droplets x 6.  May either be
                           1D, for a single droplet, or 2D for multiple droplets.
                           The parameters are, in order: radius, particle temperature,
                           solute term (NOT salt mass), air temperature, relative humidity,
                           and air density.

    Returns 2 Values:
      equilibria - Numpy array of radius/temperature equilibria. If no equilibrium is found, returns
                   the inputed radius/temperature. This follows the behavior found in NTLP.
      flags      - Whether a successful equilibrium was found. Note: if RH > 100%, the equilibrium
                   only has imaginary roots ==> there will be flags.
    """

    r           = droplet_parameters[..., 0]
    Tp          = droplet_parameters[..., 1]
    solute_term = droplet_parameters[..., 2]
    Tf          = droplet_parameters[..., 3]
    RH          = droplet_parameters[..., 4]
    rhoa        = droplet_parameters[..., 5]

    rhow = np.float64( RHOW )
    rhos = np.float64( 2000 )
    Mw   = np.float64( 0.018015 )
    Ru   = np.float64( 8.3144 )
    Gam  = np.float64( 7.28e-2 )

    flag  = np.zeros_like( r )

    einf = 610.94*np.exp((17.6257*(Tf-273.15))/(243.04+(Tf-273.15)))  #NTLP
    qinf = RH/rhoa*(einf*Mw/Ru/Tf)
    #qinf = RH/rhoa*(einf*Mw/Ru/Tf)

    pi2   = 2.0 * np.pi

    a = -(2*Mw*Gam)/(Ru*rhow*Tf)/np.log((Ru*Tf*rhoa*qinf)/(Mw*einf))
    # This line was changed from BE to account for using
    # solute term instead of salt mass. There was no rhow
    # term present in the original equation, so 1/rhow had
    # to be added to correct for the rhow coefficient in solute_term.
    c = (solute_term)/((2.0/3.0)*pi2*rhow)/np.log((Ru*Tf*rhoa*qinf)/(Mw*einf))

    Q = (a**2.0)/9.0
    R = (2.0*a**3.0+27.0*c)/54.0
    M = R**2.0-Q**3.0
    val = (R**2.0)/(Q**3.0)

    guess = Q

    mask = M<0

    theta = np.arccos(R/np.sqrt(Q**3.0))
    guess[mask] = (-(2*np.sqrt(Q)*np.cos((theta-pi2)/3.0))-a/3.0)[mask]

    guess[guess<0] = (-(2*np.sqrt(Q)*np.cos((theta+pi2)/3.0))-a/3.0)[guess < 0]

    S = -(R/np.abs(R))*(np.abs(R)+np.sqrt(M))**(1.0/3.0)
    T = Q/S
    guess[~mask] = (S + T - a/3.0)[~mask]

    flag[guess < 0] = 1
    guess[guess < 0] = r[guess < 0]

    return guess, flag

def dydt( t, y, parameters ):
    """
    Differential equations governing a water droplet's radius and temperature as
    a function of time, given the specified physical parameters.

    The droplets' inputs may be specified either as NumPy arrays or as PyTorch
    tensors.  The background parameters must be supplied as the same type.

    NOTE: Derivatives are computed with 32-bit floating point when PyTorch
          tensors are provided.  This is done to ensure they can be computed at
          reasonable speed as not all accelerators support 64-bit arithmetic and
          this function should not be used to create training data that *must*
          be calculated in the highest precision (read: we only support calling
          SciPy's solve_ivp() with NumPy arrays).

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
      y          - Array of droplet radii and temperatures.  May be specified as
                   either a 1D vector (of length 2), or a 2D array (sized
                   number_droplets x 2) though must be type and shape compatible
                   with the parameters array.  Must be either a NumPy array or a
                   PyTorch tensor.
      parameters - Array of background parameters containing solute term, air temperature,
                   relative humidity, and rhoa.  May be specified as either a
                   1D vector (of length 4), or a 2D array (sized number_droplets x 4)
                   though must be type and shape compatible with the y array.

    Returns 2 values:

      dradius_dt      - The derivative of the droplet's radius with respect to time, shaped
                        number_droplets x 1.  This has the same type and location as y.
      dtemperature_dt - The derivative of the droplet's temperature with respect to time,
                        shaped number_droplets x 1.  This has the same type and location as y.

    """

    # XXX: take an optional dictionary of constants so these aren't recreated
    #      every time
    CPP  = 4179.0  #CM1: 4190
    MW   = 0.018015
    RU   = 8.3144
    GAM  = 7.28e-2
    SHP  = 2.0
    SC   = 0.615
    PRA  = 0.715
    CPA  = 1006.0
    NUF  = 1.57e-5
    LV   = (25.0 - 0.02274*26)*10**5
    NUP  = 2

    # Get our constants
    if isinstance( y, np.ndarray ):
        xp = np

        #
        # NOTE: We work in 64-bit precision regardless of the input
        #       so we get an as accurate as possible answer.
        #
        r    = y[..., 0].astype( "float64" )
        Tp   = y[..., 1].astype( "float64" )

        solute_term = parameters[..., 0].astype( "float64" )
        Tf          = parameters[..., 1].astype( "float64" )
        RH          = parameters[..., 2].astype( "float64" )
        rhoa        = parameters[..., 3].astype( "float64" )

        rhow = np.float64( RHOW )
        Cpp  = np.float64( CPP )  #NTLP
        Mw   = np.float64( MW )
        Ru   = np.float64( RU )
        Gam  = np.float64( GAM )
        Shp  = np.float64( SHP )
        Sc   = np.float64( SC )
        Pra  = np.float64( PRA )
        Cpa  = np.float64( CPA )
        nuf  = np.float64( NUF )
        Lv   = np.float64( LV )
        Nup  = np.float64( NUP )
    else:
        xp = torch

        #
        # NOTE: We calculate in the highest precision that is portable across
        #       accelerators.  64-bit arithmetic is not always hardware
        #       accelerated (e.g. consumer-grade NVIDIA GPUs).  This is an
        #       intentional choice since we don't support PyTorch tensors to
        #       generate data which must be calculated in 64-bit arithmetic.
        #
        float_dtype = torch.float32
        device      = y.device

        r    = y[..., 0].to( device=device, dtype=float_dtype )
        Tp   = y[..., 1].to( device=device, dtype=float_dtype )

        solute_term = parameters[..., 0].to( device=device, dtype=float_dtype )
        Tf          = parameters[..., 1].to( device=device, dtype=float_dtype )
        RH          = parameters[..., 2].to( device=device, dtype=float_dtype )
        rhoa        = parameters[..., 3].to( device=device, dtype=float_dtype )

        rhow = torch.Tensor( [RHOW] ).to( device=device, dtype=float_dtype )
        Cpp  = torch.Tensor( [CPP] ).to( device=device, dtype=float_dtype )  #NTLP
        Mw   = torch.Tensor( [MW] ).to( device=device, dtype=float_dtype )
        Ru   = torch.Tensor( [RU] ).to( device=device, dtype=float_dtype )
        Gam  = torch.Tensor( [GAM] ).to( device=device, dtype=float_dtype )
        Shp  = torch.Tensor( [SHP] ).to( device=device, dtype=float_dtype )
        Sc   = torch.Tensor( [SC] ).to( device=device, dtype=float_dtype )
        Pra  = torch.Tensor( [PRA] ).to( device=device, dtype=float_dtype )
        Cpa  = torch.Tensor( [CPA] ).to( device=device, dtype=float_dtype )
        nuf  = torch.Tensor( [NUF] ).to( device=device, dtype=float_dtype )
        Lv   = torch.Tensor( [LV] ).to( device=device, dtype=float_dtype )
        Nup  = torch.Tensor( [NUP] ).to( device=device, dtype=float_dtype )

    Volp = 4/3*np.pi*r**3

    term1     = Lv*Mw/Ru*(1/Tf - 1/Tp)
    term2     = 2*Mw*Gam/Ru/rhow/r/Tp
    term3     = solute_term/(Volp*rhow)
    exp_stuff = term1 + term2  - term3

    #einf = 611.2*xp.exp(17.67*(Tf-273.15)/(Tf-29.65))  #CM1, Bolton (1980, MWR)
    einf  = 610.94*xp.exp((17.6257*(Tf-273.15))/(243.04+(Tf-273.15)))  #NTLP
    qinf  = RH/rhoa*(einf*Mw/Ru/Tf)
    qstar = einf*Mw/Ru/Tp/rhoa*xp.exp(exp_stuff)

    dy1dt = 1/2*Shp*nuf*rhoa/Sc/rhow/r*(qinf - qstar)
    dy2dt = -3/2*Nup*nuf*rhoa/Pra*Cpa/Cpp/rhow/r**2*(Tp - Tf) + 3*Lv/r/Cpp*dy1dt

    return [dy1dt, dy2dt]

def dydt_mass( t, y, parameters ):
    """
    Differential equations governing a water droplet's radius and temperature as
    a function of time, given the specified physical parameters.

    NOTE: This function does *NOT* handle PyTorch tensors like dydt() does!

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

    rhow = np.float64( RHOW )
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

def generate_random_droplets( number_droplets, max_salinity=np.inf, normalize_flag=False ):
    """
    Generates random droplet parameters.  Takes care to conditionally sample
    parameters so as to create physically realistic droplets.

    Takes 3 arguments:

      number_droplets - Number of droplets to generate.
      max_salinity    - Optional scalar floating point, in the range of [0,
                        np.inf) specifying the maximum salinity allowed for a
                        droplet.  Droplets whose salt solute mass are greater
                        than this salinity will have their mass adjusted to a
                        smaller value that satisfies this maximum.  Salinity
                        is defined as the ratio of salt solute divided by the
                        droplet's water volume times the dimensionless
                        density of fresh water.  If omitted, defaults to
                        np.inf which corresponds to wet salt and admits any
                        salt solute mass.
      normalize_flag  - Optional flag indicating the droplets generated should be
                        normalized into the range of [-1, 1].  If omitted,
                        defaults to False and the droplets are in physical
                        units.

    Returns 1 value:

      droplet_parameters - NumPy array, sized number_droplets x 6, containing
                           the generated non-temporal droplet parameters.

    """

    # Uniformly sample the normalized ranges to start.
    droplet_parameters = np.reshape( np.random.uniform( -1, 1, number_droplets*6 ),
                                     (number_droplets, 6) ).astype( "float32" )
    droplet_parameters = scale_droplet_parameters( droplet_parameters )

    # Conditionally sample the droplets' temperatures based on their surrounding
    # air temperatures.
    droplet_parameters[:, 1] = (droplet_parameters[:, 3] +
                                DROPLET_DELTA_AIR_TEMPERATURE * np.random.uniform( -1, 1, number_droplets ))

    # Ensure the salt solute terms respect the maximum salinity requested.  This
    # calculates the maximum admissible salt solute term for the input droplet
    # volumes and scales the randomly chosen solute term between the solute
    # range's minimum and the per-droplet maximum.  As a result, this samples
    # the log-uniform solute distribute limited by the caller's salinity
    # request.
    #
    # NOTE: Infinite maximum salinity (i.e. wet salt) is correctly handled.
    #
    # NOTE: Zero maximum salinity will still admit saline particles as this
    #       never produces salt solutes that are outside of their parameter
    #       range.
    #
    log_solute_range         = get_parameter_ranges()["salt_solute"]
    effective_max_log_solute = np.log10( max_salinity * RHOW *
                                         (4.0/3.0 * np.pi * (droplet_parameters[:, 0]**3)) )
    log_solute_coefficient   = 1.0 - np.maximum( (log_solute_range[1] - effective_max_log_solute) / np.diff( log_solute_range ),
                                                 0.0 )
    droplet_parameters[:, 2] = 10.0**(log_solute_coefficient * (np.log10( droplet_parameters[:, 2] ) -
                                                                log_solute_range[0] ) +
                                      log_solute_range[0])

    # Map everything into [-1, 1] if requested.
    if normalize_flag:
        droplet_parameters = normalize_droplet_parameters( droplet_parameters )

    return droplet_parameters

def get_parameter_ranges( tensor_flag=False, tensor_device=None ):
    """
    Gets the current droplet parameter ranges as a dictionary mapping
    droplet property names to a two element 1D NumPy array containing
    the property's lower and upper bounds, respectively.

    Takes 2 arguments:

      tensor_flag   - Optional Boolean flag specifying whether the parameter
                      ranges should be returned as NumPy arrays (False) or PyTorch
                      tensors (True).  If omitted, defaults to False.
      tensor_device - Optional device specifying where parameter range tensors
                      should reside and may be any value accepted by PyTorch
                      tensor's .to() method.  This argument is ignored when
                      tensor_flag is False.  If omitted, defaults to "cpu" which
                      keeps the parameter ranges in host memory.

    Returns 1 value:

      parameter_ranges - Dictionary with the following keys: "radius",
                         "temperature", "salt_solute", "air_temperature",
                         "relative_humidity", "rhoa", "time"

    """

    # Default tensors on the host to ensure they can always be used.
    if tensor_device is None:
        device_str = "cpu"
    else:
        device_str = tensor_device

    parameter_ranges = {}

    parameter_ranges["radius"]            = DROPLET_RADIUS_LOG_RANGE
    parameter_ranges["temperature"]       = DROPLET_TEMPERATURE_RANGE
    parameter_ranges["salt_solute"]       = DROPLET_SALT_SOLUTE_LOG_RANGE
    parameter_ranges["air_temperature"]   = DROPLET_AIR_TEMPERATURE_RANGE
    parameter_ranges["relative_humidity"] = DROPLET_RELATIVE_HUMIDITY_RANGE
    parameter_ranges["rhoa"]              = DROPLET_RHOA_RANGE
    parameter_ranges["time"]              = DROPLET_TIME_LOG_RANGE

    # Convert the ranges from NumPy arrays to Torch Tensors when requested.
    if tensor_flag:
        for parameter_name, parameter_range in parameter_ranges.items():
            parameter_ranges[parameter_name] = torch.from_numpy( parameter_range ).type( torch.float32 ).to( device_str )

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
                         "radius", "temperature", "salt_solute",
                         "air_temperature", "relative_humidity", "rhoa", "time"

    Returns nothing.

    """

    global DROPLET_RADIUS_LOG_RANGE, DROPLET_TEMPERATURE_RANGE, \
           DROPLET_SALT_SOLUTE_LOG_RANGE, DROPLET_AIR_TEMPERATURE_RANGE, \
           DROPLET_RELATIVE_HUMIDITY_RANGE, DROPLET_RHOA_RANGE, \
           DROPLET_TIME_LOG_RANGE

    if "radius" in parameter_ranges:
        DROPLET_RADIUS_LOG_RANGE = parameter_ranges["radius"]
    if "temperature" in parameter_ranges:
        DROPLET_TEMPERATURE_RANGE = parameter_ranges["temperature"]
    if "salt_solute" in parameter_ranges:
        DROPLET_SALT_SOLUTE_LOG_RANGE = parameter_ranges["salt_solute"]
    if "air_temperature" in parameter_ranges:
        DROPLET_AIR_TEMPERATURE_RANGE = parameter_ranges["air_temperature"]
    if "relative_humidity" in parameter_ranges:
        DROPLET_RELATIVE_HUMIDITY_RANGE = parameter_ranges["relative_humidity"]
    if "rhoa" in parameter_ranges:
        DROPLET_RHOA_RANGE = parameter_ranges["rhoa"]
    if "time" in parameter_ranges:
        DROPLET_TIME_LOG_RANGE = parameter_ranges["time"]

def normalize_droplet_parameters( droplet_parameters, parameter_ranges={} ):
    """
    Normalizes an array of droplet parameters into the range [-1, 1].  Operates
    on both droplet inputs and background parameters.

    Droplet parameters may be specified either as NumPy arrays or as PyTorch
    tensors.

    Takes 2 arguments:

      droplet_parameters - Array of parameters to normalize, shaped
                           number_droplets x number_parameters.  May either be
                           1D, for a single droplet, or 2D for multiple droplets.
                           number_parameters must be either 2 (output parameters),
                           6 (input parameters without t_final), or 7 (input
                           parameters with t_final).  Must be either a NumPy
                           array or a PyTorch tensor.
      parameter_ranges   - Optional dictionary of parameter ranges to normalize
                           with.  These must be type compatible with
                           droplet_parameters as well as residing on the same
                           device when PyTorch tensors are provided.  If
                           omitted, the default parameter ranges are used and
                           are created to match the type and device placement of
                           droplet_parameters.

    Returns 1 value:

      normalized_droplet_parameters - Array of normalized parameters, shaped
                                      number_droplets x number_parameters.  Has the
                                      same type and shape as droplet_parameters.

    """

    # Figure out what type of array we were provided and get our ranges in
    # the right type and location.
    if isinstance( droplet_parameters, np.ndarray ):
        xp = np
    else:
        xp = torch

    # Get our parameter ranges by name, either as the caller provided or
    # creating tensors from our default ranges.
    if len( parameter_ranges ) > 0:
        #
        # NOTE: There is no verification these are actually in the right format
        #       or on the same device as droplet_parameters.  We assume the
        #       caller will do the right thing if they use the advanced
        #       interface...
        #
        radius_log_range        = parameter_ranges["radius"]
        temperature_range       = parameter_ranges["temperature"]
        salt_solute_log_range   = parameter_ranges["salt_solute"]
        air_temperature_range   = parameter_ranges["air_temperature"]
        relative_humidity_range = parameter_ranges["relative_humidity"]
        rhoa_range              = parameter_ranges["rhoa"]
    else:
        radius_log_range        = DROPLET_RADIUS_LOG_RANGE
        temperature_range       = DROPLET_TEMPERATURE_RANGE
        salt_solute_log_range   = DROPLET_SALT_SOLUTE_LOG_RANGE
        air_temperature_range   = DROPLET_AIR_TEMPERATURE_RANGE
        relative_humidity_range = DROPLET_RELATIVE_HUMIDITY_RANGE
        rhoa_range              = DROPLET_RHOA_RANGE

        # Do we need to convert the NumPy ranges to tensors?
        if xp == torch:
            # Make sure we move the converted tensors to the correct device.
            droplet_device = droplet_parameters.device

            radius_log_range        = torch.from_numpy( radius_log_range ).to( droplet_device )
            temperature_range       = torch.from_numpy( temperature_range ).to( droplet_device )
            salt_solute_log_range   = torch.from_numpy( salt_solute_log_range ).to( droplet_device )
            air_temperature_range   = torch.from_numpy( air_temperature_range ).to( droplet_device )
            relative_humidity_range = torch.from_numpy( relative_humidity_range ).to( droplet_device )
            rhoa_range              = torch.from_numpy( rhoa_range ).to( droplet_device )

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

    normalized_droplet_parameters = xp.empty_like( droplet_parameters )

    # We always have radius and temperature.
    normalized_droplet_parameters[..., 0] = (xp.log10( droplet_parameters[..., 0] ) - xp.mean( radius_log_range )) / (xp.diff( radius_log_range ) / 2)
    normalized_droplet_parameters[..., 1] = (droplet_parameters[..., 1] - xp.mean( temperature_range )) / (xp.diff( temperature_range ) / 2)

    # Sometimes we have the remaining parameters.
    if number_parameters > 2:
        normalized_droplet_parameters[..., 2] = (xp.log10( droplet_parameters[..., 2] ) - xp.mean( salt_solute_log_range )) / (xp.diff( salt_solute_log_range ) / 2)
        normalized_droplet_parameters[..., 3] = (droplet_parameters[..., 3] - xp.mean( air_temperature_range )) / (xp.diff( air_temperature_range ) / 2)
        normalized_droplet_parameters[..., 4] = (droplet_parameters[..., 4] - xp.mean( relative_humidity_range )) / (xp.diff( relative_humidity_range ) / 2)
        normalized_droplet_parameters[..., 5] = (droplet_parameters[..., 5] - xp.mean( rhoa_range )) / (xp.diff( rhoa_range ) / 2)

        # Copy the evaluation times when they're present.
        if number_parameters > 6:
            normalized_droplet_parameters[..., 6] = droplet_parameters[..., 6]

    return normalized_droplet_parameters

def scale_droplet_parameters( droplet_parameters, parameter_ranges={} ):
    """
    Scales an array of droplet parameters from the range [-1, 1] to their
    expected ranges of physical values.  Operates on both input and parameters.

    Droplet parameters may be specified either as NumPy arrays or as PyTorch
    tensors.

    Takes 2 arguments:

      droplet_parameters - Array of parameters to scale, shaped
                           number_droplets x number_parameters.  May either be
                           1D, for a single droplet, or 2D for multiple droplets.
                           number_parameters must be either 2 (output parameters),
                           6 (input parameters without t_final), or 7 (input
                           parameters with t_final).  Must be either a NumPy
                           array or a PyTorch tensor.
      parameter_ranges   - Optional dictionary of parameter ranges to scale
                           with.  These must be type compatible with
                           droplet_parameters as well as residing on the same
                           device when PyTorch tensors are provided.  If
                           omitted, the default parameter ranges are used and
                           are created to match the type and device placement of
                           droplet_parameters.

    Returns 1 value:

      scaled_droplet_parameters - Array of scaled parameters, shaped
                                  number_droplets x number_parameters.  Has the
                                  same type and shape as droplet_parameters.


    """

    # Figure out what type of array we were provided and get our ranges in
    # the right type and location.
    if isinstance( droplet_parameters, np.ndarray ):
        xp = np
    else:
        xp = torch

    # Get our parameter ranges by name, either as the caller provided or
    # creating tensors from our default ranges.
    if len( parameter_ranges ) > 0:
        #
        # NOTE: There is no verification these are actually in the right format
        #       or on the same device as droplet_parameters.  We assume the
        #       caller will do the right thing if they use the advanced
        #       interface...
        #
        radius_log_range        = parameter_ranges["radius"]
        temperature_range       = parameter_ranges["temperature"]
        salt_solute_log_range   = parameter_ranges["salt_solute"]
        air_temperature_range   = parameter_ranges["air_temperature"]
        relative_humidity_range = parameter_ranges["relative_humidity"]
        rhoa_range              = parameter_ranges["rhoa"]
    else:
        radius_log_range        = DROPLET_RADIUS_LOG_RANGE
        temperature_range       = DROPLET_TEMPERATURE_RANGE
        salt_solute_log_range   = DROPLET_SALT_SOLUTE_LOG_RANGE
        air_temperature_range   = DROPLET_AIR_TEMPERATURE_RANGE
        relative_humidity_range = DROPLET_RELATIVE_HUMIDITY_RANGE
        rhoa_range              = DROPLET_RHOA_RANGE

        # Do we need to convert the NumPy ranges to tensors?
        if xp == torch:
            # Make sure we move the converted tensors to the correct device.
            droplet_device = droplet_parameters.device

            radius_log_range        = torch.from_numpy( radius_log_range ).to( droplet_device )
            temperature_range       = torch.from_numpy( temperature_range ).to( droplet_device )
            salt_solute_log_range   = torch.from_numpy( salt_solute_log_range ).to( droplet_device )
            air_temperature_range   = torch.from_numpy( air_temperature_range ).to( droplet_device )
            relative_humidity_range = torch.from_numpy( relative_humidity_range ).to( droplet_device )
            rhoa_range              = torch.from_numpy( rhoa_range ).to( droplet_device )

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

    scaled_droplet_parameters = xp.empty_like( droplet_parameters )

    # We always have radius and temperature.
    scaled_droplet_parameters[..., 0] = 10.0 ** (droplet_parameters[..., 0] * (xp.diff( radius_log_range ) / 2) + xp.mean( radius_log_range ))
    scaled_droplet_parameters[..., 1] = droplet_parameters[..., 1] * (xp.diff( temperature_range ) / 2) + xp.mean( temperature_range )

    # Sometimes we have the remaining parameters.
    if number_parameters > 2:
        scaled_droplet_parameters[..., 2] = 10.0 ** (droplet_parameters[..., 2] * (xp.diff( salt_solute_log_range ) / 2) + xp.mean( salt_solute_log_range ))
        scaled_droplet_parameters[..., 3] = droplet_parameters[..., 3] * (xp.diff( air_temperature_range ) / 2) + xp.mean( air_temperature_range )
        scaled_droplet_parameters[..., 4] = droplet_parameters[..., 4] * (xp.diff( relative_humidity_range ) / 2) + xp.mean( relative_humidity_range )
        scaled_droplet_parameters[..., 5] = droplet_parameters[..., 5] * (xp.diff( rhoa_range ) / 2) + xp.mean( rhoa_range )

        # Copy the evaluation times when they're present.
        if number_parameters > 6:
            scaled_droplet_parameters[..., 6] = droplet_parameters[..., 6]

    return scaled_droplet_parameters

def solve_ivp_float32_outputs( dydt, t_span, y0, atol=BDF_TOLERANCE_ABSOLUTE, rtol=BDF_TOLERANCE_RELATIVE, **kwargs ):
    """
    Solves an initial value problem and returns the solutions in 32-bit precision.
    Error checking is present to ensure that failed solves provide a solution
    (i.e. NaN's) and have the same caller interface as successful solves.
    In the case of solver failures a warning is issued reporting the inputs and
    parameters to aide in debugging why the solve failed.

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
      atol   - Optional floating point scalar specifying the absolute tolerance
               to use during the solve.  If omitted, defaults to BDF_TOLERANCE_ABSOLUTE.
      rtol   - Optional floating point scalar specifying the relative tolerance
               to use during the solve.  If omitted, defaults to BDF_TOLERANCE_RELATIVE.
      kwargs - Optional keyword arguments to supply to solve_ivp().

    Returns 1 value:

      solution - NumPy array, shaped len( y0 ) x 1, containing the solution of
                 dydt at t_span[1] using y0 as the initial conditions.  Will
                 contain np.nan if solve_ivp() failed to provide a solution.

    """

    # Figure out how many points we're solving for.
    if "t_eval" in kwargs:
        number_outputs = len( kwargs["t_eval"] )
    else:
        number_outputs = 1

    solution = np.empty( (number_outputs, len( y0 ),), dtype=np.float32 )

    solve_failed_flag = False
    failure_message   = None

    # Solve the ODE in the precision supplied by the caller.
    try:
        ode_solution = solve_ivp( dydt, t_span, y0, atol=atol, rtol=rtol, **kwargs )

        if ode_solution.success and isinstance( ode_solution.y, np.ndarray ):
            # If we could integrate, return the outputs as the requested
            # precision.
            solution[:] = ode_solution.y.T.astype( "float32" )
        else:
            # Otherwise, capture useful information from the solution object to
            # help in debugging.
            solve_failed_flag = True
            failure_message   = "failed integration - message: {:s}, status: {:d}, number_rhs_evals: {:d}".format(
                ode_solution.message,
                ode_solution.status,
                ode_solution.nfev )
    except Exception as e:
        # Something went sideways, so report the exception.
        solve_failed_flag = True
        failure_message   = "exception - '{:s}'".format( str( e ) )

    # Failed solves return NaNs and trigger a warning with all of the
    # important details.
    if solve_failed_flag:
        solution[:, :] = np.nan

        warnings.warn( "solve_ivp() failed for "
                           "inputs=np.array( [{:.15g}, {:.15g}] ), "
                           "parameters=np.array( [{:.15g}, {:.15g}, {:.15g}, {:.15g}] ), "
                           "time={:15g}, kwargs={}\n{:s}.".format(
                               y0[0],
                               y0[1],
                               kwargs["args"][0][0],
                               kwargs["args"][0][1],
                               kwargs["args"][0][2],
                               kwargs["args"][0][3],
                               t_span[-1],
                               kwargs,
                               failure_message ) )

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

    NOTE: This is a wrapper around solve_ivp_float32_outputs(), so it takes the
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

      solution - NumPy array, shaped len( y0 ) x 1, containing the solution of
                 dydt at t_span[1] using y0 as the initial conditions.  Will
                 contain np.nan if solve_ivp() failed to provide a solution.

    """

    return solve_ivp_float32_outputs( *args, **kwargs )
