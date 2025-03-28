#!/usr/bin/env python3

#
# NOTE: This is a snapshot of code in the Jupyter Notebook from late August
#       2024.  It should be factored so both the notebook and this script use
#       the same functions.  Be *very* careful when making updates to one or
#       the other.
#

# Usage: generate_training_data.py <output_path> <spreadsheet_path>
#
# Creates 1,048,576 (1024^2) droplet parameters and writes them to <output_path>.
# Any weird parameters are logged to <spreadsheet_path> and ignored so they are
# not written to <output_path>.

import sys

import numpy as np
from scipy.integrate import solve_ivp

import errno
import os
import signal
import functools

# Config constants
NUMBER_DROPLETS = 1024*1024*5

DROPLET_RADIUS_LOG_RANGE          = np.array( (-8, -3) )
DROPLET_TEMPERATURE_RANGE         = np.array( (273, 310) )
DROPLET_SALINITY_LOG_RANGE        = np.array( (-22, -10) )
DROPLET_AIR_TEMPERATURE_RANGE     = np.array( (273, 310) )
DROPLET_RELATIVE_HUMIDITY_RANGE   = np.array( (0.65, 1.1) )
DROPLET_RHOA_RANGE                = np.array( (0.8, 1.2) )

def dydt( t, y, parameters ):

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

def solve_ivp_float32_outputs( dydt, t_span, y0, **kwargs ):
    """
    Solves the initial value problem in the requested precision and returns
    them in 32-bit precision.
    """

    # Solve the ODE in the precision supplied by the caller.
    solution = solve_ivp( dydt, t_span, y0, **kwargs )

    # Return the outputs as the requested precision.
    solution.t = solution.t.astype( "float32" )
    solution.y = solution.y.astype( "float32" )

    return solution

# Notice that this function is never called!
def normalize_droplet_parameters( droplet_parameters ):
    """
    """

    #print( "XXX: Handle the optional log/exponentiation, and in scale() as well" )

    #
    # NOTE: We take care to handle arbitrary rank arrays!  While we expect either
    #       rank-1 or rank-2, the code is written to handle larger ranks naturally
    #       with the inner dimension representing a single droplet's parameters.
    #

    # Grabs the last dimension of the array
    number_parameters = droplet_parameters.shape[-1]

    if number_parameters != 2 and number_parameters != 6: # Either just temp/radius or full model?
        raise ValueError( "Unknown number of parameters to normalize ({:d})!".format( number_parameters ) )

    normalized_droplet_parameters = np.empty_like( droplet_parameters )

    # We always have radius and temperature.
    normalized_droplet_parameters[..., 0] = (np.log10( droplet_parameters[..., 0] ) - np.mean( DROPLET_RADIUS_LOG_RANGE )) / (np.diff( DROPLET_RADIUS_LOG_RANGE ) / 2)
    normalized_droplet_parameters[..., 1] = (droplet_parameters[..., 1] - np.mean( DROPLET_TEMPERATURE_RANGE )) / (np.diff( DROPLET_TEMPERATURE_RANGE ) / 2)

    # Sometimes we have the remaining parameters.
    if number_parameters > 2:
        normalized_droplet_parameters[..., 2] = (np.log10( droplet_parameters[..., 2] ) - np.mean( DROPLET_SALINITY_LOG_RANGE )) / (np.diff( DROPLET_SALINITY_LOG_RANGE ) / 2)
        normalized_droplet_parameters[..., 3] = (droplet_parameters[..., 3] - np.mean( DROPLET_AIR_TEMPERATURE_RANGE )) / (np.diff( DROPLET_AIR_TEMPERATURE_RANGE ) / 2)
        normalized_droplet_parameters[..., 4] = (droplet_parameters[..., 4] - np.mean( DROPLET_RELATIVE_HUMIDITY_RANGE )) / (np.diff( DROPLET_RELATIVE_HUMIDITY_RANGE ) / 2)
        normalized_droplet_parameters[..., 5] = (droplet_parameters[..., 5] - np.mean( DROPLET_RHOA_RANGE )) / (np.diff( DROPLET_RHOA_RANGE ) / 2)

    return normalized_droplet_parameters

def scale_droplet_parameters( droplet_parameters ):
    """
    """

    #
    # NOTE: We take care to handle arbitrary rank arrays!  While we expect either
    #       rank-1 or rank-2, the code is written to handle larger ranks naturally
    #       with the inner dimension representing a single droplet's parameters.
    #

    number_parameters = droplet_parameters.shape[-1]

    if number_parameters != 2 and number_parameters != 6:
        raise ValueError( "Unknown number of parameters to scale ({:d})!".format( number_parameters ) )

    scaled_droplet_parameters = np.empty_like( droplet_parameters )

    # We always have radius and temperature.
    scaled_droplet_parameters[..., 0] = 10.0 ** (droplet_parameters[..., 0] * (np.diff( DROPLET_RADIUS_LOG_RANGE ) / 2) + np.mean( DROPLET_RADIUS_LOG_RANGE ))
    scaled_droplet_parameters[..., 1] = droplet_parameters[..., 1] * (np.diff( DROPLET_TEMPERATURE_RANGE ) / 2) + np.mean( DROPLET_TEMPERATURE_RANGE )

    # Sometimes we have the remaining parameters.
    if number_parameters > 2:
        scaled_droplet_parameters[..., 2] = 10.0 ** (droplet_parameters[..., 2] * (np.diff( DROPLET_SALINITY_LOG_RANGE ) / 2) + np.mean( DROPLET_SALINITY_LOG_RANGE ))
        scaled_droplet_parameters[..., 3] = droplet_parameters[..., 3] * (np.diff( DROPLET_AIR_TEMPERATURE_RANGE ) / 2) + np.mean( DROPLET_AIR_TEMPERATURE_RANGE )
        scaled_droplet_parameters[..., 4] = droplet_parameters[..., 4] * (np.diff( DROPLET_RELATIVE_HUMIDITY_RANGE ) / 2) + np.mean( DROPLET_RELATIVE_HUMIDITY_RANGE )
        scaled_droplet_parameters[..., 5] = droplet_parameters[..., 5] * (np.diff( DROPLET_RHOA_RANGE ) / 2) + np.mean( DROPLET_RHOA_RANGE )

    return scaled_droplet_parameters

def validate_output_parameters( input_parameters, output_parameters, integration_times ):
    """
    """

    def split_inputs_and_outputs( input_parameters, output_parameters, integration_times, droplet_index ):
        """
        """

        return [input_parameters[droplet_index, :2],
                input_parameters[droplet_index, 2:],
                integration_times[droplet_index],
                output_parameters[droplet_index, :]]

    #
    # NOTE: We take care to handle arbitrary rank arrays!  While we expect either
    #       rank-1 or rank-2, the code is written to handle larger ranks naturally
    #       with the inner dimension representing a single droplet's parameters.
    #

    number_input_parameters  = input_parameters.shape[-1]
    number_output_parameters = output_parameters.shape[-1]
    input_parameters         = input_parameters.reshape( -1, number_input_parameters )
    output_parameters        = output_parameters.reshape( -1, number_output_parameters )
    number_droplets          = input_parameters.shape[0]

    invalid_outputs = {
        "radius_nonpositive":    [],
        "radius_too_small":      [],
        "radius_too_large":      [],
        "temperature_too_small": [],
        "temperature_too_large": []
    }

    if input_parameters.shape[0] != output_parameters.shape[0]:
        raise ValueError( "Must have an equal number of input and output parameters ({:d} != {:d})!".format(
            input_parameters.shape[0],
            output_parameters.shape[0] ) )

    if number_input_parameters != 6 and number_output_parameters != 2 :
        raise ValueError( "Invalid number of parameters to validate for inputs and outputs ({:d}, {:d})!".format(
            number_input_parameters,
            number_output_parameters ) )

    # Iterate through each of the droplets and verify that the values are inside
    # the expected ranges.
    for droplet_index in np.arange( number_droplets ):

        inputs_and_outputs = split_inputs_and_outputs( input_parameters,
                                                       output_parameters,
                                                       integration_times,
                                                       droplet_index )

        # Ensure that we have an output radius inside the expected range.
        #
        # NOTE: We don't check for non-positive values here as this is a post-generation
        #       validation.  Non-positive values are detected and logged elsewhere.
        #
        if output_parameters[droplet_index, 0] < 10.0**DROPLET_RADIUS_LOG_RANGE[0]:
            invalid_outputs["radius_too_small"].append( inputs_and_outputs )
        elif output_parameters[droplet_index, 0] > 10.0**DROPLET_RADIUS_LOG_RANGE[1]:
            invalid_outputs["radius_too_large"].append( inputs_and_outputs )

        # Ensure the output temperature stays within the range.
        if output_parameters[droplet_index, 1] < DROPLET_TEMPERATURE_RANGE[0]:
            invalid_outputs["temperature_too_small"].append( inputs_and_outputs )
        elif output_parameters[droplet_index, 1] > DROPLET_TEMPERATURE_RANGE[1]:
            invalid_outputs["temperature_too_large"].append( inputs_and_outputs )

    return invalid_outputs


class TimeoutError(Exception):
    pass

def timeout(seconds=10, error_message=os.strerror(errno.ETIME)):
    def decorator(func):
        def _handle_timeout(signum, frame):
            raise TimeoutError(error_message)

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            signal.signal(signal.SIGALRM, _handle_timeout)
            signal.alarm(seconds)
            try:
                result = func(*args, **kwargs)
            finally:
                signal.alarm(0)
            return result

        return wrapper

    return decorator

@timeout(5)
def timed_solve_ivp(*args, **kwargs):
    return solve_ivp(*args, **kwargs)


def create_droplet_batch( number_droplets ):
    """
    Creates a batch of random droplets.

    Takes arguments:

      number_droplets - Number of droplets to generate parameters for.

    Returns 4 values:

      random_inputs     - Array, sized number_droplets x 6, containing the droplets
                          radii, temperatures, salinity, air temperature, relative
                          humidity, and rhoa.
      random_outputs    - Array, sized number_droplets x 2, containing the droplets
                          radii and temperatures.
      integration_times - Array, sized number_droplets, containing the times corresponding
                          to the associated random_inputs and random_outputs.
      weird_indices     - List of indices into random_outputs indicating NaN
                          values.  This will be empty when all outputs are valid.

    """

    # XXX: move this
    import warnings
    import time

    warnings.simplefilter( "error", RuntimeWarning )

    random_inputs     = scale_droplet_parameters( np.reshape( np.random.uniform( -1, 1, number_droplets*6 ),
                                                             (number_droplets, 6) ).astype( "float32" ) )
    random_outputs    = np.empty_like( random_inputs, shape=(number_droplets, 2) )



    integration_times = (10.0**np.random.uniform( -3.2, 1.1, number_droplets )).astype( "float32" )

    # Hard code problematic values

    number_warnings = 0

    weird_inputs = { "evaluation_warning": [],
                     "failed_solve":       [],
                     "invalid_inputs":     [],
                     "timeout_inputs":     [] }
    weird_outputs = { "radius_nonpositive":      [],
                      "temperature_nonpositive": [] }

    for droplet_index in np.arange( number_droplets ):
        nudge_count = 0
        while True:
            y0         = random_inputs[droplet_index, :2] # Radius, Temperature
            parameters = random_inputs[droplet_index, 2:] # Environemntal Variables
            t_final    = integration_times[droplet_index] # Time to integrate

            #y0 = np.array([5.9981421145494096e-06, 277.07498168945312])
            #parameters = np.array([3.8935551476971464e-22, 307.10308837890625, 0.7068866491317749, 0.93063300848007202])
            #t_final = 10.0
            # ^^ Presumably problematic hard-coded values

            #print(f"t_final = {t_final:.17g}")

            try:
                #print(f"droplet_index = ",droplet_index)
                #print("y0 =  [" + ", ".join(f"{x:.17g}" for x in y0) + "]")
                #print("parameters = [" + ", ".join(f"{x:.17g}" for x in parameters) + "]")
                start_time = time.time()

                solution = timed_solve_ivp( dydt, [0, t_final], y0, method="BDF", t_eval=[t_final], args=(parameters,) )

                #solution = solve_ivp( dydt, [0, t_final], y0, method="Radau", t_eval=[t_final], args=(parameters,) )
                end_time = time.time()
                elapsed_time = end_time - start_time
                #if elapsed_time > 1.0:
                    #print(f"LONG TIME",elapsed_time,y0,parameters,t_final)
                    #print(f"LONG TIME",elapsed_time)
                    #print("LONG TIME 1, [" + ", ".join(f"{x:.17g}" for x in y0) + "]")
                    #print("LONG TIME 2, [" + ", ".join(f"{x:.17g}" for x in parameters) + "]")
                if solution.success:
                    random_outputs[droplet_index, 0] = solution.y[0][0]
                    random_outputs[droplet_index, 1] = solution.y[1][0]

                    # Check that we didn't get a physically impossible solution.  These
                    # will be logged and we'll reroll the dice to replace them.
                    #
                    # NOTE: We don't strictly need to check for negative temperatures as
                    #       that will get covered in validate_output_parameters() but
                    #       we leave it here so it is easy to identify physically impossible
                    #       cases.
                    #
                    if solution.y[0][0] <= 0.0:
                        weird_outputs["radius_nonpositive"].append( [y0, parameters, t_final, random_outputs[droplet_index, :]] )
                    elif solution.y[1][0] <= 0.0:
                        weird_outputs["temperature_nonpositive"].append( [y0, parameters, t_final, random_outputs[droplet_index, :]] )
                    else:
                        break
                else:
                    # Record the cases that fail to converge.
                    weird_inputs["failed_solve"].append( [y0, parameters, t_final, solution.message] )

            except RuntimeWarning as e:
                #print("Warning message:",e)
                #print("RuntimeWarning:, [" + ", ".join(f"{x:.17g}" for x in y0) + "]")
                weird_inputs["evaluation_warning"].append( [y0, parameters, t_final, str( e )] )
            except ValueError as e:
                #print("ValueError:, [" + ", ".join(f"{x:.17g}" for x in y0) + "]")
                weird_inputs["invalid_inputs"].append( [y0, parameters, t_final, str( e )] )
            except TimeoutError as e:
                #print("TimeoutError:, [" + ", ".join(f"{x:.17g}" for x in y0) + "]")
                weird_inputs["timeout_inputs"].append( [y0, parameters, t_final, str( e )] )
                nudge_count += 1
                if (nudge_count < 3):
                    # Try again with slightly different values up to three times
                    #print("OLD m_s", random_inputs[droplet_index, :][2])
                    random_inputs[droplet_index, :][2] *= 1.0 + (0.0001*np.random.choice([-1,1]))# nudge value slightly
                    #print("NEW m_s", random_inputs[droplet_index, :][2])
                    continue
                nudge_count = 0
                # Otherwise, completely reroll

            # What if we did the same for long runtimes.
            # Reroll the dice but only slightly. Thoughts?

            # We failed to create acceptable parameters.  Reroll the dice for
            # this droplet and try again.
            random_inputs[droplet_index, :] = scale_droplet_parameters(
                np.random.uniform( -1, 1, 6 ).astype( "float32" ) )

    # Sanity check that we created good parameters and integrate any issues seen
    # into our current dictionary.
    weird_outputs.update( validate_output_parameters( random_inputs,
                                                      random_outputs,
                                                      integration_times ) )

    # Warn users if there were strange droplets.  XXX

    return random_inputs, random_outputs, integration_times, weird_inputs, weird_outputs

def merge_weird_parameters( parameters_1, parameters_2 ):
    """
    """

    merged_parameters = parameters_1.copy()

    for type_name, weird_things in parameters_2.items():
        if type_name in merged_parameters:
            merged_parameters[type_name].extend( weird_things )
        else:
            merged_parameters[type_name] = weird_things

    return merged_parameters

def write_weird_parameters_to_spreadsheet( file_name, weird_inputs, weird_outputs ):
    """
    Creates an Excel file containing one spreadsheet per type of weird inputs or
    outputs.

    Takes 3 arguments:

      file_name     -
      weird_inputs  -
      weird_outputs -

    Returns nothing.

    """

    # XXX: move this
    import pandas as pd

    # Open the file for writing, overwriting an existing file if necessary.
    with pd.ExcelWriter( file_name ) as writer:

        # Create a DataFrame for each category of input parameters weirdness.
        for type_name, weird_things in weird_inputs.items():
            sheet_name = "inputs-{:s}".format( type_name )

            # Columns in the DataFrame.
            columns = ["initial_radius",
                       "initial_temperature",
                       "salinity",
                       "air_temperature",
                       "relative_humidity",
                       "rhoa",
                       "t_final",
                       "error"]

            # Convert our list of lists of NumPy arrays to lists of lists
            # so spreadsheet cells contain a single value.
            weird_things_listified = []
            for weird_thing in weird_things:
                weird_thing_listified = []

                weird_thing_listified.extend( weird_thing[0].tolist() )
                weird_thing_listified.extend( weird_thing[1].tolist() )
                weird_thing_listified.extend( weird_thing[2:] )

                weird_things_listified.append( weird_thing_listified )

            # Build a DataFrame and write it as a new sheet.
            filtered_df = pd.DataFrame( weird_things_listified, columns=columns )
            filtered_df.to_excel( writer, sheet_name=sheet_name, index=False )

        # Create a DataFrame for each category of output parameters weirdness.
        for type_name, weird_things in weird_outputs.items():
            sheet_name = "outputs-{:s}".format( type_name )

            # Columns in the DataFrame.
            columns = ["initial_radius",
                       "initial_temperature",
                       "salinity",
                       "air_temperature",
                       "relative_humidity",
                       "rhoa",
                       "t_final",
                       "radius",
                       "temperature"]

            # Convert our list of lists of NumPy arrays to lists of lists
            # so spreadsheet cells contain a single value.
            weird_things_listified = []
            for weird_thing in weird_things:
                weird_thing_listified = []

                weird_thing_listified.extend( weird_thing[0].tolist() )
                weird_thing_listified.extend( weird_thing[1].tolist() )
                weird_thing_listified.append( weird_thing[2] )
                weird_thing_listified.extend( weird_thing[3].tolist() )

                weird_things_listified.append( weird_thing_listified )

            # Build a DataFrame and write it as a new sheet.
            filtered_df = pd.DataFrame( weird_things_listified, columns=columns )
            filtered_df.to_excel( writer, sheet_name=sheet_name, index=False )

def create_training_file( file_name, number_droplets, weird_file_name=None, user_batch_size=None ):
    """
    """

    # Balance the time it takes to generate a single batch vs the efficiency of
    # writing it out.  Each droplet's parameters takes 36 bytes.
    if user_batch_size is not None:
        BATCH_SIZE = user_batch_size
    else:
        BATCH_SIZE = 1024 * 10

    number_batches = (number_droplets + BATCH_SIZE - 1) // BATCH_SIZE

    weird_inputs  = {}
    weird_outputs = {}

    with open( file_name, "wb" ) as output_fp:
        # XXX: factor this into a private function
        for batch_index in range( number_batches ):
            print(f"Currently at {batch_index*BATCH_SIZE}")

            if batch_index != (number_batches - 1):
                batch_size = BATCH_SIZE
            else:
                batch_size = number_droplets % BATCH_SIZE

            # Get the next batch of droplets.
            (inputs,
             outputs,
             times,
             batch_weird_inputs,
             batch_weird_outputs) = create_droplet_batch( batch_size )

            # Track the weirdness for post-mortem analysis.
            weird_inputs  = merge_weird_parameters( weird_inputs, batch_weird_inputs )
            weird_outputs = merge_weird_parameters( weird_outputs, batch_weird_outputs )

            inputs_outputs = np.hstack( (inputs,
                                         times.reshape( (batch_size, 1) ),
                                         outputs) )

            inputs_outputs.tofile( output_fp )

    if weird_file_name is not None:
        write_weird_parameters_to_spreadsheet( weird_file_name,
                                               weird_inputs,
                                               weird_outputs )

def main( argv ):
    """
    """

    weird_file_name = None

    if len( argv ) == 1:
        print( "Usage: <output_file> [<weird_file>]" )
        sys.exit( 1 )
    if len( argv ) > 2:
        weird_file_name = argv[2]

    create_training_file( argv[1],
                          NUMBER_DROPLETS,
                          weird_file_name=weird_file_name )

if __name__ == "__main__":
    sys.exit( main( sys.argv ) )
