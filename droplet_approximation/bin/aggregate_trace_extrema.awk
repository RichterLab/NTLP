#!/usr/bin/awk -f

# Computes extrema of NTLP particle trace records from extrema computed over
# subsets of records (i.e. individual trace files as they're being written).
# This allows reporting of parameter ranges from the raw data without 1)
# waiting for the simulation to complete or 2) post-processing into a different
# format.
#
# Intended to be run like so:
#
#  $ tail -q -n9 be_dump_*.result | ./aggregate_trace_extrema.awk
#
# Where the last 9 lines of each .result file are of the format:
#
#   Final results after 50956369 lines:
#   Column    Minimum           Maximum
#   ------    -------           -------
#        4    8.53217e-07       0.000287857
#        5    297.164           301.321
#        6    1.13935e-15       3.397e-09
#        7    297.215           301.356
#        8    0.769865          0.982808
#        9    1                 1
#
# Which are the ranges of values found in NTLP trace records, as printed out by
# print_trace_extrema.awk.  Columns are:
#
#  4. Input radii
#  5. Input temperature
#  6. Salt mass
#  7. Air temperature
#  8. Relative humidity (fractional)
#  9. Air density (fractional)
#

BEGIN {
    # Initialize global extrema arrays.
    #
    # NOTE: We use empty strings as sentinels instead of large positive/negative
    #       values since +-Inf doesn't seem to work.
    #
    for (i = 4; i <= 9; i++) {
        global_min[i] = ""
        global_max[i] = ""
    }
    total_lines = 0
    file_count = 0
}

# Process lines that contain column data (numeric column number followed by values)
/^[ ]*[4-9][ ]+/ {
    col = $1 + 0  # Column number
    min_val = $2 + 0  # Minimum value
    max_val = $3 + 0  # Maximum value

    # Update global extrema
    if (global_min[col] == "" || min_val < global_min[col]) {
        global_min[col] = min_val
    }
    if (global_max[col] == "" || max_val > global_max[col]) {
        global_max[col] = max_val
    }
}

# Extract total lines processed from each file
/Final results after [0-9]+ lines:/ {
    file_count++
    lines = $4 + 0
    total_lines += lines
}

END {
    printf "=== GLOBAL EXTREMA ACROSS ALL %d FILES ===\n", file_count
    printf "Total lines processed: %d\n\n", total_lines
    printf "Column    Global Minimum    Global Maximum\n"
    printf "------    --------------    --------------\n"

    for (i = 4; i <= 9; i++) {
        if (global_min[i] != "" && global_max[i] != "") {
            printf "%6d    %-15g   %-15g\n", i, global_min[i], global_max[i]
        }
    }
}
