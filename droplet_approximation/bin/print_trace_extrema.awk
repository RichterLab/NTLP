#!/usr/bin/awk -f

# Computes the extrema of NTLP particle trace records from a raw trace file.
# Expects each record's 9 values to be printed on a single line like the output
# from:
#
#   $ hexdump -e '2/4 "%10d " " " 7/4 "%.8g " "\n"' be_dump_0001.data
#
# Prints output like every 1,000,000 records:
#
#    Lines processed: 4000000
#    Column    Minimum           Maximum
#    ------    -------           -------
#         4    8.53125e-07       0.000230812
#         5    298.297           301.305
#         6    1.13935e-15       1.75123e-09
#         7    298.996           301.305
#         8    0.769815          0.988707
#         9    1                 1
#
# The columns printed are:
#
#  4. Input radii
#  5. Input temperature
#  6. Salt mass
#  7. Air temperature
#  8. Relative humidity (fractional)
#  9. Air density (fractional)
#
# A final summary if printed if the file does not contain exactly a multiple of
# 1,000,000 records:
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
# This script can be run on a single trace file or all trace files in parallel
# like so:
#
#   $ find . -name 'be_dump_*.data' | \
#         parallel -j 64 'echo "Processing: {/}"; \
#                        '"${DUMP_TRACE}"' {} | \
#                        '"${AWK_SCRIPT}"' > '"${OUTPUT_DIR}"'/{/.}.results; \
#                        echo "Completed: {/}"' &
#
# With DUMP_TRACE set to the command that generates the hexdump output described
# above and this script set to AWK_SCRIPT.

BEGIN {
    # Initialize our extrema arrays for each parameter tracked.
    #
    # NOTE: We use empty strings as sentinels instead of large positive/negative
    #       values since +-Inf doesn't seem to work.
    #
    for (i = 4; i <= 9; i++) {
        min[i] = ""
        max[i] = ""
    }
    line_count = 0
}

{
    line_count++

    # Process columns 4-9 (fields $4 through $9)
    for (i = 4; i <= 9; i++) {
        val = $i + 0  # Convert to numeric

        # Initialize min/max on first line
        if (min[i] == "") {
            min[i] = val
            max[i] = val
        } else {
            # Update min/max
            if (val < min[i]) min[i] = val
            if (val > max[i]) max[i] = val
        }
    }

    # Print results every 1,000,000 lines
    if (line_count % 1000000 == 0) {
        printf "Lines processed: %d\n", line_count
        printf "Column    Minimum           Maximum\n"
        printf "------    -------           -------\n"
        for (i = 4; i <= 9; i++) {
            printf "%6d    %-15g   %-15g\n", i, min[i], max[i]
        }
        printf "\n"
    }
}

END {
    # Print final results if we didn't end on a 10,000 line boundary
    if (line_count % 1000000 != 0) {
        printf "Final results after %d lines:\n", line_count
        printf "Column    Minimum           Maximum\n"
        printf "------    -------           -------\n"
        for (i = 4; i <= 9; i++) {
            printf "%6d    %-15g   %-15g\n", i, min[i], max[i]
        }
    }
}
