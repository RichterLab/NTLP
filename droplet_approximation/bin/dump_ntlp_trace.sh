#!/bin/sh

# Pretty prints one or more NTLP particle traces.  Outputs 9 columns of data:
#
#  1. Global particle index
#  2. Backward Euler status flag (zero is success, non-zero denotes the failure mode)
#  3. Simulation time (not delta t)
#  4. Input particle radius
#  5. Input particle temperature
#  6. Salt mass
#  7. Air temperature
#  8. Relative humidity
#  9. Air density
#
# Usage: dump_ntlp_trace.sh <trace_path> [...]
#

hexdump -e '2/4 "%10d " " " 7/4 "%.8g " "\n"'  $*
