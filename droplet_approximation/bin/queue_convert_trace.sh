#!/bin/sh

# Converts a single NTLP trace file using N processes in a job array.
# Instantiates a template and submits a job array task to the scheduler.

TESTING_FLAG="no"

if [ "$1" == "-t" ]; then
    TESTING_FLAG="yes"
    shift 1
fi

if [ $# != 5 ]; then
    echo "Usage: $0 <simulation_name> <particles_root> <number_processes> <trace_number> <trace_path>"
    exit 1
fi

# Pull the template parameters from our command line.
TRACE_NAME=$1
PARTICLES_ROOT=$2
NUMBER_NODES=$3
TRACE_NUMBER=$4
TRACE_PATH=$5

# Send mail to ourselves when the job starts and ends.
EMAIL="${USER}@nd.edu"

# Figure out where the job template lives.
INSTALL_DIR=`dirname $0`
CONVERT_TRACE_TEMPLATE=${INSTALL_DIR}/convert_trace.job_template

# Simply print the instantiated template instead of submitting the job
# when we're testing.
if [ "${TESTING_FLAG}" = "yes" ]; then
    QSUB=cat
else
    QSUB="qsub -V"
fi

# Change to the particles root we're populating to collect the job's output
# files.
cd ${PARTICLES_ROOT}

sed -e "s#_EMAIL_#${EMAIL}#g" \
    -e "s#_NUMBER_NODES_#${NUMBER_NODES}#g" \
    -e "s#_PARTICLES_ROOT_#${PARTICLES_ROOT}#g" \
    -e "s#_TRACE_NAME_#${TRACE_NAME}#g" \
    -e "s#_TRACE_NUMBER_#${TRACE_NUMBER}#g" \
    -e "s#_TRACE_PATH_#${TRACE_PATH}#g" \
    ${CONVERT_TRACE_TEMPLATE} | \
    ${QSUB}
