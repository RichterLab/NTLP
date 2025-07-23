#!/bin/bash

# Evaluates all of the particles in a directory hierarchy using multiple
# instances of Evaluate Particles each working on a subset of the particles.
# Each instance is submitted to the scheduler for batch processing.
#
# NOTE: This currently expects the particle evaluation script to live in the
#       same directory as this script!
#
USAGE_STR="$0 <simulation_name> <particles_root> <job_count> [<testing_flag>]"

# Parse the command line arguments.
SIMULATION_NAME=${1?${USAGE_STR}}
PARTICLES_ROOT=${2?${USAGE_STR}}
JOB_COUNT=${3?${USAGE_STR}}

if [ $# -gt 3 ]; then
    TESTING_FLAG="yes"
else
    TESTING_FLAG="no"
fi

# Takes a count and the number of chunks to break the count up into.  Echoes an
# array of (start, end) indices suitable for use in Python (open upper bounds)
# to standard output.  When the count isn't evenly divisible by the number of
# chunks the remainder is distributed across the initial chunks.
partition_indices()
{
    local COUNT=$1
    local NUMBER_CHUNKS=$2

    local CHUNK_SIZE=$((COUNT / NUMBER_CHUNKS))
    local REMAINDER=$((COUNT % NUMBER_CHUNKS))

    local START_INDEX=0
    local INDICES=()

    for ((CHUNK_INDEX = 0; CHUNK_INDEX < NUMBER_CHUNKS; CHUNK_INDEX++)); do
        local END_INDEX=$((START_INDEX + CHUNK_SIZE + (CHUNK_INDEX < REMAINDER ? 1 : 0)))

        INDICES+=("${START_INDEX}:${END_INDEX}")

        START_INDEX=${END_INDEX}
    done

    echo "${INDICES[@]}"
}

# Job parameters that are independent of the simulation and job chunk.  These
# specify:
#
#  - Where to send completion emails to (-M)
#  - Which nodes (-pe) and queue (-q) should be used
#  - Merging standard output and error together (-j)
#
# NOTE: We shouldn't quote things as some qsub options have multiple values
#       despite being whitespace-delimeted.
#
SCHEDULER_ARGUMENTS="-M ${USER}@nd.edu \
-pe mpi-64 64 \
-q *@@richter \
-j y
"

# This script's installation directory.
SCRIPTS_DIRECTORY=`dirname $0`

# Path to the evaluation script.
EVALUATE_PARTICLES="${SCRIPTS_DIRECTORY}/evaluate_particles.py"

# Specify the number of processes to use for each job.
#
# NOTE: This is hardcoded to the AMD Epyc cluster's nodes' core count.  We don't
#       compute this via lscpu as we're running on a head node which isn't
#       representative of the cluster's nodes (48 vs 64 cores).
#
NUMBER_CORES=64

# Path to the particles index.
PARTICLES_INDEX="${PARTICLES_ROOT}/particles.index"
if [ ! -d ${PARTICLES_ROOT} ]; then
    echo "'${PARTICLES_ROOT}' does not exist!" >&2
    exit 1
elif [ ! -f ${PARTICLES_INDEX} ]; then
    echo "'${PARTICLES_INDEX}' does not exist!" >&2
    exit 2
fi

# Get the number of particles from the index's size.
NUMBER_PARTICLES=`stat -c %s ${PARTICLES_INDEX} | awk '{ print $1 / 4}'`

# Compute the individual indices per job we're submitting.
INDEX_ARRAY=($(partition_indices ${NUMBER_PARTICLES} ${JOB_COUNT}))

# Report what we're doing.
echo "Evaluating ${NUMBER_PARTICLES} particles in '${PARTICLES_ROOT}' with ${JOB_COUNT} jobs."

JOB_NUMBER=1
for INDEX in "${INDEX_ARRAY[@]}"; do
    # Provide the job's command line on standard input.  Have the scheduler
    # maintain the user's submission environment (-V) for when the job is
    # launched and name it in a way that we can identify it later (-N).
    SCHEDULE_COMMAND="echo ${EVALUATE_PARTICLES} \
                              ${PARTICLES_ROOT} \
                              ${INDEX} \
                              ${NUMBER_CORES} | \
                      qsub \
                          -V \
                          ${SCHEDULER_ARGUMENTS} \
                          -N ${SIMULATION_NAME}-${JOB_NUMBER}"

    # Print or execute based on whether we're testing.
    if [ "${TESTING_FLAG}" = "yes" ]; then
        echo ${SCHEDULE_COMMAND}
    else
        eval ${SCHEDULE_COMMAND}
    fi

    JOB_NUMBER=$((JOB_NUMBER + 1))
done
