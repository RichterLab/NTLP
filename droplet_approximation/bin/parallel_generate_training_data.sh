#!/bin/sh

# Creates training data in parallel using N-many instances at a time.  Data are
# created in separate files that can be concatenated after the fact
# (e.g. cat *.data > combined.data).

# Script that generates the training data.
GENERATOR=./generate_training_data.py

# Number of cores on the local system.  This ignores hardware threads that share
# a core.  Enumerate the "CPU" resource instead of "CORE" if hardware threads
# are desired.
NUMBER_SYSTEM_CORES=`(lscpu -p=CORE 2>/dev/null || echo 1) | grep -v ^# | sort -u | wc -l`

# Prefix of the files generated.  This can include paths so data are generated
# in a different directory.
OUTPUT_PREFIX=../data/time_log_spaced

# Total number of training data files to create.
NUMBER_ITERATIONS=320

# Number of files to create in parallel.  This should be no larger than the
# number of cores on the system.
NUMBER_JOBS=${NUMBER_SYSTEM_CORES}

# Create the directory where we're generating data if it doesn't already exist.
OUTPUT_DIR=`dirname ${OUTPUT_PREFIX}`
if [ ! -d ${OUTPUT_DIR} ]; then
    mkdir -p ${OUTPUT_DIR}
fi

# Loop over all of the iterations and generate data in batches.
ITERATION_NUMBER=0
while [ ${ITERATION_NUMBER} -lt ${NUMBER_ITERATIONS} ]; do
    START_NUMBER=${ITERATION_NUMBER}
    ITERATION_NUMBER=`expr ${ITERATION_NUMBER} + ${NUMBER_JOBS}`
    END_NUMBER=`expr ${ITERATION_NUMBER} - 1`

    # Kick off the next batch of generators in the background.
    for JOB_NUMBER in `seq -f %04g ${START_NUMBER} ${END_NUMBER}`; do
        OUTPUT_FILE_NAME=${OUTPUT_PREFIX}-${JOB_NUMBER}.data
        WEIRD_FILE_NAME=${OUTPUT_PREFIX}-${JOB_NUMBER}-weird.xlsx

        echo "Generating log data #${JOB_NUMBER}"
        ${GENERATOR} ${OUTPUT_FILE_NAME} ${WEIRD_FILE_NAME} &
    done

    # Wait for each of the generators to complete before going to the next
    # batch.
    wait
done
