#!/bin/sh

# Creates training data in parallel using N-many instances at a time.  Data are
# created in separate files that can be concatenated after the fact
# (e.g. cat *.data > combined.data).

# Script that generates the training data.
GENERATOR=./generate_training_data.py

# Number of cores on the local system.  This ignores hardware threads that share
# a core.  Enumerate the "CPU" resource instead of "CORE" if hardware threads
# are desired.
#NUMBER_SYSTEM_CORES=`(lscpu -p=CORE 2>/dev/null || echo 1) | grep -v ^# | sort -u | wc -l`
NUMBER_SYSTEM_CORES=32

# Prefix of the files generated.  This can include paths so data are generated
# in a different directory.
OUTPUT_PREFIX=./test

# Total number of training data files to create.
#NUMBER_ITERATIONS=320
NUMBER_ITERATIONS=1

# Number of files to create in parallel.  This should be no larger than the
# number of cores on the system.
NUMBER_JOBS=${NUMBER_SYSTEM_CORES}

# Create the directory where we're generating data if it doesn't already exist.
OUTPUT_DIR=`dirname ${OUTPUT_PREFIX}`
if [ ! -d ${OUTPUT_DIR} ]; then
    mkdir -p ${OUTPUT_DIR}
fi


    START_NUMBER=1
    END_NUMBER=${NUMBER_SYSTEM_CORES}

    for JOB_NUMBER in `seq -f %04g ${START_NUMBER} ${END_NUMBER}`; do
        OUTPUT_FILE_NAME=${OUTPUT_PREFIX}-${JOB_NUMBER}.data

        #echo "Generating log data #${JOB_NUMBER}"
        echo "Generating file ${OUTPUT_FILE_NAME}"
        ${GENERATOR} ${OUTPUT_FILE_NAME} &
    done

#    # Wait for each of the generators to complete before going to the next
#    # batch.
    wait
#done
