if [[ "$#" -ne 2 ]]; then
	echo "Usage: $0 <job_count>"
	exit 1
fi

job_count=$1

qsub -N generate_batch -t 1-$job_count parallel_generate_training_data.sh

wait

qsub -hold_jid generate_batch train.sh
