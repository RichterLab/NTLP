if [[ "$#" -lt 1  || ( "$2" -ne "" && "$2" -ne "-j") ]]; then
	echo "Usage: $0 <job_count> [<-j>]"
	exit 1
fi

job_count=$1

qsub -N generate_batch -t 1-$job_count parallel_generate_training_data.sh $2

wait

qsub -hold_jid generate_batch train.sh
