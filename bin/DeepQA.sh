#!/bin/sh
#DeepQA prediction file for protein single-model quality assessment #
if [ $# -ne 3 ]
then
	echo "need three parameters : path of fasta sequence, directory of input pdbs, directory of output"
	exit 1
fi
perl /home/jh7x3/tools/DeepQA/script/DeepQA.pl /home/jh7x3/tools/DeepQA/script /home/jh7x3/tools/DeepQA/tools $1 $2 $3
