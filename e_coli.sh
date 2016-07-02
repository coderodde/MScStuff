#!/bin/bash 
ROOT_DIR=/cs/fs/home/tomescu/complete-contigs

### CONFIG ##############################################
DATA_DIR=/cs/fs/home/tomescu/complete-contigs/data/e_coli
INPUT_FILE=$DATA_DIR/e_coli.K12.fa
TYPE=circular
### CONFIG ##############################################

$ROOT_DIR/aligner/aligner --buildindex --$TYPE --reference $INPUT_FILE
for K in 31
do
	echo "************************** $INPUT_FILE $TYPE k $K **************************"
	$ROOT_DIR/bin/complete-contigs -i $INPUT_FILE -k $K -t 8 -g $TYPE

	# Removing non-maximal omnitigs
	$ROOT_DIR/bin/maximality $INPUT_FILE.k$K.$TYPE.omnitigs

	# Getting the stats
	$ROOT_DIR/aligner/aligner --loadindex --$TYPE --reference $INPUT_FILE --readfile $INPUT_FILE.k$K.$TYPE --outputfile $INPUT_FILE.k$K.$TYPE.stats.csv

	# Cleaning up
	rm $INPUT_FILE.k$K.$TYPE.omnitigs.rlcsa.array
	rm $INPUT_FILE.k$K.$TYPE.omnitigs.rlcsa.parameters
	#
	# rm $INPUT_FILE.k$K.$TYPE.unitigs
	rm $INPUT_FILE.k$K.$TYPE.unitigs.aligned
	# rm $INPUT_FILE.k$K.$TYPE.non-switching-contigs
	rm $INPUT_FILE.k$K.$TYPE.non-switching-contigs.aligned
	# rm $INPUT_FILE.k$K.$TYPE.YtoV-contigs
	rm $INPUT_FILE.k$K.$TYPE.YtoV-contigs.aligned
	rm $INPUT_FILE.k$K.$TYPE.omnitigs
	# rm $INPUT_FILE.k$K.$TYPE.omnitigs.maximal
	rm $INPUT_FILE.k$K.$TYPE.omnitigs.maximal.aligned
	rm $INPUT_FILE.k$K.$TYPE.safepairs
done

# # Cleaning up the index
# rm $INPUT_FILE.rlcsa.array
# rm $INPUT_FILE.rlcsa.parameters
# rm $INPUT_FILE.rlcsa.sa_samples