# create a unique id for this run, using the date command
# for MacOS, use gdate instead 
if hash gdate 2>/dev/null; then
    UNIQ_ID=$(gdate +%s%N | head -c 16)
else
    UNIQ_ID=$(date +%s%N | head -c 16)
fi

if [ "$#" -lt 3 ]; then
	
	echo "Usage: "
	echo "./run_block.sh   weights_file"
	echo "                 num_of_samples"
	echo "                 sequences_file"
	echo "                 [symmetric, width, max wins parameters]"
	exit 1
fi


INPUT_FILE=$1
NUM_OF_SAMPLES=$2
SEQUENCES_FILE=$3


PROGRAM_NAME=build/mrnai-sampler

mkdir logs
mkdir logs/"$UNIQ_ID"

echo "Rand Start :: $UNIQ_ID" > logs/"$UNIQ_ID"/info
echo ""
echo "UNIQ_ID: $UNIQ_ID"
echo ""

echo "Executing program: $PROGRAM_NAME -s $NUM_OF_SAMPLES -f structs_bip.out " "${@:4}"


echo "Sampling..."
cat "$INPUT_FILE"\
	 | "$PROGRAM_NAME" -s "$NUM_OF_SAMPLES" -f structs_bip.out ${@:4}\
	 > dump_bip.dump

cp structs_bip.out  logs/"$UNIQ_ID"/
cp dump_bip.dump logs/"$UNIQ_ID"/



cat structs_bip.out \
	| distance/sort_terminal.py \
	> mystructs_sorted_by_terminal.out

# cat mystructs_sorted_by_terminal.out \
# 	| distance/break_wins.py "$SEQUENCES_FILE" \
# 	> mystructs_broken.out

# cat mystructs_broken.out \
# 	| distance/dist_bipartite.py \
# 	> matrix.out

cat mystructs_sorted_by_terminal.out \
	| distance/dist_bipartite.py \
	> matrix.out


Rscript clustering/script_generic.r matrix.out

python clustering/read_clusters.py -c clusters.out -s mystructs_sorted_by_terminal.out \
	> representatives.out

cp mystructs_sorted_by_terminal.out logs/"$UNIQ_ID"/
# cp mystructs_broken.out logs/"$UNIQ_ID"/
# cp matrix logs/"$UNIQ_ID"/
cp representatives.out logs/"$UNIQ_ID"/

rm filename.eps
rm *.out
rm *.dump

