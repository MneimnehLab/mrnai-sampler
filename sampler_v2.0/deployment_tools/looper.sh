# need ABSOLUTE PATH to weights file:
INPUT_FILE=/Users/aliahmed/research/research2/compiled_yeast_trunc_weights

# need ABSOLUTE PATH to sequences file:
SEQUENCE_FILE=/Users/aliahmed/research/research2/yeast_trunc

# specify paramters (e.g. "-m MAX_WINS", "-e", "-w MAX_WIN_SIZE")
PARAM_MAX_WINS=''
PARAM_MAX_WIN_SIZE=''
PARAM_SYMMETRIC_WINS=''

NAME=0509-yeast-1-runeq
NUM_SAMPLES_EACH_RUN=1000
NUM_RUNS=3

rm sample_counter_$NAME

for i in `seq 1 $NUM_RUNS`
do
	echo $i
	./run.sh $INPUT_FILE $NUM_SAMPLES_EACH_RUN $SEQUENCE_FILE $PARAM_MAX_WINS $PARAM_MAX_WIN_SIZE $PARAM_SYMMETRIC_WINS
done
