# mrnai-sampler
Sampling multiple RNA interaction structures

Running the program:

```
./run.sh data/compiled_yeast_trunc.weights 100 data/yeast_trunc.seq -e
```

Parameters to the main binary:

```
-s N : num of samples 
-i N : num of intervals
-b N : burn in steps
-r N : random start steps
-w N : max window size
-m N : max number of "good" windows to pick
-e : whether symmetric or not
```