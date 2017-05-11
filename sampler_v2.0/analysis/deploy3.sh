SCRIPT_TO_RUN=looper.sh
# RANDOM
# this should unzip into a folder "windows"
# next steps will break without this
unzip pkg.zip

cd sampler_v2.0

make;

cd ..

for i in {0..3}
do
	# unzip compiled.zip -d "instance$i";
	echo "Copy $i"
	cp -r sampler_v2.0 instance$i
	cp ~/$SCRIPT_TO_RUN instance$i/.
done;


for i in {0..3}
do
	cd instance$i
	chmod 755 $SCRIPT_TO_RUN
	echo "Running instance$i"
	# ./$SCRIPT_TO_RUN corrected_cops >> ins"$i"out  &
	./$SCRIPT_TO_RUN  &
	cd ..

	sleep 1s
done
