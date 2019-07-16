./prog;
./prog;
./prog;

for i in $(seq 10000); do
	TIMEFORMAT="%E";
	time { ./prog 1>&2; };
done 2> output.txt;
