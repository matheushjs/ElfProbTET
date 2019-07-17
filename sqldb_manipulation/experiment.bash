python3 ./prog;
python3 ./prog;
python3 ./prog;

for i in $(seq 10000); do
	TIMEFORMAT="%E";
	time { python3 ./prog 1>&2; };
done 2> output.txt;
