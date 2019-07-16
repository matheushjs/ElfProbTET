./prog;
./prog;
./prog;

for i in $(seq 1000); do
	TIMEFORMAT="%E";
	time { ./prog; };
done 2> output.txt;
