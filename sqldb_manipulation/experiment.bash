python3 ./prog.py;
python3 ./prog.py;
python3 ./prog.py;

for i in $(seq 10000); do
	TIMEFORMAT="%E";
	time { python3 ./prog.py; };
done 2> output.txt;
