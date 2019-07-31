function mtime {
	perf stat $@ 2>&1 > /dev/null | grep "task-clock" | sed -e "s/^ \+//g" | cut -f1 -d" ";
}

python3 ./prog.py;
python3 ./prog.py;
python3 ./prog.py;

for i in $(seq 1000); do
	TIMEFORMAT="%E";
	mtime python3 ./prog.py;
done > output.txt;
