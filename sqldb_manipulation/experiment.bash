function mtime {
	perf stat { ./prog &>/dev/null; } 2>&1 | grep "task-clock" | sed -e "s/^ \+//g" | cut -f1 -d" "
}

python3 ./prog.py;
python3 ./prog.py;
python3 ./prog.py;

for i in $(seq 1000); do
	TIMEFORMAT="%E";
	mtime python3 ./prog.py;
done 2> output.txt;
