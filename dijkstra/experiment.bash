function mtime {
	perf stat { ./prog &>/dev/null; } 2>&1 | grep "task-clock" | sed -e "s/^ \+//g" | cut -f1 -d" ";
};

./prog;
./prog;
./prog;

for i in $(seq 1000); do
	TIMEFORMAT="%E";
	mtime ./prog;
done 2> output.txt;
