function mtime {
	perf stat $@ 2>&1 > /dev/null | grep "task-clock" | sed -e "s/^ \+//g" | cut -f1 -d" ";
}

./prog;
./prog;
./prog;

for i in $(seq 1000); do
	TIMEFORMAT="%E";
	mtime ./prog;
done > output.txt;
