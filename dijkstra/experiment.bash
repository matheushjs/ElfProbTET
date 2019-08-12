function mtime {
	perf stat $@ 2>&1 > /dev/null | grep "seconds time elapsed" | sed -e "s/^ \+//g" | cut -f1 -d" ";
}

./prog;
./prog;
./prog;

for psize in 500000 1000000 2000000 10000000; do
	for i in $(seq 1000); do
		mtime ./prog $psize;
	done;
done > output.txt;
