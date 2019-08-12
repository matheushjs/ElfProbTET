function mtime {
	perf stat $@ 2>&1 > /dev/null | grep "seconds time elapsed" | sed -e "s/^ \+//g" | cut -f1 -d" ";
}

python3 ./prog.py;
python3 ./prog.py;
python3 ./prog.py;

for psize in 100 150 300 1500; do
	for i in $(seq 1000); do
		mtime python3 ./prog.py $psize;
	done;
done > output.txt;
