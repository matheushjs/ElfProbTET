allDirs="bfgs-estimate-c bfgs-subtract-min bfgs-subtract-estimated-min bfgs-no-subtract-c"

for dir in $allDirs; do
	ofile="$dir".csv
	files="$(ls $dir/*.csv)"

	arr=($files)
	firstFile=${arr[0]}
	cat $firstFile | head -1 | sed -e "s/^/\"directory\",/g" > $ofile

	for i in $files; do
		cat $i | tail +2 | sed -e "s/^/\"$dir\",/g";
	done >> $ofile

	ofile="$dir".pdf
	files="$(ls $dir/*.png)"
	magick $files $ofile
done

arr=($allDirs)
firstFile=${arr[0]}.csv
cat $firstFile | head -1 > all.csv

for dir in $allDirs; do
	cat $dir.csv | tail +2
done >> all.csv
