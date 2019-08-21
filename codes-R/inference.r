
source("./kwcwg.r")
source("./weibull.r")
source("./gamma.r")
source("./norm.r")

# Returns the sample sets in each file
# Each file has 1000 samples for each experiment made, so each sample set has 1000 entries
get_sample_sets = function(file){
	data = read.csv(file, header=F)
	data = as.data.frame(matrix(t(data), nrow=1000))
	return(data)
}

# Already prepares all sample sets, with all information we might need
# Returns a data frame with all information we need
#    - origin file name
#    - index of first sample of this experiment
#    - problem size of this experiment
#    - machine in which the experiment was executed
#    - algorithm in question
experiment.files = function(){
	cmd = paste("ls", "../experiments")
	files = system(cmd, intern=T)

	fullDataset = list()

	for(file in files){
		entry = list()

		# remove extension
		fileRoot = strsplit(file, ".txt")[[1]]
		fileMetadata = strsplit(fileRoot, "_")[[1]]

		algor = fileMetadata[1]
		machine = fileMetadata[2]
		psizes = fileMetadata[3:length(fileMetadata)]

		samples = get_sample_sets(paste("../experiments/", file, sep=""))

		entry$filename = file
		entry$algorithm = algor
		entry$machine = machine
		entry$psizes = psizes
		entry$fileroot = fileRoot
		entry$samples = samples

		fullDataset[[length(fullDataset)+1]] = entry
	}

	return (fullDataset)
}

# Plots a "histogram" of the dataset
samples.hist = function(samples, breaks=20){
	minVal = min(samples) * 0.9;
	maxVal = max(samples) * 1.1;

	# Plot the resulting pdf
	histData = hist(samples, breaks=breaks, prob=T, col="peachpuff", xlab="Execution Time (s)")
}

# Process all files and save all plots
generate.plots = function(fullDataset){
	for(i in 1:length(fullDataset)){
		dataset = fullDataset[[i]]

		for(j in 1:length(dataset$psizes)){
			psize = dataset$psizes[j]
			samples = dataset$samples[,j]

			samples.hist(samples)

			params = kwcwg.infer(samples)
			params = as.matrix(params[nrow(params),])
			kwcwg.lines(samples, params, lty=1, col=1, lwd=3)

			params = gamma.infer(samples)
			params = as.matrix(params[nrow(params),])
			gamma.lines(samples, params, lty=2, col=2, lwd=3)

			params = weibull.infer(samples)
			params = as.matrix(params[nrow(params),])
			weibull.lines(samples, params, lty=3, col=3, lwd=3)

			outputName = paste(dataset$fileroot, "-", psize, ".png", sep="")
			savePlot(outputName, type="png")
		}
	}
}

# main()

fullDataset = experiment.files()

cat("Variable: fullDataset\n")
cat("# entries: ", length(fullDataset), "\n")
cat("Each entry has:", paste(names(fullDataset[[1]])), sep="\n\t$")

hints = function(){
	cat("
	dataset = fullDataset[[1]]
	samples.hist(dataset$samples[,1])
	params = kwcwg.infer(dataset$samples[,1]);
	params = as.matrix(params[nrow(params),])
	kwcwg.lines(dataset$samples[,1], params[1], params[2], params[3], params[4], params[5])
	
	generate.plots(fullDataset)", "\n")
}

cat("INFO: I've lodaded the dataset in the variable 'fullDataset'", "\n")
cat("`hints()` to get help", "\n")
