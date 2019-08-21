
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
		outFile = paste(fileRoot, ".png", sep="")

		algor = fileMetadata[1]
		machine = fileMetadata[2]
		psizes = fileMetadata[3:length(fileMetadata)]

		samples = get_sample_sets(paste("../experiments/", file, sep=""))

		entry$filename = file
		entry$algorithm = algor
		entry$machine = machine
		entry$psizes = psizes
		entry$outputFile = outFile
		entry$samples = samples

		fullDataset[[length(fullDataset)+1]] = entry
	}

	return (fullDataset)
}

# Plots a "histogram" of the dataset, and the inferred PDF function
plotme = function(dataset, params){
	# Plot the resulting pdf
	histData = hist(dataset, breaks=20, prob=T, col="peachpuff")

	x = seq(min(histData$mids)*0.9, max(histData$mids)*1.1, length=200)
	y = kwcwg.pdf(x, params[1], params[2], params[3], params[4], params[5])
	
	lines(x, y, "l")
	scan()
}

# Process all files and save all plots
main = function(){
	counter = 0

	# Each file has data from multiple experiments
	for(file in files){
		samples = get_sample_sets(file)

		# For each experiment, perform the inference and show the result
		for(i in 1:ncol(samples)){
			dataset = samples[,i]

			cat("Processing set number #", counter, "\n", sep="")
			counter = counter + 1

			# print(samples[,i])
			table = kwcwg.infer(dataset)
			print(table)

			# Get last row
			params = table[nrow(table),]
			
			# Get first 5 columns
			params = as.matrix(params[,1:5])
			cat("params:", params, "\n")

			plotme(dataset, params)
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
	params = kwcwg.infer(dataset$samples[,1]);
	params = as.matrix(params[nrow(params),])
	kwcwg.plot(dataset$samples[,1], params[1], params[2], params[3], params[4], params[5])", "\n")
}

cat("INFO: I've lodaded the dataset in the variable 'fullDataset'", "\n")
cat("`hints()` to get help", "\n")
