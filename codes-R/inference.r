
source("./kwcwg.r")
source("./weibull.r")
source("./gamma.r")
source("./norm.r")
source("./gengamma.r")
source("./expweibull.r")

# Returns the sample sets in each file
# Each file has 1000 samples for each experiment made, so each sample set has 1000 entries
# If zeroPositioning is TRUE, we subtract the lowest execution time from the sample, making the empirical distribution begin at zero.
get_sample_sets = function(file, zeroPositioning=FALSE){
	data = read.csv(file, header=F)
	data = as.data.frame(matrix(t(data), nrow=1000))

	if(zeroPositioning){
		data = apply(data, 2, function(col){ col - min(col) })
	}

	return(data)
}

paste.vector = function(vec){
	if(length(vec) == 0)
		return("")

	retval = as.character(vec[1])
	for(item in vec[2:length(vec)]){
		retval = paste(retval, item)
	}
	return(retval)
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

		samples = get_sample_sets(paste("../experiments/", file, sep=""), zeroPositioning=TRUE)

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
			df = data.frame()

			samples.hist(samples)

			elapsed = system.time({ retval = gamma.infer(samples) })["elapsed"]
			retval = retval[nrow(retval),]
			params = as.numeric(retval[1:length(retval)-1])
			gamma.lines(samples, params, lty=2, col=2, lwd=3)
			df = rbind(df, c("Gamma", paste.vector(params), paste(-retval["value"]), paste(elapsed)), stringsAsFactors=FALSE)
			colnames(df) = c("model", "estimates", "log-likelihood", "elapsed.time")

			elapsed = system.time({ retval = weibull.infer(samples) })["elapsed"]
			retval = retval[nrow(retval),]
			params = as.numeric(retval[1:length(retval)-1])
			weibull.lines(samples, params, lty=3, col=3, lwd=3)
			df = rbind(df, c("Weibull", paste.vector(params), paste(-retval["value"]), paste(elapsed)), stringsAsFactors=FALSE)

			elapsed = system.time({ retval = norm.infer(samples) })["elapsed"]
			retval = retval[nrow(retval),]
			params = as.numeric(retval[1:length(retval)-1])
			norm.lines(samples, params, lty=4, col=4, lwd=3)
			df = rbind(df, c("Normal", paste.vector(params), paste(-retval["value"]), paste(elapsed)), stringsAsFactors=FALSE)

			elapsed = system.time({ retval = kwcwg.infer(samples) })["elapsed"]
			retval = retval[nrow(retval),]
			params = as.numeric(retval[1:length(retval)-1])
			kwcwg.lines(samples, params, lty=1, col=1, lwd=3)
			df = rbind(df, c("KW-CWG", paste.vector(params), paste(-retval["value"]), paste(elapsed)), stringsAsFactors=FALSE)

			elapsed = system.time({ retval = gengamma.infer(samples) })["elapsed"]
			retval = retval[nrow(retval),]
			params = as.numeric(retval[1:length(retval)-1])
			gengamma.lines(samples, params, lty=5, col=5, lwd=3)
			df = rbind(df, c("G.Gamma", paste.vector(params), paste(-retval["value"]), paste(elapsed)), stringsAsFactors=FALSE)

			elapsed = system.time({ retval = expweibull.infer(samples) })["elapsed"]
			retval = retval[nrow(retval),]
			params = as.numeric(retval[1:length(retval)-1])
			expweibull.lines(samples, params, lty=6, col=6, lwd=3)
			df = rbind(df, c("E.Weibull", paste.vector(params), paste(-retval["value"]), paste(elapsed)), stringsAsFactors=FALSE)

			print(df)
			print(df["model"])

			legend("topright", unname(unlist(df["model"])), lty=1:nrow(df), col=1:nrow(df), lwd=3)

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
	kwcwg.lines(dataset$samples[,1], params)
	
	generate.plots(fullDataset)", "\n")
}

cat("INFO: I've lodaded the dataset in the variable 'fullDataset'", "\n")
cat("`hints()` to get help", "\n")
