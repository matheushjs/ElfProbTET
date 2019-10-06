
source("./kwcwg.r")
source("./weibull.r")
source("./gamma.r")
source("./norm.r")
source("./gengamma.r")
source("./expweibull.r")

# Returns the sample sets in each file
# Each file has 1000 samples for each experiment made, so each sample set has 1000 entries
get_sample_sets = function(file){
	data = read.csv(file, header=F)
	data = as.data.frame(matrix(t(data), nrow=1000))

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
samples.hist = function(samples, breaks=20, main="", xmin=NULL){
	delta = diff(quantile(samples, c(0.05, 0.95)))
	limits = as.numeric(quantile(samples, c(0.05, 0.95))) + c(-0.95*delta, 1.05*delta)

	if(is.numeric(xmin))
		limits[1] = xmin

	# Plot the resulting pdf
	histData = hist(samples, breaks=breaks, prob=T, col="peachpuff", xlab="Execution Time (s)", main=main, xlim=limits)
}

# Process all files and save all plots
# If zeroPositioning is TRUE, we subtract the lowest execution time from the sample, making the empirical distribution begin at zero.
generate.plots = function(fullDataset, zeroPositioning=FALSE, useHeuristic=FALSE, useC=FALSE){
	for(i in 1:length(fullDataset)){
		dataset = fullDataset[[i]]

		for(j in 1:length(dataset$psizes)){
			psize = dataset$psizes[j]
			samples = dataset$samples[,j]
			sampleSize = length(samples)
			df = data.frame()
			histMinX = NULL

			if(zeroPositioning){
				samples = samples - min(samples)
				histMinX = 0
			}

			title = paste(capitalize(dataset$algorithm), capitalize(dataset$machine), psize, sep="-")

			samples.hist(samples, main=title, xmin=histMinX)

			elapsed = system.time({ retval = gamma.infer(samples, useHeuristic, useC) })["elapsed"]
			retval = retval[nrow(retval),]
			params = as.numeric(retval[1:(length(retval)-2)])
			nParams = length(params)
			errors = retval["convergence"] != 0
			errorRatio = sum(errors) / length(errors)
			minus2l = -2*retval["value"]
			gamma.lines(samples, params, useC, lty=1, col=1, lwd=3)
			df = rbind(df, c(title, "Gamma", paste.vector(params), 
							 paste(minus2l), paste(minus2l + 2*nParams), paste(minus2l + 2*nParams*sampleSize/(sampleSize - nParams - 1)),
							 paste(minus2l + 2*nParams*log(log(sampleSize))), paste(minus2l + nParams*log(sampleSize)), paste(elapsed), paste(errorRatio)),
					   stringsAsFactors=FALSE)
			colnames(df) = c("title", "model", "estimates", "-2l", "AIC", "CAIC", "BIC", "HQIC", "elapsed.time", "optErrorRatio")

			elapsed = system.time({ retval = weibull.infer(samples, useHeuristic, useC) })["elapsed"]
			retval = retval[nrow(retval),]
			params = as.numeric(retval[1:(length(retval)-2)])
			nParams = length(params)
			minus2l = -2*retval["value"]
			weibull.lines(samples, params, useC, lty=2, col=2, lwd=3)
			df = rbind(df, c(title, "Weibull", paste.vector(params),
							 paste(minus2l), paste(minus2l + 2*nParams), paste(minus2l + 2*nParams*sampleSize/(sampleSize - nParams - 1)),
							 paste(minus2l + 2*nParams*log(log(sampleSize))), paste(minus2l + nParams*log(sampleSize)), paste(elapsed), paste(errorRatio)),
					   stringsAsFactors=FALSE)

			elapsed = system.time({ retval = norm.infer(samples, useHeuristic, useC) })["elapsed"]
			retval = retval[nrow(retval),]
			params = as.numeric(retval[1:(length(retval)-2)])
			nParams = length(params)
			minus2l = -2*retval["value"]
			norm.lines(samples, params, useC, lty=3, col=3, lwd=3)
			df = rbind(df, c(title, "T.Normal", paste.vector(params),
							 paste(minus2l), paste(minus2l + 2*nParams), paste(minus2l + 2*nParams*sampleSize/(sampleSize - nParams - 1)),
							 paste(minus2l + 2*nParams*log(log(sampleSize))), paste(minus2l + nParams*log(sampleSize)), paste(elapsed), paste(errorRatio)),
					   stringsAsFactors=FALSE)

			elapsed = system.time({ retval = kwcwg.infer(samples, useHeuristic, useC) })["elapsed"]
			retval = retval[nrow(retval),]
			params = as.numeric(retval[1:(length(retval)-2)])
			nParams = length(params)
			minus2l = -2*retval["value"]
			kwcwg.lines(samples, params, useC, lty=4, col=4, lwd=3)
			df = rbind(df, c(title, "KW-CWG", paste.vector(params),
							 paste(minus2l), paste(minus2l + 2*nParams), paste(minus2l + 2*nParams*sampleSize/(sampleSize - nParams - 1)),
							 paste(minus2l + 2*nParams*log(log(sampleSize))), paste(minus2l + nParams*log(sampleSize)), paste(elapsed), paste(errorRatio)),
					   stringsAsFactors=FALSE)

			elapsed = system.time({ retval = gengamma.infer(samples, useHeuristic, useC) })["elapsed"]
			retval = retval[nrow(retval),]
			params = as.numeric(retval[1:(length(retval)-2)])
			nParams = length(params)
			minus2l = -2*retval["value"]
			gengamma.lines(samples, params, useC, lty=5, col=5, lwd=3)
			df = rbind(df, c(title, "G.Gamma", paste.vector(params),
							 paste(minus2l), paste(minus2l + 2*nParams), paste(minus2l + 2*nParams*sampleSize/(sampleSize - nParams - 1)),
							 paste(minus2l + 2*nParams*log(log(sampleSize))), paste(minus2l + nParams*log(sampleSize)), paste(elapsed), paste(errorRatio)),
					   stringsAsFactors=FALSE)

			elapsed = system.time({ retval = expweibull.infer(samples, useHeuristic, useC) })["elapsed"]
			retval = retval[nrow(retval),]
			params = as.numeric(retval[1:(length(retval)-2)])
			nParams = length(params)
			minus2l = -2*retval["value"]
			expweibull.lines(samples, params, useC, lty=6, col=6, lwd=3)
			df = rbind(df, c(title, "E.Weibull", paste.vector(params),
							 paste(minus2l), paste(minus2l + 2*nParams), paste(minus2l + 2*nParams*sampleSize/(sampleSize - nParams - 1)),
							 paste(minus2l + 2*nParams*log(log(sampleSize))), paste(minus2l + nParams*log(sampleSize)), paste(elapsed), paste(errorRatio)),
					   stringsAsFactors=FALSE)

			print(df, width=150)
			write.csv(df, file=paste(dataset$fileroot, "-", psize, ".csv", sep=""))

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
