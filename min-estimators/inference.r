require(viridis)

source("./allModels.r");

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

#mycolors = inferno(100)[c(1, 25, 45, 66, 95)];
mycolors = c(
	"#000000FF",
	"#FF0000FF",
	3,
	"#0000FFFF",
	"#00FFFFFF",
	"#FF00FFFF",
	"#FFFF00FF",
	"#999999FF"
);

mylty = c(
	"solid",
	"62",
	"43",
	"24",
	"12"
);

getlwd = function(idx){
	return( 7 - 0.5*(idx-1) );
	#rep(3, length(idx));
}

#x = seq(0, 10, length=1000);
#idx=1; plot (x**(idx+1), lty=mylty[idx], col=mycolors[idx], lwd=7);
#idx=2; lines(x**(idx+1), lty=mylty[idx], col=mycolors[idx], lwd=7);
#idx=3; lines(x**(idx+1), lty=mylty[idx], col=mycolors[idx], lwd=7);
#idx=4; lines(x**(idx+1), lty=mylty[idx], col=mycolors[idx], lwd=7);
#idx=5; lines(x**(idx+1), lty=mylty[idx], col=mycolors[idx], lwd=7);
#quit();

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

	# I want 'breaks' bins in the interval shown, not in the whole interval
	breaks = ceiling(breaks * (max(samples) - min(samples)) / (limits[2] - limits[1]) )

	# Plot the resulting pdf
	histData = hist(samples, breaks=breaks, prob=T, col="peachpuff", xlab="Execution Time (s)", main=main, xlim=limits)
}

# Process all files and save all plots
# If zeroPositioning is TRUE, we subtract the lowest execution time from the sample, making the empirical distribution begin at zero.
generate.plots = function(fullDataset, zeroPositioning=FALSE, useC=FALSE, useMinEstimator=FALSE){
	# Research Questions:
	# 1. Even if the model has a long left tail, wont the likelihood be higher if we use the c + X model?
	# 2. If it is actually a c + X model, how to estimate c? Can we also achieve better likelihood?

	# Cautions:
	# 1. We assume the optimization finds the global maximum likelihood

	for(i in 1:100){
		data = rgamma(n=1000, shape=10, scale=10);
		print(sum(dgamma(data, shape=10, scale=10, log=T)));

		models = allModels(data);
		res = infer(models[[1]], data);
		print(res$results[nrow(res$results),"value"]);
		print(res$cross);

		factor = 1 - sd(data)/(mean(data) * log10(length(data)));
		data = data - factor*min(data);
		models = allModels(data);
		res = infer(models[[1]], data);
		print(res$results[nrow(res$results),"value"]);
		print(res$cross);
		print("======")
	}

	return();

	for(i in 1:length(fullDataset)){
		dataset = fullDataset[[i]]

		for(j in 1:length(dataset$psizes)){
			psize = dataset$psizes[j]
			samples = dataset$samples[,j]
			sampleSize = length(samples)
			df = data.frame()
			histMinX = NULL

			if(zeroPositioning){
				if(useMinEstimator){
					factor = 1 - sd(samples)/(mean(samples) * log10(sampleSize))
				} else {
					factor = 1
				}
				samples = samples - factor*min(samples)
				cat("factor: ", factor, "\n")
				histMinX = 0
			}

			graphics.off();
			dev.new(width=1*12, height=1*6)
			par(mfrow=c(1, 2));
			samples.hist(samples, xmin=histMinX)
			plotCount = 1;

			title = paste(capitalize(dataset$algorithm), capitalize(dataset$machine), psize, sep="-")
			title(title, line = -3, outer = TRUE)

			models = allModels(samples);
			for(i in 1:length(models)){
				model = models[[i]];

				elapsed = system.time({ retval = infer(model, samples, useC) })["elapsed"]
				results = retval$results;
				results = results[nrow(results),]
				params = as.numeric(results[1:(length(results)-2)])
				nParams = length(params)
				errors = results["convergence"] != 0
				errorRatio = sum(errors) / length(errors)
				minus2l = -2*results["value"]
				lines(model, samples, params, useC, lty=mylty[plotCount], col=mycolors[plotCount], lwd=getlwd(plotCount))
				df = rbind(df, c(title, model$name, paste.vector(params),
								 retval$cross,
								 paste(minus2l), paste(minus2l + 2*nParams), paste(minus2l + 2*nParams*sampleSize/(sampleSize - nParams - 1)),
								 paste(minus2l + 2*nParams*log(log(sampleSize))), paste(minus2l + nParams*log(sampleSize)), paste(elapsed), paste(errorRatio)),
						   stringsAsFactors=FALSE)
				plotCount = plotCount + 1;

				colnames(df) = c("title", "model", "estimates", "crossValid", "-2l", "AIC", "CAIC", "BIC", "HQIC", "elapsed.time", "optErrorRatio")

				if(plotCount == 6){
					legend("topright", unname(unlist(df["model"]))[1:5], lty=mylty[1:5], col=mycolors[1:5], lwd=getlwd(1:5), seg.len=5);
					samples.hist(samples, xmin=histMinX)
					plotCount = 1;
				}
			}
			legend("topright", unname(unlist(df["model"]))[6:9], lty=mylty[1:4], col=mycolors[1:4], lwd=getlwd(1:4), seg.len=5);

			print(df, width=150)
			write.csv(df, file=paste(dataset$fileroot, "-", psize, ".csv", sep=""))

			outputName = paste(dataset$fileroot, "-", psize, ".png", sep="")
			savePlot(outputName, type="png")
			graphics.off();
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
