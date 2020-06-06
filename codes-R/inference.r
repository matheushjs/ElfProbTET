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
	"#000000CC",
	"#FF0000CC",
	"#00CC00CC",
	"#0000FFCC",
	"#00FFFFCC",
	"#FF00FFCC",
	"#FFFF00CC",
	"#999999CC"
);

mypch = c(19, 3, 17, 4, 15, 2);
mycex = seq(1.5, 1, length=6);

#mycolors = c(
#	"#ff0000", # Red
#	"#0000ff", # Blue
#	"#00cc00", # Green
#	"#ff00ff", # Magenta
#	"#00e0e0", # Cyan
#	"#7f00ff", # Purple
#	"#dddd00", # Yellow
#	"#007ca8", # Slate-blue (similar)
#	"#ff7000", # Orange
#	"#b2a200"  # Darker yellow
#);

#mycolors = colorspace::qualitative_hcl(5, palette="Dark 2");
#mycolors = viridis(5);
#mycolors = RColorBrewer::brewer.pal(5, "Set1");
#mycolors = c(
#	rgb(0.2980392156862745, 0.4470588235294118, 0.6901960784313725),
#	rgb(0.8666666666666667, 0.5176470588235295, 0.3215686274509804),
#	rgb(0.3333333333333333, 0.6588235294117647, 0.40784313725490196),
#	rgb(0.7686274509803922, 0.3058823529411765, 0.3215686274509804),
#	rgb(0.5058823529411764, 0.4470588235294118, 0.7019607843137254),
#	rgb(0.5764705882352941, 0.47058823529411764, 0.3764705882352941),
#	rgb(0.8549019607843137, 0.5450980392156862, 0.7647058823529411),
#	rgb(0.5490196078431373, 0.5490196078431373, 0.5490196078431373),
#	rgb(0.8, 0.7254901960784313, 0.4549019607843137),
#	rgb(0.39215686274509803, 0.7098039215686275, 0.803921568627451)
#)

#Show colors
dev.new(width=1*12, height=1*6)
par(mfrow=c(1, 2));
x = seq(0, 10, length=10000);
plot(-1, -1, xlim=c(0, 10), ylim=c(-1, 1));
for(i in seq_along(mycolors)){
	lines(x, sin(x/2 - i/5), lwd=6, col=mycolors[i]);
}
x = seq_along(mycolors);
plot(x, rep(1, length(x)), pch=16, cex=7, col=mycolors);

#

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
	#rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "#e3e3e3")
	#histData = hist(samples, breaks=breaks, prob=T, col="peachpuff", xlab="Execution Time (s)", main=main, xlim=limits, add=T)
}

# Process all files and save all plots
# If minEstimator is TRUE, we subtract the lowest execution time from the sample, making the empirical distribution begin at zero.
# @param minEstimator can be FALSE if data should not be subtracted from some estimated populational minimum.
#                        otherwise it can be "c1", "c2", "c3", "c4" depending on the estimator you want to use.
# @param plotType can be "pdf" for histogram and pdf plotting, or "pp" for the PP-plot
generate.plots = function(fullDataset, minEstimator=FALSE, useC=FALSE, iteratedC=FALSE, plotType="pdf"){
	for(i in 1:length(fullDataset)){
		dataset = fullDataset[[i]]

		for(j in 1:length(dataset$psizes)){
			psize = dataset$psizes[j]
			samples = dataset$samples[,j]
			sampleSize = length(samples)
			df = data.frame()
			histMinX = NULL

			if(is.character(minEstimator)){
				hatMin = switch(minEstimator,
					c1 = min(samples) - abs(min(samples) * sd(samples) / mean(samples)) / log10(sampleSize),
					c2 = min(samples) - sd(samples) / log10(sampleSize),
					c3 = min(samples) - sd(samples) * sqrt(log(log(sampleSize)) / (2*sampleSize)),
					c4 = min(samples) - sd(samples) * sqrt(-log(0.05/2) / (2*sampleSize))
					)
				histMinX = 0
				samples = samples - hatMin;
				cat("hatMin: ", hatMin, "\n")
			}
			if(useC || iteratedC){
				c2 = min(samples) - sd(samples) / log10(sampleSize);
				models = allModels(samples - c2);
			} else {
				models = allModels(samples);
			}



			graphics.off();
			dev.new(width=1*12, height=1*6)
			par(mfrow=c(1, 2));

			if(plotType == "pdf"){
				samples.hist(samples, xmin=histMinX);
			} else if(plotType == "pp"){
				plot(-1, -1, xlim=c(0, sqrt(2)), ylim=c(-0.15, 0.15), xlab="", ylab="", axes=F);
				axis(2);
				box();
				title(ylab="disparity", line=2);
			} else {
				stop("Invalid plotType");
			}
			plotCount = 1;

			title = paste(capitalize(dataset$algorithm), capitalize(dataset$machine), psize, sep="-")
			title(title, line = -3, outer = TRUE)

			for(i in 1:length(models)){
				model = models[[i]];

				elapsed = system.time({ retval = infer(model, samples, useC, iteratedC) })["elapsed"]
				results = retval$results;
				results = results[nrow(results),]
				params = as.numeric(results[1:(length(results)-2)])
				nParams = length(params)
				errors = results["convergence"] != 0
				errorRatio = sum(errors) / length(errors)
				minus2l = -2*results["value"]

				if(plotType == "pdf"){
					lines(model, samples, params, useC, iteratedC, lty=mylty[plotCount], col=mycolors[plotCount], lwd=getlwd(plotCount))
				} else if(plotType == "pp"){
					ppplot(model, samples, params, useC, iteratedC, col=mycolors[plotCount], pch=mypch[plotCount], cex=mycex[plotCount], lwd=2);
				}

				df = rbind(df, c(title, model$name, paste.vector(params),
								 -2*retval$cross,
								 paste(minus2l), paste(minus2l + 2*nParams), paste(minus2l + 2*nParams*sampleSize/(sampleSize - nParams - 1)),
								 paste(minus2l + 2*nParams*log(log(sampleSize))), paste(minus2l + nParams*log(sampleSize)), paste(elapsed), paste(errorRatio)),
						   stringsAsFactors=FALSE)
				plotCount = plotCount + 1;

				colnames(df) = c("title", "model", "estimates", "crossValid", "-2l", "AIC", "CAIC", "BIC", "HQIC", "elapsed.time", "optErrorRatio")

				if(plotCount == 6){
					if(plotType == "pdf"){
						legend("topright", unname(unlist(df["model"]))[1:5], lty=mylty[1:5], col=mycolors[1:5], lwd=getlwd(1:5), seg.len=5);
						samples.hist(samples, xmin=histMinX)
					} else if(plotType == "pp"){
						legend("topleft", unname(unlist(df["model"]))[1:5], col=mycolors[1:5], pch=mypch[1:5], pt.lwd=2);
						plot(-1, -1, xlim=c(0, sqrt(2)), ylim=c(-0.15, 0.15), xlab="", ylab="", axes=F);
						axis(2);
						box();
						title(ylab="disparity", line=2);
					}

					plotCount = 1;
				}
			}

			if(plotType == "pdf"){
				legend("topright", unname(unlist(df["model"]))[6:9], lty=mylty[1:4], col=mycolors[1:4], lwd=getlwd(1:4), seg.len=5);
			} else if(plotType == "pp"){
				legend("topleft", unname(unlist(df["model"]))[6:9], col=mycolors[1:4], pch=mypch[1:4], pt.lwd=2);
			}

			outputName = paste(dataset$fileroot, "-", psize, sep="");
			if(minEstimator != FALSE){
				outputName = paste(outputName, "-", minEstimator, sep="");
			}
			if(useC != FALSE){
				outputName = paste(outputName, "-inferC", sep="");
			}
			if(iteratedC != FALSE){
				outputName = paste(outputName, "-iteratedC", sep="");
			}
			outputName = paste(outputName, "-", plotType, sep="");


			print(df, width=150)
			write.csv(df, file=paste(outputName, ".csv", sep=""))
			savePlot(paste(outputName, ".png", sep=""), type="png")
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
	models = allModels(dataset$samples[,1]);
	params = infer(models[[1]], dataset$samples[,1]);
	params = as.matrix(params$results[nrow(params$results),])
	lines(models[[1]], dataset$samples[,1], params)
	
	generate.plots(fullDataset)", "\n")
}

cat("INFO: I've lodaded the dataset in the variable 'fullDataset'", "\n")
cat("`hints()` to get help", "\n")


# All experiments

for(type in c("pdf", "pp")){
	generate.plots(fullDataset, plotType=type);
	generate.plots(fullDataset, useC=TRUE, plotType=type);
	generate.plots(fullDataset, iteratedC=TRUE, plotType=type);

	for(estim in paste("c", 1:4, sep=""))
		generate.plots(fullDataset, minEstimator=estim, plotType=type);
}
