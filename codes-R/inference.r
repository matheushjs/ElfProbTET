
source("./kwcwg.r")

# Files that we will process
files = c(
	"../sqldb_manipulation/output-100-150-300-1500.txt",
	"../mandelbrot/output_5000_10000_30000.txt",
	"../dijkstra/output_500K1M2M10B.txt"
)

# Plots a "histogram" of the dataset, and the inferred PDF function
plotme = function(dataset, params){
	# Plot the resulting pdf
	histData = hist(dataset, breaks=25, plot=F)
	plot(histData$mids, histData$density)
	#plot(density(dataset))

	x = seq(min(histData$mids)*0.9, max(histData$mids)*1.1, by=0.005)
	lines(x, kwcwg.pdf(x, params[1], params[2], params[3], params[4], params[5]), "l")
	scan()
}

# Returns the sample sets in each file
# Each file has 1000 samples for each experiment made, so each sample set has 1000 entries
get_sample_sets = function(file){
	data = read.csv(file, header=F)
	data = matrix(t(data), nrow=1000)
	return(data)
}

# Process all files and save all plots
main = function(){
	counter = 0

	# Each file has data from multiple experiments
	for(file in files){
		samples = get_sample_sets(file)

		# For each experiment, perform the inference and show the result
		for(i in 1:ncol(samples)){
			cat("Processing set number #", counter, "\n", sep="")
			counter = counter + 1

			# print(samples[,i])
			params = kwcwg.infer(samples[,i])
			plotme(samples[,i], params[nrow(params),1:5])

			break;
		}
		break;
	}
}

main()
