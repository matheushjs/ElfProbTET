
source("./kwcwg.r")

# Files that we will process
files = c(
	"../sqldb_manipulation/output-100-150-300-1500.txt",
	"../mandelbrot/output_5000_10000_30000.txt",
	"../dijkstra/output_500K1M2M10M.txt"
)

# Returns the sample sets in each file
# Each file has 1000 samples for each experiment made, so each sample set has 1000 entries
get_sample_sets = function(file){
	data = read.csv(file, header=F)
	data = matrix(t(data), nrow=1000)
	return(data)
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


# Already prepares all sample sets, with all information we might need
# Make them available in the prompt, for the user
fullDataset = list()
counter = 1
for(file in files){
	samples = get_sample_sets(file)
	for(i in 1:ncol(samples)){
		dataset = list()
		dataset$samples = samples[,i]
		dataset$srcFile = file
		dataset$destFile = paste(strsplit(dataset$srcFile, "/")[[1]][2], "-", round(mean(dataset$samples), 2), ".png", sep="")

		cat("File:", file, "\n")

		fullDataset[[counter]] = dataset
		counter = counter + 1
	}
}

for(i in 1:length(fullDataset)){
	dataset = fullDataset[[i]]
	cat("Source:", dataset$srcFile, "\n")
	cat("Dest png:", dataset$destFile, "\n")
	cat("Sample count:", length(dataset$samples), "\n")
}

hints = function(){
	cat("
	dataset = fullDataset[[1]]
	hist(dataset$samples, prob=T, xlim=c(0,3))
	x = seq(0, 3, length=200)
	y = kwcwg.pdf(x, 0.5, 20, 2, 1, 1)
	lines(x, y)", "\n")
}

cat("INFO: I've lodaded the dataset in the variable 'fullDataset'", "\n")
cat("`hints()` to get help", "\n")
