
norm.infer = function(samples){
	# The likelihood function
	likelihood = function(params){
		allLogs = log(dnorm(samples, mean=params[1], sd=params[2]))

		# Set all NA and Inf to 0
		allLogs[!is.finite(allLogs)] = 0

		theSum = -sum(allLogs)
		
		return(theSum)
	}

	sampleMean = mean(samples)
	sampleSd   = sd(samples)

	retval = NULL

	for(sd in c(sampleSd*0.666, sampleSd, sampleSd*1.5)){
		for(mean in c(sampleMean* 0.666, sampleMean, sampleMean*1.5)){
			params = c(mean, sd)
			# print(params)
			
			cat("Optimizing with initial params:", params, "\n")
			result = optim(params, likelihood)
			result = optim(result$par, likelihood)
			params = result$par
			val = result$value
			cat("Got params:", params, "\n")

			retval = rbind(retval, c(params, val))
		}
	}

	retval = as.data.frame(retval)
	colnames(retval) = c("mean", "stdDev", "value")

	# We sort it by value
	sortedIdx = sort.list(retval$value, decreasing=TRUE)
	retval = retval[sortedIdx,]

	return(retval)
}

norm.lines = function(samples, params, ...){
	minVal = min(samples) * 0.9;
	maxVal = max(samples) * 1.1;

	x = seq(minVal, maxVal, length=200)
	y = dnorm(x, mean=params[1], sd=params[2], ...)
	lines(x, y)
}
