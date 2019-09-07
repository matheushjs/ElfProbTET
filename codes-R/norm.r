require(GenSA)

norm.infer = function(samples, useHeuristic=FALSE){
	# The likelihood function
	likelihood = function(params){
		allLogs = log(dnorm(samples, mean=params[1], sd=params[2]))

		# Set all NA and Inf to an irrelevant value
		allLogs[!is.finite(allLogs)] = log(1e-300)

		theSum = -sum(allLogs)
		
		return(theSum)
	}

	sampleMean = mean(samples)
	sampleSd   = sd(samples)

	retval = NULL

	lower = c(-Inf, 1e-10)
	upper = c(Inf, Inf)

	for(sd in c(sampleSd*0.666, sampleSd, sampleSd*1.5))
	for(mean in c(sampleMean* 0.666, sampleMean, sampleMean*1.5)){
		params = c(mean, sd)
		# print(params)
		
		# cat("Optimizing with initial params:", params, "\n")
		if(useHeuristic == FALSE){
			result = optim(params, likelihood, method="BFGS")
			result = optim(result$par, likelihood, method="BFGS")
		} else {
			result = GenSA(params, likelihood, lower=lower, upper=upper)
		}

		params = result$par
		val = result$value
		# cat("Got params:", params, "\n")

		retval = rbind(retval, c(params, val))
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
	y = dnorm(x, mean=params[1], sd=params[2])
	lines(x, y, ...)
}
