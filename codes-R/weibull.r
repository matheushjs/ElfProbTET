require(GenSA)

weibull.infer = function(samples, useHeuristic=FALSE){
	# This quantile is apparently a good estimator for the beta parameter
	estimatedBeta = quantile(samples, p=.632)

	# The likelihood function
	likelihood = function(params){
		allLogs = log(dweibull(samples, shape=params[1], scale=params[2]))

		# Set all NA and Inf to an irrelevant value
		allLogs[!is.finite(allLogs)] = log(1e-300)

		theSum = -sum(allLogs)
		
		return(theSum)
	}

	retval = NULL

	lower = c(1e-10, 1e-10)
	upper = c(Inf, Inf)

	for(shape in c(0.5, 2, 5, 10, 20, 30, 40))
	for(scale in c(estimatedBeta*0.666, estimatedBeta, estimatedBeta*1.5)){
		params = c(shape, scale)
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
	colnames(retval) = c("shape", "scale", "value")

	# We sort it by value
	sortedIdx = sort.list(retval$value, decreasing=TRUE)
	retval = retval[sortedIdx,]

	return(retval)
}

weibull.lines = function(samples, params, ...){
	delta = diff(quantile(samples, c(0.05, 0.95)))
	minVal = min(samples) - 0.95*delta;
	maxVal = max(samples) + 1.05*delta;

	if(minVal <= 0)
		minVal = 1e-100

	x = seq(minVal, maxVal, length=1000)
	y = dweibull(x, shape=params[1], scale=params[2])
	lines(x, y, ...)
}
