require(rmutil)
require(GenSA)

expweibull.infer = function(samples, useHeuristic=FALSE){
	# The likelihood function
	likelihood = function(params){
		# dgweibull does not accept a sample of value 0, so we fix this here
		samples[samples == 0] = 1e-100

		# shape s > 0, scale m > 0, family f > 0
		allLogs = dgweibull(samples, s=params[1], m=params[2], f=params[3], log=TRUE)

		# Set all NA and Inf to 0
		allLogs[!is.finite(allLogs)] = log(1e-300)

		theSum = -sum(allLogs)
		
		return(theSum)
	}

	retval = NULL

	lower = c(1e-10, 1e-10, 1e-10)
	upper = c(Inf, Inf, Inf)

	for(shape in c(0.5, 2, 5, 10))
	for(scale in c(0.5, 2, 5, 10))
	for(family in c(0.5, 2, 5, 10)){
		params = c(shape, scale, family)
		# print(params)
		
		# cat("Optimizing with initial params:", params, "\n")
		if(useHeuristic == FALSE){
			result = optim(params, likelihood, lower=lower, upper=upper, method="BFGS")
			result = optim(result$par, likelihood, lower=lower, upper=upper, method="BFGS")
		} else {
			result = GenSA(params, likelihood, lower=lower, upper=upper)
		}
		
		params = result$par
		val = result$value
		# cat("Got params:", params, "\n")

		retval = rbind(retval, c(params, val))
	}

	retval = as.data.frame(retval)
	colnames(retval) = c("shape", "scale", "family", "value")

	# We sort it by value
	sortedIdx = sort.list(retval$value, decreasing=TRUE)
	retval = retval[sortedIdx,]

	return(retval)
}

expweibull.lines = function(samples, params, ...){
	delta = diff(quantile(samples, c(0.05, 0.95)))
	minVal = min(samples) - 0.95*delta;
	maxVal = max(samples) + 1.05*delta;

	x = seq(minVal, maxVal, length=1000)
	x[x == 0] = 1e-100
	y = dgweibull(x, s=params[1], m=params[2], f=params[3])
	lines(x, y, ...)
}
