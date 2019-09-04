
gamma.infer = function(samples){
	# The likelihood function
	likelihood = function(params){
		allLogs = dgamma(samples, shape=params[1], scale=params[2], log=TRUE)
		allLogs[!is.finite(allLogs)] = log(1e-300)

		theSum = -sum(allLogs)
		
		return(theSum)
	}

	retval = NULL

	for(shape in c(0.5, 2, 5, 10, 20, 100, 200, 500, 1000, 5000, 10000, 20000)){
		# The mean of a gamma is shape*scale
		# So we estimate the scale as being scale = mean(samples) / shape
		estimatedScale = mean(samples) / shape

		for(scale in c(estimatedScale * 0.666, estimatedScale, estimatedScale*1.5)){
			params = c(shape, scale)
			# print(params)
			
			# cat("Optimizing with initial params:", params, "\n")
			result = optim(params, likelihood, method="BFGS")
			result = optim(result$par, likelihood, method="BFGS")
			params = result$par
			val = result$value
			# cat("Got params:", params, "\n")

			retval = rbind(retval, c(params, val))
		}
	}

	retval = as.data.frame(retval)
	colnames(retval) = c("shape", "scale", "value")

	# We sort it by value
	sortedIdx = sort.list(retval$value, decreasing=TRUE)
	retval = retval[sortedIdx,]

	return(retval)
}

gamma.lines = function(samples, params, ...){
	minVal = min(samples) * 0.9;
	maxVal = max(samples) * 1.1;

	x = seq(minVal, maxVal, length=200)
	y = dgamma(x, shape=params[1], scale=params[2])
	lines(x, y, ...)
}
