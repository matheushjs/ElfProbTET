require(GenSA)

gamma.infer = function(samples, useHeuristic=FALSE){
	# The likelihood function
	likelihood = function(params){
		allLogs = dgamma(samples, shape=params[1], scale=params[2], log=TRUE)
		allLogs[!is.finite(allLogs)] = log(1e-300)

		theSum = -sum(allLogs)
		
		return(theSum)
	}

	retval = NULL

	lower = c(1e-10, 1e-10)
	upper = c(2000, 2000)

	if(useHeuristic == FALSE){
		shapeList = c(0.5, 2, 5, 10, 20, 100, 200, 500, 1000, 5000, 10000, 20000)
	} else {
		shapeList = c(5)
	}

	for(shape in shapeList){
		# The mean of a gamma is shape*scale
		# So we estimate the scale as being scale = mean(samples) / shape
		estimatedScale = mean(samples) / shape

		if(useHeuristic == FALSE){
			scaleList = c(estimatedScale * 0.666, estimatedScale, estimatedScale*1.5)
		} else {
			scaleList = c(estimatedScale)
		}

		for(scale in scaleList){
			params = c(shape, scale)
			# print(params)
			
			#cat("Optimizing with initial params:", params, "\n")
			if(useHeuristic == FALSE){
				result = optim(params, likelihood, method="BFGS")
				result = optim(result$par, likelihood, method="BFGS")
			} else {
				result = GenSA(params, likelihood, lower=lower, upper=upper)
			}

			params = result$par
			val = result$value
			#cat("Got params:", params, "\n")

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
	delta = diff(quantile(samples, c(0.05, 0.95)))
	minVal = min(samples) - 0.95*delta;
	maxVal = max(samples) + 1.05*delta;

	if(minVal <= 0)
		minVal = 1e-100

	x = seq(minVal, maxVal, length=1000)
	y = dgamma(x, shape=params[1], scale=params[2])
	lines(x, y, ...)
}
