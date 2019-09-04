require(flexsurv)

gengamma.infer = function(samples){
	# The likelihood function
	likelihood = function(params){
		# shape > 0, scale > 0, k > 0
		allLogs = dgengamma.orig(samples, shape=params[1], scale=params[2], k=params[3], log=TRUE)

		# Set all NA and Inf to an irrelevant value
		allLogs[!is.finite(allLogs)] = log(1e-300)

		theSum = -sum(allLogs)
		
		return(theSum)
	}

	retval = NULL

	for(shape in c(0.5, 2, 5, 10))
	for(scale in c(0.5, 2, 5, 10))
	for(k in c(0.5, 2, 5, 10)){
		params = c(shape, scale, k)
		# print(params)
		
		# cat("Optimizing with initial params:", params, "\n")
		result = optim(params, likelihood, method="BFGS")
		result = optim(result$par, likelihood, method="BFGS")
		params = result$par
		val = result$value
		# cat("Got params:", params, "\n")

		retval = rbind(retval, c(params, val))
	}

	retval = as.data.frame(retval)
	colnames(retval) = c("shape", "scale", "k", "value")

	# We sort it by value
	sortedIdx = sort.list(retval$value, decreasing=TRUE)
	retval = retval[sortedIdx,]

	return(retval)
}

gengamma.lines = function(samples, params, ...){
	minVal = min(samples) * 0.9;
	maxVal = max(samples) * 1.1;

	x = seq(minVal, maxVal, length=200)
	y = dgengamma.orig(x, shape=params[1], scale=params[2], k=params[3])
	lines(x, y, ...)
}
