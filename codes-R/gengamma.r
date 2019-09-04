require(flexsurv)

gengamma.infer = function(samples){
	# The likelihood function
	likelihood = function(params){
		# mu > 0, sigma > 0, Q real
		allLogs = dgengamma(samples, mu=params[1], sigma=params[2], Q=params[3], log=TRUE)

		# Set all NA and Inf to 0
		allLogs[!is.finite(allLogs)] = 0

		theSum = -sum(allLogs)
		
		return(theSum)
	}

	lower = c(0, 0, -Inf)
	upper = Inf

	retval = NULL

	for(mu in c(0.5, 2, 5, 10))
	for(sigma in c(0.5, 2, 5, 10))
	for(Q in c(-10, -5, -2, -0.5, 0.5, 2, 5, 10)){
		params = c(mu, sigma, Q)
		# print(params)
		
		cat("Optimizing with initial params:", params, "\n")
		result = optim(params, likelihood, lower=lower, upper=upper)
		result = optim(result$par, likelihood, lower=lower, upper=upper)
		params = result$par
		val = result$value
		cat("Got params:", params, "\n")

		retval = rbind(retval, c(params, val))
	}

	retval = as.data.frame(retval)
	colnames(retval) = c("mu", "sigma", "Q", "value")

	# We sort it by value
	sortedIdx = sort.list(retval$value, decreasing=TRUE)
	retval = retval[sortedIdx,]

	return(retval)
}

gengamma.lines = function(samples, params, ...){
	minVal = min(samples) * 0.9;
	maxVal = max(samples) * 1.1;

	x = seq(minVal, maxVal, length=200)
	y = dgengamma(x, mu=params[1], sigma=params[2], Q=params[3])
	lines(x, y, ...)
}
