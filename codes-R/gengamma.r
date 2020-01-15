require(ggamma)
require(GenSA)

source("myoptim.r")

# @param useHeuristic Tells us to use genetic algorithm as optimization function.
# @param useC Tells us to also estimate parameter C, which is the amount to subtract from the samples.
gengamma.infer = function(samples, useHeuristic=FALSE, useC=FALSE){
	estimatedC    = min(samples) * 0.995

	# The likelihood function
	likelihood = function(params, C=0){
		ourSamples = samples - C

		isZero = which(ourSamples == 0)
		if(sum(isZero) > 0)
			ourSamples[isZero] = min(ourSamples[-isZero])

		# shape > 0, scale > 0, k > 0
		allLogs = dggamma(ourSamples, a=params[1], b=params[2], k=params[3], log=TRUE)

		problems = which(!is.finite(allLogs))
		if(length(problems) > 0 && length(problems) <= 5){
			allLogs[problems] = min(allLogs[-problems])
		} else {
			allLogs[problems] = log(1e-300)
		}

		if(length(problems) > 0 && length(problems) < 5){
			warning(paste("gengamma: Low amount (<5) of warnings at points:", ourSamples[problems]), call.=FALSE)
		}

		theSum = -sum(allLogs) # Invert so that minimization yields a maximum
		
		return(theSum)
	}

	retval = NULL

	if(useC == FALSE){
		lower = c(1e-7, 1e-7, 1e-7)
		upper = c(Inf, Inf, Inf)
	} else {
		lower = c(1e-7, 1e-7, 1e-7, 1e-7)
		upper = c(Inf, Inf, Inf, min(samples) - 1e-7)
	}

	for(shape in c(0.5, 2, 5, 10))
	for(scale in c(0.5, 2, 5, 10))
	for(k in c(0.5, 2, 5, 10)){
		params = c(shape, scale, k)
		if(useC)
			params = c(params, estimatedC)
		# print(params)
		
		# cat("Optimizing with initial params:", params, "\n")
		if(useHeuristic == FALSE){
			if(useC == FALSE){
				result = myoptim(params, function(p) likelihood(p), lower=lower, upper=upper, method="L-BFGS-B")
				result = myoptim(result$par, function(p) likelihood(p), lower=lower, upper=upper, method="L-BFGS-B")
			} else {
				result = myoptim(params, function(p) likelihood(p, p[length(p)]), lower=lower, upper=upper, method="L-BFGS-B")
				result = myoptim(result$par, function(p) likelihood(p, p[length(p)]), lower=lower, upper=upper, method="L-BFGS-B")
			}
		} else {
			result = GenSA(params, likelihood, lower=lower, upper=upper)
		}

		params = result$par
		val = -result$value # Undo signal invertion in the likelihood function
		convergence = result$convergence
		# cat("Got params:", params, "\n")

		retval = rbind(retval, c(params, val, convergence))
	}

	retval = as.data.frame(retval)
	if(useC == FALSE){
		colnames(retval) = c("a", "b", "k", "value", "convergence")
	} else {
		colnames(retval) = c("a", "b", "k", "c", "value", "convergence")
	}

	# We sort it by value
	sortedIdx = sort.list(retval$value, decreasing=FALSE)
	retval = retval[sortedIdx,]

	return(retval)
}

gengamma.lines = function(samples, params, useC=FALSE, ...){
	delta = diff(quantile(samples, c(0.05, 0.95)))
	minVal = min(samples) - 0.95*delta;
	maxVal = max(samples) + 1.05*delta;

	if(minVal <= 0)
		minVal = 1e-100

	if(useC){
		minVal = minVal - params[length(params)]
		maxVal = maxVal - params[length(params)]
	}

	x = seq(minVal, maxVal, length=1000)
	y = dggamma(x, a=params[1], b=params[2], k=params[3])
	
	if(useC)
		x = x + params[length(params)]

	lines(x, y, ...)
}
