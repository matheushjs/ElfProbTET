require(rmutil)
require(GenSA)

# @param useHeuristic Tells us to use genetic algorithm as optimization function.
# @param useC Tells us to also estimate parameter C, which is the amount to subtract from the samples.
expweibull.infer = function(samples, useHeuristic=FALSE, useC=FALSE){
	estimatedC    = min(samples) * 0.995

	# dgweibull does not accept a sample of value 0, so we fix this here
	isZero = which(samples == 0)
	samples[isZero] = min(samples[-isZero])

	# The likelihood function
	likelihood = function(params, C=0){
		ourSamples = samples - C;
		# shape s > 0, scale m > 0, family f > 0
		allLogs = dgweibull(ourSamples, s=params[1], m=params[2], f=params[3], log=TRUE)

		problems = which(!is.finite(allLogs))
		if(length(problems) > 0 && length(problems) <= 5){
			allLogs[problems] = min(allLogs[-problems])
		} else {
			allLogs[problems] = log(1e-300)
		}

		if(length(problems) > 0 && length(problems) < 5){
			warning(paste("expweibul: Low amount (<5) of warnings at points:", ourSamples[problems]), call.=FALSE)
		}

		theSum = -sum(allLogs) # Invert so that minimization yields a maximum
		
		return(theSum)
	}

	retval = NULL

	if(useC == FALSE){
		lower = c(1e-10, 1e-10, 1e-10)
		upper = c(Inf, Inf, Inf)
	} else {
		lower = c(1e-10, 1e-10, 1e-10, 1e-10)
		upper = c(Inf, Inf, Inf, min(samples) - 1e-10)
	}

	for(shape in c(0.5, 2, 5, 10))
	for(scale in c(0.5, 2, 5, 10))
	for(family in c(0.5, 2, 5, 10)){
		params = c(shape, scale, family)
		if(useC)
			params = c(params, estimatedC)
		# print(params)
		
		# cat("Optimizing with initial params:", params, "\n")
		if(useHeuristic == FALSE){
			if(useC == FALSE){
				result = optim(params, function(p) likelihood(p), lower=lower, upper=upper, method="BFGS")
				result = optim(result$par, function(p) likelihood(p), lower=lower, upper=upper, method="BFGS")
			} else {
				result = optim(params, function(p) likelihood(p, p[length(p)]), lower=lower, upper=upper, method="BFGS")
				result = optim(result$par, function(p) likelihood(p, p[length(p)]), lower=lower, upper=upper, method="BFGS")
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
		colnames(retval) = c("shape", "scale", "family", "value", "convergence")
	} else {
		colnames(retval) = c("shape", "scale", "family", "c", "value", "convergence")
	}

	# We sort it by value
	sortedIdx = sort.list(retval$value, decreasing=FALSE)
	retval = retval[sortedIdx,]

	return(retval)
}

expweibull.lines = function(samples, params, useC=FALSE, ...){
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
	x[x == 0] = 1e-100
	y = dgweibull(x, s=params[1], m=params[2], f=params[3])

	if(useC)
		x = x + params[length(params)]

	lines(x, y, ...)
}
