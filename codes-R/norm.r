require(GenSA)

# @param useHeuristic Tells us to use genetic algorithm as optimization function.
# @param useC Tells us to also estimate parameter C, which is the amount to subtract from the samples.
norm.infer = function(samples, useHeuristic=FALSE, useC=FALSE){
	estimatedC    = min(samples) * 0.995

	isZero = which(samples == 0)
	samples[isZero] = min(samples[-isZero])

	# The likelihood function
	likelihood = function(params, C=0){
		ourSamples = samples - C

		allLogs = dnorm(ourSamples, mean=params[1], sd=params[2], log=TRUE)

		problems = which(!is.finite(allLogs))
		if(length(problems) > 0 && length(problems) <= 5){
			allLogs[problems] = min(allLogs[-problems])
		} else {
			allLogs[problems] = log(1e-300)
		}

		if(length(problems) > 0 && length(problems) < 5){
			warning(paste("norm: Low amount (<5) of warnings at points:", ourSamples[problems]), call.=FALSE)
		}

		theSum = -sum(allLogs) # Invert so that minimization yields a maximum
		
		return(theSum)
	}

	sampleMean = mean(samples)
	sampleSd   = sd(samples)

	retval = NULL

	if(useC == FALSE){
		lower = c(-Inf, 1e-10)
		upper = c(Inf, Inf)
	} else {
		lower = c(-Inf, 1e-10, 1e-10)
		upper = c(Inf, Inf, min(samples) - 1e-10)
	}

	for(sd in c(sampleSd*0.666, sampleSd, sampleSd*1.5))
	for(mean in c(sampleMean* 0.666, sampleMean, sampleMean*1.5)){
		params = c(mean, sd)
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
		colnames(retval) = c("mean", "stdDev", "value", "convergence")
	} else {
		colnames(retval) = c("mean", "stdDev", "c", "value", "convergence")
	}

	# We sort it by value
	sortedIdx = sort.list(retval$value, decreasing=FALSE)
	retval = retval[sortedIdx,]

	return(retval)
}

norm.lines = function(samples, params, useC=FALSE, ...){
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
	y = dnorm(x, mean=params[1], sd=params[2])

	if(useC)
		x = x + params[length(params)]

	lines(x, y, ...)
}
