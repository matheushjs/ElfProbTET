require(GenSA)

# @param useHeuristic Tells us to use genetic algorithm as optimization function.
# @param useC Tells us to also estimate parameter C, which is the amount to subtract from the samples.
norm.infer = function(samples, useHeuristic=FALSE){
	isZero = which(samples == 0)
	samples[isZero] = min(samples[-isZero])

	# This will be shared between likelihood.base and likelihood.c
	ourSamples = samples

	# The likelihood function
	likelihood = function(params){
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
	delta = diff(quantile(samples, c(0.05, 0.95)))
	minVal = min(samples) - 0.95*delta;
	maxVal = max(samples) + 1.05*delta;

	if(minVal <= 0)
		minVal = 1e-100

	x = seq(minVal, maxVal, length=1000)
	y = dnorm(x, mean=params[1], sd=params[2])
	lines(x, y, ...)
}
