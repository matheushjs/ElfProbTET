require(truncnorm)

source("myoptim.r")

# @param useC Tells us to also estimate parameter C, which is the amount to subtract from the samples.
tnorm.infer = function(samples, useC=FALSE){
	estimatedC    = min(samples) * 0.995

	# The likelihood function
	likelihood = function(ourSamples, params, C=0){
		ourSamples = ourSamples - C

		allLogs = log(dtruncnorm(ourSamples, a=0, mean=params[1], sd=params[2]))

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

	if(useC == FALSE){
		sampleMean = mean(samples)
		sampleSd   = sd(samples)
	} else {
		ourSamples = samples - estimatedC
		sampleMean = mean(ourSamples)
		sampleSd   = sd(ourSamples)
	}

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
		if(useC == FALSE){
			result = myoptim(params, function(p) likelihood(samples, p), lower=lower, upper=upper, method="L-BFGS-B")
			#result = myoptim(result$par, function(p) likelihood(samples, p), lower=lower, upper=upper, method="L-BFGS-B")
		} else {
			result = myoptim(params, function(p) likelihood(samples, p, p[length(p)]), lower=lower, upper=upper, method="L-BFGS-B")
			#result = myoptim(result$par, function(p) likelihood(samples, p, p[length(p)]), lower=lower, upper=upper, method="L-BFGS-B")
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

	# Now perform cross validation using the best parameters as initial conditions
	bestResult = retval[nrow(retval),];
	if(useC == FALSE){
		params = bestResult[1,1:2];
		cross  = cross.validation(params, samples, likelihood, C=0, lower=lower, upper=upper, method="L-BFGS-B");
	} else {
		params = bestResult[1,1:3];
		cross  = cross.validation(params, samples, likelihood, C=p[length(p)], lower=lower, upper=upper, method="L-BFGS-B");
	}

	return(list(results=retval, cross=cross));
}

tnorm.lines = function(samples, params, useC=FALSE, ...){
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
	y = dtruncnorm(x, a=0, mean=params[1], sd=params[2])

	if(useC)
		x = x + params[length(params)]

	lines(x, y, ...)
}
