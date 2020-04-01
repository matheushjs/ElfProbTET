require(ollggamma)

source("myoptim.r")

# @param useC Tells us to also estimate parameter C, which is the amount to subtract from the samples.
ollgengamma.infer = function(samples, useC=FALSE){
	estimatedC    = min(samples) * 0.995

	# The likelihood function
	likelihood = function(ourSamples, params, C=0){
		ourSamples = ourSamples - C

		# shape > 0, scale > 0, k > 0, lambda > 0
		allLogs = dollggamma(ourSamples, a=params[1], b=params[2], k=params[3], lambda=params[4], log=TRUE)

		problems = which(!is.finite(allLogs))
		allLogs[problems] = log(1e-300); # Merely to force optim to continue optimizing
		if(length(problems) > 0 && length(problems) < 5){
			warning(paste("oll-gengamma: Low amount (<5) of warnings at points:", ourSamples[problems]), call.=FALSE)
		}

		theSum = -sum(allLogs) # Invert so that minimization yields a maximum
		
		return(theSum)
	}

	retval = NULL

	if(useC == FALSE){
		lower = c(1e-7, 1e-7, 1e-7, 1e-7)
		upper = c(Inf, Inf, Inf, Inf)
	} else {
		lower = c(1e-7, 1e-7, 1e-7, 1e-7, 1e-7)
		upper = c(Inf, Inf, Inf, Inf, min(samples) - 1e-7)
	}

	for(shape in c(0.5, 2, 5))
	for(scale in c(0.5, 2, 5))
	for(k in c(0.5, 5))
	for(lambda in c(0.25, 1.5)){
		params = c(shape, scale, k, lambda)
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
		colnames(retval) = c("a", "b", "k", "lambda", "value", "convergence")
	} else {
		colnames(retval) = c("a", "b", "k", "lambda", "c", "value", "convergence")
	}

	# We sort it by value
	sortedIdx = sort.list(retval$value, decreasing=FALSE)
	retval = retval[sortedIdx,]

	# Now perform cross validation using the best parameters as initial conditions
	bestResult = retval[nrow(retval),];
	if(useC == FALSE){
		params = bestResult[1,1:4];
		cross  = cross.validation(params, samples, likelihood, useC=F, lower=lower, upper=upper, method="L-BFGS-B");
	} else {
		params = bestResult[1,1:5];
		cross  = cross.validation(params, samples, likelihood, useC=T, lower=lower, upper=upper, method="L-BFGS-B");
	}

	return(list(results=retval, cross=cross));
}

ollgengamma.lines = function(samples, params, useC=FALSE, ...){
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
	y = dollggamma(x, a=params[1], b=params[2], k=params[3], lambda=params[4])
	
	if(useC)
		x = x + params[length(params)]

	lines(x, y, ...)
}