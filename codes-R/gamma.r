source("myoptim.r")

# @param useC Tells us to also estimate parameter C, which is the amount to subtract from the samples.
gamma.infer = function(samples, useC=FALSE){
	estimatedC    = min(samples) * 0.995

	# The likelihood function
	likelihood = function(ourSamples, params, C=0){
		ourSamples = ourSamples - C

		allLogs = dgamma(ourSamples, shape=params[1], scale=params[2], log=TRUE)

		problems = which(!is.finite(allLogs))
		if(length(problems) > 0 && length(problems) < 5){
			warning(paste("gamma: Low amount (<5) of warnings at points:", ourSamples[problems]), call.=FALSE)
		}

		theSum = -sum(allLogs) # Invert so that minimization yields a maximum
		
		return(theSum)
	}

	retval = NULL

	if(useC == FALSE){
		lower = c(1e-10, 1e-10)
		upper = c(Inf, Inf)
	} else {
		lower = c(1e-10, 1e-10, 1e-10)
		upper = c(Inf, Inf, min(samples) - 1e-10)
	}

	shapeList = c(0.5, 2, 5, 10, 20, 100, 200, 500, 1000, 5000, 10000, 20000)
	for(shape in shapeList){
		# The mean of a gamma is shape*scale
		# So we estimate the scale as being scale = mean(samples) / shape
		estimatedScale = mean(samples) / shape
		scaleList = c(estimatedScale * 0.666, estimatedScale, estimatedScale*1.5)

		for(scale in scaleList){
			params = c(shape, scale)
			if(useC)
				params = c(params, estimatedC)
			# print(params)
			
			#cat("Optimizing with initial params:", params, "\n")
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
			#cat("Got params:", params, "\n")

			retval = rbind(retval, c(params, val, convergence))
		}
	}

	retval = as.data.frame(retval)
	if(useC == FALSE){
		colnames(retval) = c("shape", "scale", "value", "convergence")
	} else {
		colnames(retval) = c("shape", "scale", "c", "value", "convergence")
	}

	# We sort it by value
	sortedIdx = sort.list(retval$value, decreasing=FALSE)
	retval = retval[sortedIdx,]

	# Now perform cross validation using the best parameters as initial conditions
	bestResult = retval[nrow(retval),];
	if(useC == FALSE){
		params = bestResult[1,1:2];
		cross  = cross.validation(params, samples, likelihood, useC=F, lower=lower, upper=upper, method="L-BFGS-B");
	} else {
		params = bestResult[1,1:3];
		cross  = cross.validation(params, samples, likelihood, useC=T, lower=lower, upper=upper, method="L-BFGS-B");
	}

	return(list(results=retval, cross=cross));
}

gamma.lines = function(samples, params, useC=FALSE, ...){
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
	y = dgamma(x, shape=params[1], scale=params[2])

	if(useC)
		x = x + params[length(params)]

	lines(x, y, ...)
}
