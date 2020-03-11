require(rmutil)

source("myoptim.r")

# Parameters are shape s > 0, scale m > 0, family f > 0
# @param useC Tells us to also estimate parameter C, which is the amount to subtract from the samples.
expweibull.infer = function(samples, useC=FALSE){
	estimatedC    = min(samples) * 0.995

	# The likelihood function
	likelihood = function(ourSamples, params, C=0){
		ourSamples = ourSamples - C;

		# dgweibull does not accept a sample of value 0
		retval = try(rmutil::dgweibull(ourSamples[ourSamples > 0], s=params[1], m=params[2], f=params[3], log=TRUE));
		if(is.list(retval)){
			allLogs = retval;
		} else {
			allLogs = rep(log(0), length(ourSamples));
		}

		problems = which(!is.finite(allLogs))
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
		colnames(retval) = c("shape", "scale", "family", "value", "convergence")
	} else {
		colnames(retval) = c("shape", "scale", "family", "c", "value", "convergence")
	}

	# We sort it by value
	sortedIdx = sort.list(retval$value, decreasing=FALSE)
	retval = retval[sortedIdx,]

	# Now perform cross validation using the best parameters as initial conditions
	bestResult = retval[nrow(retval),];
	if(useC == FALSE){
		params = bestResult[1,1:3];
		cross  = cross.validation(params, samples, likelihood, useC=F, lower=lower, upper=upper, method="L-BFGS-B");
	} else {
		params = bestResult[1,1:4];
		cross  = cross.validation(params, samples, likelihood, useC=T, lower=lower, upper=upper, method="L-BFGS-B");
	}

	return(list(results=retval, cross=cross));
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
	x[x <= 1e-10] = 1e-10
	y = try(rmutil::dgweibull(x, s=params[1], m=params[2], f=params[3]));
	if(!is.list(y))
		y = rep(0, length(x));

	if(useC)
		x = x + params[length(params)]

	lines(x, y, ...)
}
