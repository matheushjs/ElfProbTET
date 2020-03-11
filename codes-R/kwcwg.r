require(elfDistr)

source("myoptim.r")

#kwcwg.pdf = function(x, alpha, beta, gamma, a, b){
#	#cat("Called with", "(", x, ")", alpha, beta, gamma, a, b, "\n")
#	
#	# If parameters are not within valid range, return 0
#	if(  sum(x < 0) > 0
#	  || alpha < 0 || alpha > 1
#	  || beta < 0
#	  || gamma < 0
#	  || a < 0
#	  || b < 0
#	){
#		return(rep(0, length(x)));
#	}
#	
#	# Original function
#	# return(
#	#	alpha**a * beta * gamma * a * b * (gamma * x)**(beta-1) * exp(-(gamma*x)**beta) *
#	#	(
#	#		(1 - exp(-(gamma*x)**beta))**(a-1) /
#	#		(alpha + (1 - alpha)*exp(-(gamma*x)**beta))**(a+1)
#	#	) *
#	#	(
#	#		1 -
#	#		(alpha**a*(1 - exp(-(gamma*x)**beta))**a) /
#	#		(alpha + (1-alpha)*exp(-(gamma*x)**beta))**a
#	#	)**(b-1)
#	#)
#	
#	# A common term in the equation
#	aux1 = exp(-(gamma*x)**beta)
#	
#	# Here we will factor f(x) as being A * (B/C) * (1 - D/E)^(b-1)
#	A = alpha**a * beta * gamma * a * b * (gamma * x)**(beta-1) * aux1
#	B = 1 - aux1
#	C = alpha + (1 - alpha)*aux1
#	D = (alpha**a * (1 - aux1)**a)
#	E = (alpha + (1-alpha)*aux1)**a
#	
#	# cat("A", A, "\n")
#	# cat("B", B, "\n")
#	# cat("C", C, "\n")
#	# cat("D", D, "\n")
#	# cat("E", E, "\n")
#	
#	result = A * (B**(a-1)/C**(a+1)) * (1 - D/E)**(b-1)
#	
#	return(result)
#	
#	# Old way
#	# logA = log(A)
#	# logB = (a-1)*log(B)
#	# logC = (a+1)*log(C)
#	# logDE = (b-1)*log(1 - D/(E))
#	
#	# result = exp(logA + logB - logC + logDE)
#	
#	# Set to 0 all results that are not finite
#	# result[!is.finite(result)] = 0
#	
#	#return(result)
#}

# Get the parameters of the kwcwg after fitting the samples
# @param useC Tells us to also estimate parameter C, which is the amount to subtract from the samples.
kwcwg.infer = function(samples, useC=FALSE){
	# This quantile is apparently a good estimator for the beta parameter
	estimatedBeta = quantile(samples, p=.632)
	estimatedC    = min(samples) * 0.995

	# The likelihood function
	likelihood = function(ourSamples, params, C=0){
		ourSamples = ourSamples - C;

		allLogs = dkwcwg(ourSamples, params[1], params[2], params[3], params[4], params[5], log=T)

		problems = which(!is.finite(allLogs))
		if(length(problems) > 0 && length(problems) < 5){
			warning("kwcwg(", params, "): Low amount (<5) of warnings at points: ", ourSamples[problems], call.=FALSE)
		}

		theSum = -sum(allLogs) # Invert so that minimization yields a maximum
		
		return(theSum)
	}

	retval = NULL

	if(useC == FALSE){
		lower = c(1e-10, 1e-10, 1e-10, 1e-10, 1e-10)
		upper = c(1, Inf, Inf, Inf, Inf)
	} else {
		lower = c(1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10)
		upper = c(1, Inf, Inf, Inf, Inf, min(samples) - 1e-10)
	}

	# We use a grid of initial values, and take the best of them
	#for(alpha in c(0.1, 0.2, 0.4, 0.6, 0.8, 0.9))
	#for(beta in c(estimatedBeta*0.666, estimatedBeta, estimatedBeta*1.5))
	#for(gamma in c(0.1, 2, 6, 10))
	#for(a in c(0.1, 2, 6, 10))
	#for(b in c(0.1, 2, 6, 10)){
	for(alpha in c(0.1, 0.9))
	for(beta in c(estimatedBeta*0.666, estimatedBeta, estimatedBeta*1.5))
	for(gamma in c(2, 10))
	for(a in c(2, 10))
	for(b in c(0.1, 10)){
		params = c(alpha, beta, gamma, a, b)
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
		colnames(retval) = c("alpha", "beta", "gamma", "a", "b", "value", "convergence")
	} else {
		colnames(retval) = c("alpha", "beta", "gamma", "a", "b", "c", "value", "convergence")
	}

	# We sort it by value
	sortedIdx = sort.list(retval$value, decreasing=FALSE)
	retval = retval[sortedIdx,]

	# Now perform cross validation using the best parameters as initial conditions
	bestResult = retval[nrow(retval),];
	if(useC == FALSE){
		params = bestResult[1,1:5];
		cross  = cross.validation(params, samples, likelihood, useC=F, lower=lower, upper=upper, method="L-BFGS-B");
	} else {
		params = bestResult[1,1:6];
		cross  = cross.validation(params, samples, likelihood, useC=T, lower=lower, upper=upper, method="L-BFGS-B");
	}

	return(list(results=retval, cross=cross));
}

kwcwg.lines = function(samples, params, useC=FALSE, ...){
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
	y = dkwcwg(x, params[1], params[2], params[3], params[4], params[5])
	y[!is.finite(y)] = 0;
	
	if(useC)
		x = x + params[length(params)]

	lines(x, y, ...)
}
