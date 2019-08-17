
kwcwg.pdf = function(x, alpha, beta, gamma, a, b){
	#cat("Called with", "(", x, ")", alpha, beta, gamma, a, b, "\n")

	# If parameters are not within valid range, return 0
	if(  sum(x < 0) > 0
	  || alpha < 0 || alpha > 1
	  || beta < 0
	  || gamma < 0
	  || a < 0
	  || b < 0
	){
		return(rep(0, length(x)));
	}

	# Original function
	# return(
	#	alpha**a * beta * gamma * a * b * (gamma * x)**(beta-1) * exp(-(gamma*x)**beta) *
	#	(
	#		(1 - exp(-(gamma*x)**beta))**(a-1) /
	#		(alpha + (1 - alpha)*exp(-(gamma*x)**beta))**(a+1)
	#	) *
	#	(
	#		1 -
	#		(alpha**a*(1 - exp(-(gamma*x)**beta))**a) /
	#		(alpha + (1-alpha)*exp(-(gamma*x)**beta))**a
	#	)**(b-1)
	#)

	# A common term in the equation
	aux1 = exp(-(gamma*x)**beta)
	
	# Here we will factor f(x) as being A * (B/C) * (1 - D/E)^(b-1)
	A = alpha**a * beta * gamma * a * b * (gamma * x)**(beta-1) * aux1
	B = 1 - aux1
	C = alpha + (1 - alpha)*aux1
	D = (alpha**a * (1 - aux1)**a)
	E = (alpha + (1-alpha)*aux1)**a

	# cat("A", A, "\n")
	# cat("B", B, "\n")
	# cat("C", C, "\n")
	# cat("D", D, "\n")
	# cat("E", E, "\n")

	result = A * (B**(a-1)/C**(a+1)) * (1 - D/E)**(b-1)

	return(result)

	# Old way
	# logA = log(A)
	# logB = (a-1)*log(B)
	# logC = (a+1)*log(C)
	# logDE = (b-1)*log(1 - D/(E))

	# result = exp(logA + logB - logC + logDE)

	# Set to 0 all results that are not finite
	# result[!is.finite(result)] = 0

	#return(result)
}

# Get the parameters of the kwcwg after fitting the samples
kwcwg.infer = function(samples){
	# This quantile is apparently a good estimator for the beta parameter
	estimatedBeta = quantile(samples, p=.632)

	# The likelihood function
	likelihood = function(params){
		allLogs = log(kwcwg.pdf(samples, params[1], params[2], params[3], params[4], params[5]))

		# Set all NA and Inf to 0
		allLogs[!is.finite(allLogs)] = 0

		theSum = -sum(allLogs)
		
		return(theSum)
	}

	retval = NULL

	# We use a grid of initial values, and take the best of them
	#for(alpha in c(0.1, 0.2, 0.4, 0.6, 0.8, 0.9)){
	# for(beta in c(estimatedBeta*0.666, estimatedBeta, estimatedBeta*1.5)){
	#for(gamma in c(0.1, 2, 6, 10)){
	#for(a in c(0.1, 2, 6, 10)){
	#for(b in c(0.1, 2, 6, 10)){
	for(alpha in c(0.1, 0.9)){
	for(beta in c(estimatedBeta*0.666, estimatedBeta, estimatedBeta*1.5)){
	for(gamma in c(2, 10)){
	for(a in c(2, 10)){
	for(b in c(0.1, 10)){
		params = c(alpha, beta, gamma, a, b)
		# print(params)
		
		cat("Optimizing with initial params:", params, "\n")
		result = optim(params, likelihood)
		result = optim(result$par, likelihood)
		params = result$par
		val = result$value
		cat("Got params:", params, "\n")

		retval = rbind(retval, c(params, val))
	} } } } }

	retval = as.data.frame(retval)
	colnames(retval) = c("alpha", "beta", "gamma", "a", "b", "value")

	# We sort it by value
	sortedIdx = sort.list(retval$value, decreasing=TRUE)
	retval = retval[sortedIdx,]

	return(retval)
}

kwcwg.plot = function(samples, alpha, beta, gamma, a, b){
	minVal = min(samples) * 0.9;
	maxVal = max(samples) * 1.1;

	hist(samples, prob=T, xlim=c(minVal, maxVal))
	x = seq(minVal, maxVal, length=200)
	y = kwcwg.pdf(x, alpha, beta, gamma, a, b)
	print(y)
	lines(x, y)
}
