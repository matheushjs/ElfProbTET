
kwcwg.pdf = function(x, alpha, beta, gamma, a, b){
	#cat("Called with", "(", x, ")", alpha, beta, gamma, a, b, "\n")

	if(  sum(x < 0) > 0
	  || alpha < 0 || alpha > 1
	  || beta < 0
	  || gamma < 0
	  || a < 0
	  || b < 0
	){
#		cat("Must respect:
#		x > 0
#		0 < alpha < 1
#		beta > 0
#		gamma > 0
#		a > 0
#		b > 0\n");
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
	
	# Here we will factor f(x) as being A * (B/C) * (1 - D/E)^(b-1)
	A = alpha**a * beta * gamma * a * b * (gamma * x)**(beta-1) * exp(-(gamma*x)**beta)
	B = 1 - exp(-(gamma*x)**beta)
	C = alpha + (1 - alpha)*exp(-(gamma*x)**beta)
	D = (alpha**a*(1 - exp(-(gamma*x)**beta))**a)
	E = (alpha + (1-alpha)*exp(-(gamma*x)**beta))**a

	#cat("A", A, "\n")
	#cat("B", B, "\n")
	#cat("C", C, "\n")
	#cat("D", D, "\n")
	#cat("E", E, "\n")

	logA = log(A)
	logB = (a-1)*log(B)
	logC = (a+1)*log(C)
	logDE = (b-1)*log(1 - D/(E))

	result = exp(logA + logB - logC + logDE)

	result[!is.finite(result)] = 0

	return(result)
}

# Get the parameters of the kwcwg after fitting the dataset
kwcwg.infer = function(dataset){

	# The likelihood function
	likelihood = function(params){
		return(-sum(log(
			kwcwg.pdf(dataset, params[1], params[2], params[3], params[4], params[5]) + 1e-100
		)))
	}

	retval = NULL

	# We use a grid of initial values, and take the best of them
#	for(alpha in c(0.1, 0.2, 0.4, 0.6, 0.8, 0.9)){
#	for(beta in c(0.1, 2, 6, 10)){
#	for(gamma in c(0.1, 2, 6, 10)){
#	for(a in c(0.1, 2, 6, 10)){
#	for(b in c(0.1, 2, 6, 10)){
	for(alpha in c(0.1, 0.4, 0.9)){
	for(beta in c(6, 10)){
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


