
kwcwg.pdf = function(x, alpha, beta, gamma, a, b){
	if(   x < 0
	   || alpha < 0 || alpha > 1
	   || beta < 0
	   || gamma < 0
	   || a < 0
	   || b < 0
	){
		cat("Must respect:
		x > 0
		0 < alpha < 1
		beta > 0
		gamma > 0
		a > 0
		b > 0");
		return(0);
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
	logDE = (b-1)*log(1 - D/(E + 1e-10))

	return(exp(logA + logB - logC + logDE))
}
