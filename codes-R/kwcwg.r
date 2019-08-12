
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

	return(
		alpha**a * beta * gamma * a * b * (gamma * x)**(beta-1) * exp(-(gamma*x)**beta) *
		(
			(1 - exp(-(gamma*x)**beta))**(a-1) /
			(alpha + (1 - alpha)*exp(-(gamma*x)**beta))**(a+1)
		) *
		(
			1 -
			(alpha**a*(1 - exp(-(gamma*x)**beta))**a) /
			(alpha + (1-alpha)*exp(-(gamma*x)**beta))**a
		)**(b-1)
	)
}
