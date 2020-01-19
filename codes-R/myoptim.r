
myoptim = function(params, f, ...){
	result = try(optim(params, f, ...))

	if(is.list(result))
		return(result)
	else
		return(list(par=rep(0, length(params)), value=1e300, convergence=0))
}

cross.validation = function(params, samples, likelihood, C=C, folds=5, ...){
	foldLen = length(samples) / folds;
	testResults = rep(0, folds);

	for(i in 1:folds){
		foldIdx = (foldLen*(i-1) + 1):(foldLen*i);
		test = samples[foldIdx];
		train = samples[-foldIdx];

		result = myoptim(params, function(p) likelihood(train, p, C=C), ...);
		testResults[i] = likelihood(test, result$par, C=C);
	}

	mean(testResults);
}
