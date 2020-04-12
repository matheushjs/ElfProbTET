
myoptim = function(params, f, ...){
	result = try(optim(params, f, ...))

	if(is.list(result))
		return(result)
	else
		return(list(par=rep(0, length(params)), value=1e300, convergence=0))
}

cross.validation = function(params, samples, likelihood, useC=F, folds=5, ...){
	foldLen = length(samples) / folds;
	testResults = rep(0, folds);
	samples = sample(samples);

	for(i in 1:folds){
		foldIdx = (foldLen*(i-1) + 1):(foldLen*i);
		test = samples[foldIdx];
		train = samples[-foldIdx];

		if(useC == FALSE){
			result = myoptim(params, function(p) likelihood(train, p, C=0), ...);
			testResults[i] = -likelihood(test, result$par, C=0);
		} else {
			result = myoptim(params, function(p) likelihood(train, p, C=p[length(p)]), ...);
			testResults[i] = -likelihood(test, result$par, C=result$par[length(result$par)]);
		}
	}

	mean(testResults);
}
