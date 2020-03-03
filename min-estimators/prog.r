
source("myoptim.r");

estimator.log10 = function(data){
	(1 - (sd(data) / (mean(data) * log10(length(data))))) * min(data)
}

estimator.mcdiarmid = function(data){
	sqrt(-log(0.01 / 2) / (2*n));
}

gamma.test = function(estimator, sampleSize=100, numIter=100){
	results = NULL;

	for(shape in c(0.25, 0.5, 1, 2, 5, 10, 20)){
	for(scale in c(0.25, 0.5, 1, 2, 5, 10, 20)){
		squares = rep(0, numIter);
		for(i in seq(numIter)){
			realMin = runif(n=1, min=0, max=100);
			data = realMin + rgamma(n=sampleSize, shape=shape, scale=scale);
			estimatedMin = estimator.log10(data);
			squares[i] = (realMin - estimatedMin)**2;
			print(c(realMin, realMin, estimatedMin));
		}
		results = rbind(results, c(shape, scale, sqrt(mean(squares))));
		print(tail(results));
	}
	}
	results;
}

results = gamma.test(estimator.mcdiarmid);
