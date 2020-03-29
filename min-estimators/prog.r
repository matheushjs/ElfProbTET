
source("myoptim.r");

estimator.log10 = function(data){
	(1 - (sd(data) / (mean(data) * log10(length(data))))) * min(data)
}

estimator.mcdiarmid = function(data){
	min(data) * (1 - sqrt(-log(0.01 / 2) / (2*length(data))));
}

gamma.test = function(estimator, sampleSize=100, numIter=100){
	results = NULL;

	# The likelihood function
	likelihood = function(ourSamples, params, C=0){
		ourSamples = ourSamples - C

		allLogs = dgamma(ourSamples, shape=params[1], scale=params[2], log=TRUE)

		problems = which(!is.finite(allLogs))
		allLogs[problems] = log(1e-300); # Merely to force optim to continue optimizing
		if(length(problems) > 0 && length(problems) < 5){
			warning(paste("gamma: Low amount (<5) of warnings at points:", ourSamples[problems]), call.=FALSE)
		}

		theSum = -sum(allLogs) # Invert so that minimization yields a maximum
		
		return(theSum)
	}

	for(shape in c(0.25, 0.5, 1, 2, 5, 10, 20, 40)){
	for(scale in c(0.25, 0.5, 1, 2, 5, 10, 20, 40)){
		table = NULL;
		for(i in 1:numIter){
			data = rgamma(n=sampleSize, shape=shape, scale=scale);
			estimatedMin = estimator(data);
			subData = data - estimatedMin;

			estimatedShape = shape
			estimatedScale = mean(subData) / estimatedShape

			result1 = myoptim(c(estimatedShape, estimatedScale), function(p) likelihood(subData, p), lower=c(0, 0), upper=c(Inf, Inf), method="L-BFGS-B")
			result2 = myoptim(c(shape, scale), function(p) likelihood(data, p), lower=c(0, 0), upper=c(Inf, Inf), method="L-BFGS-B")

			#print(c(estimatedMin, -result$value, result$par));
			table = rbind(table, c(estimatedMin, -result1$value, -result2$value));
		}
		results = rbind(results, colMeans(table));
		print(tail(results));
	}
	}
	results;
}

results = gamma.test(estimator.mcdiarmid);
