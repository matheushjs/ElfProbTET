
source("myoptim.r")

InfModel = function(
	name,          # = "Gamma"
	paramNames,    # = c("shape", "scale")
	lowerBounds,   # = c(1e-10, 1e-10)
	upperBounds,   # = c(Inf, Inf)
	initialParams, # = list(c(0.02, 0.5, 2, 5, 10), c(0.02, 0.5, 2, 5, 10))
	param_pdf      # = function(samples, p, ...) dgamma(samples, shape=p[1], scale=p[2], ...)
){
	if(!is.character(name)) stop("name must be string.");
	if(!is.character(paramNames)) stop("paramNames must be string.");
	if(!is.numeric(lowerBounds) || !is.numeric(upperBounds)) stop("lower and upperBounds must be numeric vectors.");
	if(length(paramNames) != length(lowerBounds) || length(lowerBounds) != length(upperBounds)) stop("Length of parameter vectors differ.");
	if(!is.list(initialParams)) stop("initialParams must be list.");
	for(i in 1:length(initialParams))
		if(!is.numeric(initialParams[[1]])) stop("InitialParams must be list of numerical vectors.");
	if(!is.function(param_pdf)) stop("param_pdf must be function.");

	structure(list(name=name, paramNames=paramNames, lowerBounds=lowerBounds,
				   upperBounds=upperBounds, initialParams=initialParams,
				   param_pdf=param_pdf),
			  class="InfModel");
};

# Generic for inference
infer = function(x, ...) UseMethod("infer")

infer.InfModel = function(model, samples, useC=FALSE){
	estimatedC = min(samples) * 0.995;

	# The likelihood function
	likelihood = function(ourSamples, params, C=0){
		ourSamples = ourSamples - C;

		# dgweibull does not accept a sample of value 0
		allLogs = model$param_pdf(ourSamples, params, log=TRUE);

		problems = which(!is.finite(allLogs))
		allLogs[problems] = log(1e-300); # Merely to force optim to continue optimizing
		if(length(problems) > 0 && length(problems) < 5){
			warning(paste(model$name, ": Low amount (<5) of warnings at points:", ourSamples[problems]), call.=FALSE)
		}

		theSum = -sum(allLogs) # Invert so that minimization yields a maximum
		
		return(theSum)
	}

	retval = NULL;

	if(useC){
		model$paramNames = c(model$paramNames, "c");
		model$lowerBounds = c(model$lowerBounds, 1e-10);
		model$upperBounds = c(model$upperBounds, min(samples) - 1e-10);
		model$initialParams[[length(model$initialParams) + 1]] = c(estimatedC);
	}

	# Recursive function to iterate over the initial parameters
	func = function(initialParams, curParams, curIdx){
		if(curIdx > length(initialParams)){
			print(curParams);

			if(useC == FALSE){
				result = myoptim(curParams, function(p) likelihood(samples, p),
								 lower=model$lowerBounds, upper=model$upperBounds, method="L-BFGS-B");
			} else {
				result = myoptim(curParams, function(p) likelihood(samples, p, p[length(p)]),
								 lower=model$lowerBounds, upper=model$upperBounds, method="L-BFGS-B");
			}
			
			params = result$par;
			val = -result$value; # Undo signal invertion in the likelihood function
			convergence = result$convergence;
			# cat("Got params:", params, "\n")

			# Access to parent scope
			retval <<- rbind(retval, c(params, val, convergence));
		} else {
			initParams = initialParams[[curIdx]];
			if(is.function(initParams))
				initParams = initParams(curParams);
			
			for(p in initParams)
				func(initialParams, c(curParams, p), curIdx+1);
		}
	}
	func(model$initialParams, NULL, 1);

	retval = as.data.frame(retval);
	colnames(retval) = c(model$paramNames, "value", "convergence")

	# We sort it by value
	sortedIdx = sort.list(retval$value, decreasing=FALSE);
	retval = retval[sortedIdx,];

	# Now perform cross validation using the best parameters as initial conditions
	bestResult = retval[nrow(retval),];
	params = bestResult[1,1:length(model$paramNames)];
	cross  = cross.validation(params, samples, likelihood, useC=useC,
							  lower=model$lowerBounds, upper=model$upperBounds, method="L-BFGS-B");

	return(list(results=retval, cross=cross));
}

# Lines is already a generic
lines.InfModel = function(model, samples, params, useC=FALSE, ...){
	delta = diff(quantile(samples, c(0.05, 0.95)));
	minVal = min(samples) - 0.95*delta;
	maxVal = max(samples) + 1.05*delta;
	
	if(minVal <= 0)
		minVal = 1e-100;

	if(useC){
		minVal = minVal - params[length(params)];
		maxVal = maxVal - params[length(params)];
	}

	x = seq(minVal, maxVal, length=1000);
	x[x <= 1e-10] = 1e-10;
	y = model$param_pdf(x, params);
	if(!is.numeric(y))
		y = rep(0, length(x));

	if(useC)
		x = x + params[length(params)];

	lines(x, y, ...);
}

# For debugging:
model = InfModel("Gamma", c("shape", "scale"), c(1-10, 1-10), c(Inf, Inf), list(0.02, 0.02), function(p, ...) dgamma(shape=p[1], shape=p[2], ...))
