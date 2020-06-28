
source("myoptim.r")

InfModel = function(
	name,          # = "Gamma"
	paramNames,    # = c("shape", "scale")
	lowerBounds,   # = c(1e-10, 1e-10)
	upperBounds,   # = c(Inf, Inf)
	initialParams, # = list(c(0.02, 0.5, 2, 5, 10), c(0.02, 0.5, 2, 5, 10))
	param_pdf,     # = function(samples, p, ...) dgamma(samples, shape=p[1], scale=p[2], ...)
	param_cdf,     # = function(probs, p, ...) pgamma(samples, shape=p[1], scale=p[2], ...)
	param_q        # = function(q, params) qgamma(q, shape=params[1], scale=params[2])
){
	if(!is.character(name)) stop("name must be string.");
	if(!is.character(paramNames)) stop("paramNames must be string.");
	if(!is.numeric(lowerBounds) || !is.numeric(upperBounds)) stop("lower and upperBounds must be numeric vectors.");
	if(length(paramNames) != length(lowerBounds) || length(lowerBounds) != length(upperBounds)) stop("Length of parameter vectors differ.");
	if(!is.list(initialParams)) stop("initialParams must be list.");
	for(i in 1:length(initialParams))
		if(!is.numeric(initialParams[[1]])) stop("InitialParams must be list of numerical vectors.");
	if(!is.function(param_pdf)) stop("param_pdf must be function.");
	if(!is.function(param_cdf)) stop("param_cdf must be function.");
	if(!is.function(param_q)) stop("param_q must be function.");

	structure(list(name=name, paramNames=paramNames, lowerBounds=lowerBounds,
				   upperBounds=upperBounds, initialParams=initialParams,
				   param_pdf=param_pdf, param_cdf=param_cdf, param_q=param_q),
			  class="InfModel");
};

# Generic for inference
infer = function(x, ...) UseMethod("infer");

# Generic for qqplot
ppplot = function(x, ...) UseMethod("ppplot");

infer.InfModel = function(model, samples, useC=FALSE, iteratedC=FALSE){
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

	bestInitParams = NULL;
	bestLikelihood = -Inf;

	# Recursive function to iterate over the initial parameters
	func = function(initialParams, curParams, curIdx){
		if(curIdx > length(initialParams)){
			if(useC == TRUE){
				elapsed = system.time({
					result = myoptim(curParams, function(p) likelihood(samples, p, p[length(p)]),
									 lower=model$lowerBounds, upper=model$upperBounds, method="L-BFGS-B");
				})["elapsed"]
			} else if(iteratedC == TRUE){
				elapsed = system.time({
					result = myoptim(curParams,
									 function(p){
										 c = min(samples) - model$param_q(1 - (1 - 0.5)**(1/length(samples)), p);
										 likelihood(samples, p, c)
									 },
									 lower=model$lowerBounds, upper=model$upperBounds, method="L-BFGS-B");
				})["elapsed"]
			} else {
				elapsed = system.time({
					result = myoptim(curParams, function(p) likelihood(samples, p),
									 lower=model$lowerBounds, upper=model$upperBounds, method="L-BFGS-B");
				})["elapsed"]
			}
			
			params = result$par;
			val = -result$value; # Undo signal invertion in the likelihood function
			convergence = result$convergence;
			# cat("Got params:", params, "\n")

			if(val > bestLikelihood){
				# access to parent scope
				bestInitParams <<- curParams;
				bestLikelihood <<- val;
			}

			# Access to parent scope
			retval <<- rbind(retval, c(params, val, convergence, elapsed));
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
	colnames(retval) = c(model$paramNames, "value", "convergence", "elapsed.inf")

	# We sort it by value
	sortedIdx = sort.list(retval$value, decreasing=FALSE);
	retval = retval[sortedIdx,];

	# Now perform cross validation using the best parameters as initial conditions
	bestResult = retval[nrow(retval),];
	params = bestResult[1,1:length(model$paramNames)];
	cross  = cross.validation(params, samples, likelihood, useC=useC,
							  lower=model$lowerBounds, upper=model$upperBounds, method="L-BFGS-B");

	slice = retval[c("value", "elapsed.inf")];
	slice = slice[which(slice$value / bestResult$value > 0.3),]; # Selects the set of best parameters (or useful parameters)

	# results["value"]: maximized likelihood
	# results["convergence"]: convergence status obtained in each inference
	# results["elapsed.inf"]: time taken to perform each inference
	# inf.time: summary of time taken for the set of all useful parameters we found (those that
	#           yielded a minimally reasonable value of likelihood).
	return(list(
		results  = retval,
		cross    = cross,
		inf.time = list(
				mean = mean(slice$elapsed.inf),
				sd   = sd(slice$elapsed.inf),
				n    = nrow(slice)
		)));
}

# Lines is already a generic
lines.InfModel = function(model, samples, params, useC=FALSE, iteratedC=FALSE, ...){
	delta = diff(quantile(samples, c(0.05, 0.95)));
	minVal = min(samples) - 0.95*delta;
	maxVal = max(samples) + 1.05*delta;
	
	if(minVal <= 0)
		minVal = 1e-100;

	if(useC){
		minVal = minVal - params[length(params)];
		maxVal = maxVal - params[length(params)];
	}
	if(iteratedC){
		c = min(samples) - model$param_q(1 - (1 - 0.5)**(1/length(samples)), params);
		minVal = minVal - c;
		maxVal = maxVal - c;
	}

	x = seq(minVal, maxVal, length=1000);
	if(model$name != "Normal")
		x[x <= 1e-10] = 1e-10;
	y = model$param_pdf(x, params);
	if(!is.numeric(y))
		y = rep(0, length(x));

	if(useC)
		x = x + params[length(params)];
	if(iteratedC){
		c = min(samples) - model$param_q(1 - (1 - 0.5)**(1/length(samples)), params);
		x = x + c;
	}

	lines(x, y, ...);
}

ppplot.InfModel = function(model, samples, params, useC=FALSE, iteratedC=FALSE, ...){
	if(useC == TRUE){
		samples = samples - params[length(params)];
	}
	if(iteratedC){
		c = min(samples) - model$param_q(1 - (1 - 0.5)**(1/length(samples)), params);
		samples = samples - c;
	}

	samples = sort(samples);

	x = model$param_cdf(samples, params);
	y = 1/length(samples) * seq(0, length(samples), length=length(samples)); # ECDF

	ang = -pi/4;
	rot = rbind(
	  c(cos(ang), -sin(ang)),
	  c(sin(ang), cos(ang))
	);

	data = cbind(x, y) %*% t(rot);

	#x = x + 2*(x - y);

	# Reduce number of points
	#data = data[seq(ceiling(runif(n=1, max=15)), length(samples), by=25),]
	fun = splinefun(data[,1], data[,2]);

	x = seq(0, sqrt(2), length=sample(18:22)[1]);
	y = fun(x);

	points(x, y, ...);
	abline(a=0, b=0);
}

# For debugging:
model = InfModel("Gamma", c("shape", "scale"), c(1-10, 1-10), c(Inf, Inf), list(0.02, 0.02), function(isamples, p, ...) dgamma(samples, shape=p[1], scale=p[2], ...), function(samples, p) pgamma(samples, shape=p[1], scale=p[2]), function(q, p) qgamma(q, shape=p[1], scale=p[2]));
