require(rmutil);

source("inf_model.r");

allModels = list(
	InfModel(
		name        = "ExpWeibull",
		paramNames  = c("shape", "scale", "family"),
		lowerBounds = c(1e-10, 1e-10, 1e-10),
		upperBounds = c(Inf, Inf, Inf),
		initialParams = list(c(0.02, 0.5, 2, 5, 10),
		                     c(0.02, 0.5, 2, 5, 10),
		                     c(0.02, 0.5, 2, 5, 10)),
		param_pdf = function(samples, p, log, ...){
			allLogs = rep(log(1e-300), length(samples));
			retval = try(rmutil::dgweibull(samples, s=p[1], m=p[2], f=p[3], log=log, ...));
			if(is.numeric(retval))
				allLogs[samples > 0] = retval;
			if(log == FALSE)
				allLogs = exp(allLogs);
			allLogs;
		}
	)
);
