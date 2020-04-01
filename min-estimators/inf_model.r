
InfModel = function(
	name,          # = "Gamma"
	paramNames,    # = c("shape", "scale")
	lowerBounds,   # = c(1e-10, 1e-10)
	upperBounds,   # = c(Inf, Inf)
	initialParams, # = list(c(0.02, 0.5, 2, 5, 10), c(0.02, 0.5, 2, 5, 10))
	param_pdf      # = function(p, ...) dgamma(shape=p[1], scale=p[2], ...)
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
infer = function(x) UseMethod("infer")

infer.InfModel = function(samples, useC=FALSE){

}

# Lines is already a generic
lines.InfModel = function(samples, params, useC=FALSE, ...){

}

# For debugging:
# model = InfModel("Gamma", c("shape", "scale"), c(1-10, 1-10), c(Inf, Inf), list(0.02, 0.02), function(p, ...) dgamma(shape=p[1], shape=p[2], ...))
