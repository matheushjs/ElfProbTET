
myoptim = function(params, f, ...){
	result = try(optim(params, f, ...))

	if(is.list(result))
		return(result)
	else
		return(list(par=rep(0, length(params)), value=1e300, convergence=0))
}
