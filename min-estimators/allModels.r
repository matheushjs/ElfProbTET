require(rmutil);
require(truncnorm);
require(ollggamma);
require(elfDistr);
require(ggamma);

source("inf_model.r");

allModels = function(samples){
	list(
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
		),
		InfModel(
			name        = "Gamma",
			paramNames  = c("shape", "scale"),
			lowerBounds = c(1e-10, 1e-10),
			upperBounds = c(Inf, Inf),
			initialParams = list(c(0.5, 2, 5, 10, 20),
								 mean(samples) / c(0.5, 2, 5, 10, 20)),
			param_pdf = function(samples, p, ...) dgamma(samples, shape=p[1], scale=p[2], ...)
		),
		InfModel(
			name        = "Weibull",
			paramNames  = c("shape", "scale"),
			lowerBounds = c(1e-7, 1e-7),
			upperBounds = c(Inf, Inf),
			initialParams = list(c(0.5, 2, 5, 10, 20, 30, 40),
								 quantile(samples, p=.632) * c(0.666, 1, 1.5)),
			param_pdf = function(samples, p, ...) dweibull(samples, shape=p[1], scale=p[2], ...)
		),
		InfModel(
			name        = "Normal",
			paramNames  = c("mean", "sd"),
			lowerBounds = c(-Inf, 1e-10),
			upperBounds = c(Inf, Inf),
			initialParams = list(mean(samples) * c(0.666, 1, 1.5),
			                     sd(samples) * c(0.666, 1, 1.5)),
			param_pdf = function(samples, p, ...) dnorm(samples, mean=p[1], sd=p[2], ...)
		),
		InfModel(
			name        = "Lognormal",
			paramNames  = c("mean", "sd"),
			lowerBounds = c(-Inf, 1e-10),
			upperBounds = c(Inf, Inf),
			initialParams = list(mean(log(samples)) * c(0.666, 1, 1.5),
			                     sd(log(samples)) * c(0.666, 1, 1.5)),
			param_pdf = function(samples, p, ...) dlnorm(samples, meanlog=p[1], sdlog=p[2], ...)
		),
		InfModel(
			name        = "T.Normal",
			paramNames  = c("mean", "sd"),
			lowerBounds = c(-Inf, 1e-10),
			upperBounds = c(Inf, Inf),
			initialParams = list(mean(samples) * c(0.666, 1, 1.5),
			                     sd(samples) * c(0.666, 1, 1.5)),
			param_pdf = function(samples, p, log, ...){
				res = dtruncnorm(samples, a=0, mean=p[1], sd=p[2], ...);
				if(log) res = log(res);
				res;
			}
		),
		InfModel(
			name        = "G.Gamma",
			paramNames  = c("shape", "scale", "k"),
			lowerBounds = c(1e-7, 1e-7, 1e-7),
			upperBounds = c(Inf, Inf, Inf),
			initialParams = list(c(0.5, 2, 5, 10),
			                     c(0.5, 2, 5, 10),
			                     c(0.5, 2, 5, 10)),
			param_pdf = function(samples, p, ...) ggamma::dggamma(samples, a=p[1], b=p[2], k=p[3], ...)
		),
		InfModel(
			name        = "Kw-CWG",
			paramNames  = c("alpha", "beta", "gamma", "a", "b"),
			lowerBounds = c(1e-10, 1e-10, 1e-10, 1e-10, 1e-10),
			upperBounds = c(1, Inf, Inf, Inf, Inf),
			initialParams = list(c(0.1, 0.9),
			                     quantile(samples, p=.632) * c(0.666, 1, 1.5),
			                     c(2, 10),
			                     c(2, 10),
			                     c(0.1, 10)),
			param_pdf = function(samples, p, ...) dkwcwg(samples, p[1], p[2], p[3], p[4], p[5], ...)
		),
		InfModel(
			name        = "OLL-GG",
			paramNames  = c("a", "b", "k", "lambda"),
			lowerBounds = c(1e-7, 1e-7, 1e-7, 1e-7, 1e-7),
			upperBounds = c(Inf, Inf, Inf, Inf),
			initialParams = list(c(0.5, 2, 5),
			                     c(0.5, 2, 5),
			                     c(0.5, 5),
			                     c(0.25, 1.5)),
			param_pdf = function(samples, p, ...) dollggamma(samples, a=p[1], b=p[2], k=p[3], lambda=p[4], ...)
		)
	);
}
