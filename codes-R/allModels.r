require(rmutil);
require(truncnorm);
require(ollggamma);
require(elfDistr);
require(ggamma);

source("inf_model.r");

allModels = function(samples){
	list(
		InfModel(
			name        = "Gamma",
			paramNames  = c("shape", "scale"),
			lowerBounds = c(1e-10, 1e-10),
			upperBounds = c(Inf, Inf),
			initialParams = list(c(0.5, 2, 5, 10, 20, 100, 200, 500, 1000, 5000, 10000, 20000),
							     function(curParams) (mean(samples) / curParams[1]) * c(0.666, 1, 1.5)),
			param_pdf = function(samples, p, ...) dgamma(samples, shape=p[1], scale=p[2], ...),
			param_cdf = function(samples, p, ...) pgamma(samples, shape=p[1], scale=p[2], ...),
			param_q   = function(q, p, ...) qgamma(q, shape=p[1], scale=p[2], ...)
		),
		InfModel(
			name        = "Weibull",
			paramNames  = c("shape", "scale"),
			lowerBounds = c(1e-7, 1e-7),
			upperBounds = c(Inf, Inf),
			initialParams = list(c(0.5, 2, 5, 10, 20, 30, 40),
								 quantile(samples, p=.632) * c(0.666, 1, 1.5)),
			param_pdf = function(samples, p, ...) dweibull(samples, shape=p[1], scale=p[2], ...),
			param_cdf = function(samples, p, ...) pweibull(samples, shape=p[1], scale=p[2], ...),
			param_q   = function(q, p, ...) qweibull(q, shape=p[1], scale=p[2], ...)
		),
		InfModel(
			name        = "Normal",
			paramNames  = c("mean", "sd"),
			lowerBounds = c(-Inf, 1e-10),
			upperBounds = c(Inf, Inf),
			initialParams = list(mean(samples) * c(0.666, 1, 1.5),
			                     sd(samples) * c(0.666, 1, 1.5)),
			param_pdf = function(samples, p, ...) dnorm(samples, mean=p[1], sd=p[2], ...),
			param_cdf = function(samples, p, ...) pnorm(samples, mean=p[1], sd=p[2], ...),
			param_q   = function(q, p, ...) qnorm(q, mean=p[1], sd=p[2], ...)
		),
		InfModel(
			name        = "T.Normal",
			paramNames  = c("mean", "sd"),
			lowerBounds = c(-Inf, 1e-10),
			upperBounds = c(Inf, Inf),
			initialParams = list(mean(samples) * c(0.666, 1, 1.5),
			                     sd(samples) * c(0.666, 1, 1.5)),
			param_pdf = function(samples, p, log=FALSE, ...){
				res = dtruncnorm(samples, a=0, mean=p[1], sd=p[2], ...);
				if(log) res = log(res);
				res;
			},
			param_cdf = function(samples, p, log=FALSE, ...){
				res = ptruncnorm(samples, a=0, mean=p[1], sd=p[2], ...);
				if(log) res = log(res);
				res;
			},
			param_q = function(q, p, ...) qtruncnorm(q, a=0, mean=p[1], sd=p[2], ...)
		),
		InfModel(
			name        = "Lognormal",
			paramNames  = c("mean", "sd"),
			lowerBounds = c(-Inf, 1e-10),
			upperBounds = c(Inf, Inf),
			initialParams = list(mean(log(samples)) * c(0.666, 1, 1.5),
			                     sd(log(samples)) * c(0.666, 1, 1.5)),
			param_pdf = function(samples, p, ...) dlnorm(samples, meanlog=p[1], sdlog=p[2], ...),
			param_cdf = function(samples, p, ...) plnorm(samples, meanlog=p[1], sdlog=p[2], ...),
			param_q = function(q, p, ...) qlnorm(q, meanlog=p[1], sdlog=p[2], ...)
		),
		InfModel(
			name        = "OLL-GG",
			paramNames  = c("a", "b", "k", "lambda"),
			lowerBounds = c(1e-7, 1e-7, 1e-7, 1e-7),
			upperBounds = c(Inf, Inf, Inf, Inf),
			initialParams = list(c(0.5, 2, 5),
			                     c(0.5, 2, 5),
			                     c(0.5, 5),
			                     c(0.25, 1.5)),
			param_pdf = function(samples, p, ...) dollggamma(samples, a=p[1], b=p[2], k=p[3], lambda=p[4], ...),
			param_cdf = function(samples, p, ...) pollggamma(samples, a=p[1], b=p[2], k=p[3], lambda=p[4], ...),
			param_q   = function(q, p, ...) qollggamma(q, a=p[1], b=p[2], k=p[3], lambda=p[4], ...)
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
			param_pdf = function(samples, p, ...) dkwcwg(samples, p[1], p[2], p[3], p[4], p[5], ...),
			param_cdf = function(samples, p, ...) pkwcwg(samples, p[1], p[2], p[3], p[4], p[5], ...),
			param_q   = function(q, p, ...) qkwcwg(q, p[1], p[2], p[3], p[4], p[5], ...)
		),
		InfModel(
			name        = "G.Gamma",
			paramNames  = c("shape", "scale", "k"),
			lowerBounds = c(1e-7, 1e-7, 1e-7),
			upperBounds = c(Inf, Inf, Inf),
			initialParams = list(c(0.5, 2, 5, 10),
			                     c(0.5, 2, 5, 10),
			                     c(0.5, 2, 5, 10)),
			param_pdf = function(samples, p, ...) ggamma::dggamma(samples, a=p[1], b=p[2], k=p[3], ...),
			param_cdf = function(samples, p, ...) ggamma::pggamma(samples, a=p[1], b=p[2], k=p[3], ...),
			param_q   = function(q, p, ...) ggamma::qggamma(q, a=p[1], b=p[2], k=p[3], ...)
		),
		InfModel(
			name        = "E.Weibull",
			paramNames  = c("shape", "scale", "family"),
			lowerBounds = c(1e-10, 1e-10, 1e-10),
			upperBounds = c(Inf, Inf, Inf),
			initialParams = list(c(0.02, 0.5, 2, 5, 10),
								 c(0.02, 0.5, 2, 5, 10),
								 c(0.02, 0.5, 2, 5, 10)),
			param_pdf = function(samples, p, log=FALSE, ...){
				allLogs = rep(log(1e-300), length(samples));
				retval = try(rmutil::dgweibull(samples[samples > 0], s=p[1], m=p[2], f=p[3], log=TRUE, ...));
				if(is.numeric(retval))
					allLogs[samples > 0] = retval;
				if(log == FALSE)
					allLogs = exp(allLogs);
				allLogs;
			},
			param_cdf = function(samples, p, ...){
				result = rep(0, length(samples));
				probs = try(rmutil::pgweibull(samples[samples > 0], s=p[1], m=p[2], f=p[3], ...));
				if(is.numeric(probs))
					result[samples > 0] = probs[samples > 0];
				result;
			},
			param_q = function(q, p, ...){
				result = rep(0, length(q));
				quantiles= try(rmutil::qgweibull(q, s=p[1], m=p[2], f=p[3], ...));
				if(is.numeric(quantiles))
					result = quantiles;
				result;
			}
		)
	);
}
