
source("myoptim.r");

gamma.test = function(){
	for(shape in c(0.25, 0.5, 1, 2, 5, 10, 20)){
	for(scale in c(0.25, 0.5, 1, 2, 5, 10, 20)){
		realMin = runif(n=1, min=0, max=100);
		data = realMin + rgamma(n = 100000, shape=shape, scale=scale);
		estimatedMin = (1 - (sd(data) / mean(data) / log10(length(data)))) * min(data);
		print(c(realMin, min(data), estimatedMin));
	}
	}
}

gamma.test();
