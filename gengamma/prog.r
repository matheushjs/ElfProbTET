graphics.off();
dev.new(height=6, width=8);

ALPHA = 1
BETA  = 0.7
GAMMA = 10

means = c();
stds  = c();

means2 = c();
stds2  = c();

allSizes = 2**seq(3, 17, length=100);

for(size in allSizes){
	estims = c();
	estims2 = c();
	for(iter in 1:100){
		sample = rggamma(size, ALPHA, BETA, GAMMA);
		estim = min(sample) - sd(sample) * sqrt(log(log(size)) / (2*size));
		estims = c(estims, estim);
		if(estim > 0){
			estims2 = c(estims2, pggamma(estim, ALPHA, BETA, GAMMA));
		} else {
			estims2 = c(estims2, 0);
		}
	}
	means = c(means, mean(estims));
	stds  = c(stds, sd(estims));

	means2 = c(means2, mean(estims2));
	stds2  = c(stds2, sd(estims2));
}

geom = rbind(
	cbind(allSizes, smooth(means + stds)),
	cbind(rev(allSizes), smooth(rev(means - stds)))
	);

geom2 = rbind(
	cbind(allSizes, smooth(means2 + stds2)),
	cbind(rev(allSizes), smooth(rev(means2 - stds2)))
	);

plot(allSizes, means, pch=19, ylim=c(0, max(means) + 2), cex=1.2, xlab="Sample Size", ylab="Estimator Average", log="x");

polygon(geom, col="#0000FF44", border="blue");
savePlot("consistency-illustration3.png");

dev.new(height=6, width=8);

plot(allSizes, means2, ylim=c(0 - 0.02, max(means2[is.finite(means2)]) + 0.02), pch=19, cex=1.2, xlab="Sample Size", ylab="Cumulative Probability", log="x");
polygon(geom2, col="#FF000044", border="red");
abline(h=0);

savePlot("consistency-illustration3-cdf.png");
