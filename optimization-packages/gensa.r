require(GenSA)

f = function(x){
	polynomial = (x+8)*(x+2)*x*(x-4)*(x-8)

	return(-polynomial**2 * exp(-x**2/10))
}

x = seq(-10, 10, length=2000)

plot(x, f(x), type="l")


initialValues = c(-8, -4, -1.5, 0, 0.5, 3, 5)

for(i in 1:length(initialValues)){
	init = initialValues[i]
	result = GenSA(c(init), f, c(-200), c(200))

	print(names(result))

	points(init, f(init), col=i, pch=19)
	points(result$par[1], result$value, col=i, pch=19)

	arrows(init, f(init), result$par[1], result$value, lwd=2, length=0.2)
}

savePlot("gensa.png")
