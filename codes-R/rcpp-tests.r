require(Rcpp)

setOption("width", 150)

f = function(x){
	A = x**2
	B = exp(A)
	C = sqrt(B)
	D = log(C)
	E = D**10
	f = log2(E)
	return(f);
}

Rcpp::cppFunction('
	NumericVector g(NumericVector v){
		for(auto &num: v){
			double A = std::pow(num, 2);
			double B = std::exp(A);
			double C = std::sqrt(B);
			double D = std::log(C);
			double E = std::pow(D, 10);
			double F = std::log2(E);
			num = F;
		}
		return v;
	}', verbose=T)

h = function(x){
	for(i in 1:length(x)){
		A = x[i]**2
		B = exp(A)
		C = sqrt(B)
		D = log(C)
		E = D**10
		x[i] = log2(E)
	}
	return(x);
}

# This is terrible.
#i = function(x){
#	return(apply(as.matrix(x), 1, function(elem){
#			A = elem**2
#			B = exp(A)
#			C = sqrt(B)
#			D = log(C)
#			E = D**10
#			return(log2(E))
#		}))
#}

x = 1:1e5
result = microbenchmark::microbenchmark(f(x), g(x), h(x), unit="ms")


print(result)
