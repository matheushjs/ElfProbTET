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
	}')

x = 1:1e6
result = microbenchmark::microbenchmark(f(x), g(x), unit="ms")


print(result)
