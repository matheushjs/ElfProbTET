require(microbenchmark)
require(elfDistr)
source("kwcwg.r")

setOption("width", 150)

x = seq(0.1, 2, length=150000)

result = microbenchmark(
	kwcwg.pdf(1, 0.5, 1, 1, 1, 1),
	dkwcwg(1, 0.5, 1, 1, 1, 1),
	unit="ms",
	control=list(warmup=1000)
)
cat("\n\n============================\n")
print(result)

result = microbenchmark(
	kwcwg.pdf(x, 0.5, 1, 1, 1, 1),
	dkwcwg(x, 0.5, 1, 1, 1, 1),
	unit="s",
	control=list(warmup=10)
)
cat("\n\n============================\n")
print(result)
