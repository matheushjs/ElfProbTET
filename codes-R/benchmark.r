require(elfDistr)
source("kwcwg.r")

time = system.time( kwcwg.pdf(seq(0.1, 2, length=5000000), 0.5, 1, 1, 1, 1) )
cat("R Implementation: ", time["elapsed"], "\n")
time = system.time( dkwcwg(seq(0.1, 2, length=5000000), 0.5, 1, 1, 1, 1) )
cat("C++ Implementation: ", time["elapsed"], "\n")
