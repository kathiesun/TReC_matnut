args <- commandArgs(trailingOnly = TRUE)
x = rnorm(n=as.numeric(args[1]), mean=as.numeric(args[2]))
print(x)
