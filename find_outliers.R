outliers <- function(x,thres=10) {m=median(x); ut = m + thres*median(abs(x-m)); x>ut;}

args <- commandArgs(TRUE)
infile = args[1]
thres = args[2]

x <- read.table(infile)$V1
x[outliers(x,thres=15)]
