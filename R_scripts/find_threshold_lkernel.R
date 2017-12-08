suppressMessages(require(BMS))

threshold <- function(y,e=0.05) {x=y[y>0];exp(quantile(density(log(x),adjust=1),p=1-e));}


args = commandArgs(TRUE)

datafile = args[1]
e = as.numeric(args[2])

d = read.table(datafile)
threshold(d$V1,e=e)
