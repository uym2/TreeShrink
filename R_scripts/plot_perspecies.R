args = commandArgs(TRUE)
infile = args[1]
outfile = args[2]

require(ggplot2)

d = read.table(infile,header=F)
p <- ggplot(d,aes(x=log(V1))) + 
    geom_histogram(bins=100,aes(y=..density..)) +
    xlab("signature") + theme_classic()
ggsave(outfile,p)
