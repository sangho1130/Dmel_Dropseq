# From the <Drop-seq alignment cookbook> provided by McCarroll Lab.

a <- read.table("read count table", header = FALSE, stringsAsFactors = FALSE)
x <- cumsum(a$V1)
x <- x/max(x)
plot(1:length(x), x, type='l', col="blue", 
     xlab = "cell barcodes sorted by number of reads\n[descending]",
     ylab = "cumulative fraction of reads", 
     xlim = c(1,5000)) 
abline(v = XXX, col = 'red', lty = 2) # XXX to a cutoff value
