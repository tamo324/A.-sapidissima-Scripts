setwd("~/Desktop/shads")

install.packages("vcfR")
library(vcfR)

vcf <- read.vcfR("filteredshad.vcf")
head(vcf)

#chromR object 
chromr <- create.chromR(vcf, name = "CHROM", seq = NULL, ann = NULL, verbose = TRUE)


#plot 1 first look
plot(chromr) #error figure margins too large 

#plot 2 a few filters 
chromfiltered <- masker(chromr, min_QUAL = 1, max_DP = 50, min_MQ = 40,  max_MQ = 61)
plot(chromfiltered)

#plot 3 include variants 
chromvariants <- proc.chromR(chromfiltered, verbose=TRUE)
plot(chromvariants)

#plot4 
chromoqc(chromvariants, dp.alpha = 66) #Error in graphics::par(userpar) : invalid value specified for graphical parameter "pin"
plot(chromoqc)


chromoqc(chromprocessed, dp.alpha = 66)














