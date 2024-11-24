setwd("~/Desktop/shads")
library(ape)
library(vcfR)
library(tidyverse)


#preliminary data viewing with vcfR:
bobvcf2 <- read.vcfR("shad.filtered.ann.vcf.gz")
#chromR object 
bobchromr <- create.chromR(bobvcf2, name = "CHROM", seq = NULL, ann = NULL, verbose = TRUE)
#plot 1 first look
plot(bobchromr) #error figure margins too large 
#plot 2 a few filters 
chromfilteredbob <- masker(bobchromr, min_QUAL = 1, max_DP = 50, min_MQ = 40,  max_MQ = 61)
plot(chromfilteredbob)
#plot 3 include variants 
chromvariantsbob <- proc.chromR(chromfilteredbob, verbose=TRUE)
plot(chromvariantsbob)
#plot4 
chromoqc(chromvariantsbob, dp.alpha = 66) #Error in graphics::par(userpar) : invalid value specified for graphical parameter "pin"
plot(chromoqc)

#SNP Density plot 
#snp variants from samtools
snpsxscaf <- read_xlsx("samtoolsdata.xlsx")
View(snpsxscaf)
# rename
snpschromo <- snpsxscaf %>%
  mutate(Chromosome = case_when(
    Chromosome == "NC_014690.1" ~ "MT",
    Chromosome == "NC_055957.1" ~ "1",
    Chromosome == "NC_055958.1" ~ "2",
    Chromosome == "NC_055959.1" ~ "3",
    Chromosome == "NC_055960.1" ~ "4",
    Chromosome == "NC_055961.1" ~ "5",
    Chromosome == "NC_055962.1" ~ "6",
    Chromosome == "NC_055963.1" ~ "7",
    Chromosome == "NC_055964.1" ~ "8",
    Chromosome == "NC_055965.1" ~ "9",
    Chromosome == "NC_055966.1" ~ "10",
    Chromosome == "NC_055967.1" ~ "11",
    Chromosome == "NC_055968.1" ~ "12",
    Chromosome == "NC_055969.1" ~ "13",
    Chromosome == "NC_055970.1" ~ "14",
    Chromosome == "NC_055971.1" ~ "15",
    Chromosome == "NC_055972.1" ~ "16",
    Chromosome == "NC_055973.1" ~ "17",
    Chromosome == "NC_055974.1" ~ "18",
    Chromosome == "NC_055975.1" ~ "19",
    Chromosome == "NC_055976.1" ~ "20",
    Chromosome == "NC_055977.1" ~ "21",
    Chromosome == "NC_055978.1" ~ "22",
    Chromosome == "NC_055979.1" ~ "23",
    Chromosome == "NC_055980.1" ~ "24",
    TRUE ~ Chromosome 
  ))
View(snpschromo)
#this looks good
ggplot(snpschromo, aes(x = Chromosome, y = Variants)) +
  geom_bar(stat = "identity", fill = "navy") +
  theme_classic() +
  labs(
    title = "Variants along Chromosome",
    x = "Chromosome",
    y = "Variants") + 
  theme(axis.text.x = element_text(angle = 45)) +
  theme(plot.title = element_text(hjust=0.5))















