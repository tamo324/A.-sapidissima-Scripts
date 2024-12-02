setwd("~/Desktop/shads/shadsprojectdata")
library(ape)
library(vcfR)
library(tidyverse)
library(readxl)
library(RcppRoll)

bobvcf2 <- read.vcfR("shad.filtered.ann.vcf.gz")

#chromR object 
bobchromr <- create.chromR(bobvcf2, name = "CHROM", seq = NULL, ann = NULL, verbose = TRUE)


#plot 1 first look
plot(bobchromr) 

#plot 2 a few filters 
chromfilteredbob <- masker(bobchromr, min_QUAL = 1, max_DP = 50, min_MQ = 40,  max_MQ = 61)
plot(chromfilteredbob)

#plot 3 include variants 
chromvariantsbob <- proc.chromR(chromfilteredbob, verbose=TRUE)
plot(chromvariantsbob)

#plot4 
chromoqc(chromvariantsbob, dp.alpha = 66) 
plot(chromoqc)



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
snpschromo$Chromosome <-factor(snpschromo$Chromosome, 
          levels = as.character(sort(as.numeric
          (unique(snpschromo$Chromosome)
          ))))

str(snpschromo)

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
 ### SNP DENSITY GOOD TO GO

snpdens <- read_tsv("snpdens.tsv")
str(snpdens)
unique(snpdens$CHROM)
snpdens$CHROM <- as.factor(snpdens$CHROM)
any(is.na(snpdens$CHROM))
snpdens <- as.data.frame(snpdens)

snpdenchr <- snpdens %>%
  filter(!str_starts(CHROM, "NW")) %>%
  mutate(CHROM = case_when(
    CHROM == "NC_014690.1" ~ "MT",
    CHROM == "NC_055957.1" ~ "01",
    CHROM == "NC_055958.1" ~ "02",
    CHROM == "NC_055959.1" ~ "03",
    CHROM == "NC_055960.1" ~ "04",
    CHROM == "NC_055961.1" ~ "05",
    CHROM == "NC_055962.1" ~ "06",
    CHROM == "NC_055963.1" ~ "07",
    CHROM == "NC_055964.1" ~ "08",
    CHROM == "NC_055965.1" ~ "09",
    CHROM == "NC_055966.1" ~ "10",
    CHROM == "NC_055967.1" ~ "11",
    CHROM == "NC_055968.1" ~ "12",
    CHROM == "NC_055969.1" ~ "13",
    CHROM == "NC_055970.1" ~ "14",
    CHROM == "NC_055971.1" ~ "15",
    CHROM == "NC_055972.1" ~ "16",
    CHROM == "NC_055973.1" ~ "17",
    CHROM == "NC_055974.1" ~ "18",
    CHROM == "NC_055975.1" ~ "19",
    CHROM == "NC_055976.1" ~ "20",
    CHROM == "NC_055977.1" ~ "21",
    CHROM == "NC_055978.1" ~ "22",
    CHROM == "NC_055979.1" ~ "23",
    CHROM == "NC_055980.1" ~ "24",
    TRUE ~ CHROM 
  )) %>%  filter(CHROM != "MT")
snpplot <- ggplot(snpdenchr, aes(x = BIN_START, y = `VARIANTS/KB`, color = CHROM)) +
  geom_col() +
  scale_color_viridis_d(name = "CHROM") + 
  labs(
      x = "Position", 
      y = "Variant Density",
      title = "Variant Density Across Genome") +
  theme_bw () +
  facet_wrap(~CHROM, scales = "free_x")

print(snpplot)

#Calculate mean from SNP density column, take take lowest 2.5 and highest 97.5 % 

rowMeans(snpdenchr)
str(snpdenchr)
meansnps<- mean(snpdenchr$SNP_COUNT, na.rm = TRUE)
print(meansnps) #81.92581 
outliers <- quantile(snpdenchr$SNP_COUNT, probs = c(0.025, 0.975))
print(outliers) #2.5% 6, 97.5% 162
top_97_5 <- snpdenchr[snpdenchr$SNP_COUNT > outliers[2], ]
top <- ggplot(top_97_5, aes(x = BIN_START, y = `VARIANTS/KB`, color = CHROM)) +
         geom_col() +
         scale_color_viridis_d(name = "CHROM") + 
         labs(
           x = "Position", 
           y = "Variant Density",
           title = "Top 97.5%") +
         theme_bw () +
         facet_wrap(~CHROM, scales = "free_x")

bot_2.5 <- snpdenchr[snpdenchr$SNP_COUNT < outliers[1], ]
bot <- ggplot(bot_2.5, aes(x = BIN_START, y = `VARIANTS/KB`, color = CHROM)) +
  geom_col() +
  scale_color_viridis_d(name = "CHROM") + 
  labs(
    x = "Position", 
    y = "Variant Density",
    title = "Bottom 2.5%") +
  theme_bw () +
  facet_wrap(~CHROM, scales = "free_x")

#export top 97 and bottom 2.5 as csv files
write.csv(top_97_5, "top_97_5.csv", row.names = FALSE)
write.csv(bot_2.5, "bot_2_5.csv", row.names = FALSE)




