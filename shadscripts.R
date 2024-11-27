setwd("~/Desktop/shads")
library(ape)
library(vcfR)
library(tidyverse)
library(readxl)


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



#################
#trying to make a chromsomewide plot heterzygosity (nHet) from the tsv file 
#https://cran.r-project.org/web/packages/RcppRoll/RcppRoll.pdf <-rollmean
install.package ("RcppRoll")
het <- read_tsv("heteroz.tsv")
head(het)
colnames(het)
hetfilter <- het %>%
  filter(!str_starts(CHR, "NW"))
#dplyr sliding window for this. plot het (y-axs) along position (x axis)

#trying to make windows and apply 
View(hetfilter)
het2 <- hetfilter %>%
  group_by(POS) %>%
  arrange(CHR) %>% # 
  mutate(roll_mean_nHet = roll_mean(x = nHet, n = 5000, by = 5000, fill = NA, align = "center")
 

# Remove NA values introduced by rolling mean
het2 <- het2 %>%
  filter(!is.na(roll_mean_nHet))

# Scatterplot of rolling means - this is not working/ just looks horroble 
ggplot(het2, aes(x = POS, y = roll_mean_nHet, color = CHR)) +
  geom_line() +
  labs(
    x = "Genomic Position",
    y = "Rolling Mean of Heterozygosity (nHet)",
    title = "Scatterplot of Rolling Mean Heterozygosity"
  ) +
  theme_minimal()

hetwindow <- roll_mean(hetfilter, n = 5000, weights = NULL, by = 5000, fill = numeric(0),
          partial = FALSE, align = c("center"), normalize = TRUE,
          na.rm = FALSE)
#heterozygosity
het <- read_tsv("heteroz.tsv")
View(het)
head(het)
colnames(het)
hetfilter <- het %>%
  filter(!str_starts(CHR, "NW")) %>%
  arrange(CHR, POS) %>%
  mutate(CHR = case_when(
  CHR == "NC_014690.1" ~ "MT",
  CHR == "NC_055957.1" ~ "1",
  CHR == "NC_055958.1" ~ "2",
  CHR == "NC_055959.1" ~ "3",
  CHR == "NC_055960.1" ~ "4",
  CHR == "NC_055961.1" ~ "5",
  CHR == "NC_055962.1" ~ "6",
  CHR == "NC_055963.1" ~ "7",
  CHR == "NC_055964.1" ~ "8",
  CHR == "NC_055965.1" ~ "9",
  CHR == "NC_055966.1" ~ "10",
  CHR == "NC_055967.1" ~ "11",
  CHR == "NC_055968.1" ~ "12",
  CHR == "NC_055969.1" ~ "13",
  CHR == "NC_055970.1" ~ "14",
  CHR == "NC_055971.1" ~ "15",
  CHR == "NC_055972.1" ~ "16",
  CHR == "NC_055973.1" ~ "17",
  CHR == "NC_055974.1" ~ "18",
  CHR == "NC_055975.1" ~ "19",
  CHR == "NC_055976.1" ~ "20",
  CHR == "NC_055977.1" ~ "21",
  CHR == "NC_055978.1" ~ "22",
  CHR == "NC_055979.1" ~ "23",
  CHR == "NC_055980.1" ~ "24",
  TRUE ~ CHR 
))


#filter out the reads not aligned
#dplyr sliding window for this. plot het (y-axs) along position (x axis)

#trying to make windows and apply (was told to use roll_mean_)
#roll_mean(x, n = 1L, weights = NULL, by = 1L, fill = numeric(0),
         # partial = FALSE, align = c("center", "left", "right"), normalize = TRUE,
         # na.rm = FALSE)
View(hetfilter)
het2 <- hetfilter %>%
  group_by(CHR) %>%
  mutate(roll_mean_nHet = roll_mean(x = nHet, n = 3000, fill = NA, align = "center")) %>%
  ungroup()
 

# Scatterplot of rolling means - this is not working/ just looks horroble 
hetplot <- ggplot(het2, aes(x = CHR, y = roll_mean_nHet, color = CHR)) +
  geom_line() +
  scale_color_viridis_d(name = "Chromosome") + 
  labs(
    x = "Genomic Position",
    y = "Rolling Mean of Heterozygosity (nHet)",
    title = "Scatterplot of Rolling Mean Heterozygosity"
  ) + theme_gray()
print(hetplot)















