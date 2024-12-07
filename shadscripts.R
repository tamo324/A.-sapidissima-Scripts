setwd("~/Desktop/shads/shadsprojectdata")
library(tidyverse)
library(readxl)

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
#this looks good (laurens ggplot)
 
ggplot(snpschromo, aes(x = Chromosome, y = Variants)) +  
  geom_bar(stat = "identity", fill = "navy", color = "black") +  # Add black outline 
  theme_classic() + 
  scale_x_discrete(limits = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "MT")) + 
  labs(title = "SNP Density", x = "Chromosome", y = "Variants") + 
  theme( 
    plot.title = element_text(hjust = 0.5, vjust = -3, size = 26),  # Lower title position and increase font size 
    axis.title.x = element_text(size = 20),  # Increase x-axis title font size 
    axis.title.y = element_text(size = 20),  # Increase y-axis title font size 
    axis.text.x = element_text(size = 15),   # Increase x-axis text font size 
    axis.text.y = element_text(size = 15)    # Increase y-axis text font size 
  ) 

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

#Tara insert data 
library(tidyverse)
library(vcfR)
#initial work
#This code involves initial visualization of SNPeff gene data
#Import SNPeff gene dataset
#import NCBI gene dataset that provides location and additional info on genes from https://www.ncbi.nlm.nih.gov/datasets/gene/GCF_018492685.1/
#helps add context to SNPeff gene data
#sample ncbi dataset for unique genes, no transcript
ncbisamles <- ncbi_dataset_1_ %>%
  group_by(GeneId) %>%
  sample_n(size=1)

#randomly sample one transcript per gene
genesamples <- snpEff_genes %>%
  group_by(GeneId) %>%
  sample_n(size=1)

#read in SNPeff annotated VCF
VCFann<- read.vcfR("C:/Users/Tara/Downloads/shad.filtered.ann.vcf.gz")
#fix has the info we want(annotations)
fix <- as.data.frame(VCFann@fix)
#narrow data frame to specifically what we need
dat <- fix %>% select("CHROM", "POS", "INFO")
#remove anything before annotations in the INFO column
dat$INFO <- str_replace(dat$INFO, '(.*?)ANN(.*?)', '')
#separate INFO column by specific delimiter used in that column
separate <- dat %>% separate_wider_delim(INFO, delim= "|", 
                                         names= c("allele", "anno", "impact", "gene", "geneid", "ftype", "transcriptid", "biotype", "rank", "hgvsc", "hgvsp", "cdna", "cds", "propos", "distofeat", "errors"), too_many = "merge", too_few= "debug")
#select what we need from the new columns
separate <- separate %>% select("geneid", "CHROM", "POS")
#make data frame we are working with smaller by selecting unique genes
  attempt <- separate %>%
  group_by(geneid) %>%
  sample_n(size=1)
#some genes denoted as start-XXXXX, get rid of anything before XXXXX
attempt$geneid <- str_replace(attempt$geneid, '(.*?)-(.*?)', '')
#this made more redundant gene names, sample again
final <- attempt %>% group_by(geneid) %>%
  sample_n(size=1)
#rename column so easier to join to gene sampe dataset
colnames(final)[7]<- "GeneId"
#join datatset
maybe1<- left_join(genesamples, final, by= "GeneId")
#save as csv
write.csv(maybe1, "maybe1.csv")

library(ggrepel)
#limit to top 50
highest<- maybe1 %>% arrange(desc(variants_impact_HIGH))
#renaming
highest$CHROM[highest$CHROM== "NC_055957.1"] <-"1"
highest$CHROM[highest$CHROM== "NC_055958.1"] <-"2"
highest$CHROM[highest$CHROM== "NC_055959.1"] <-"3"
highest$CHROM[highest$CHROM== "NC_055960.1"] <-"4"
highest$CHROM[highest$CHROM== "NC_055961.1"] <-"5"
highest$CHROM[highest$CHROM== "NC_055962.1"] <-"6"
highest$CHROM[highest$CHROM== "NC_055963.1"] <-"7"
highest$CHROM[highest$CHROM== "NC_055964.1"] <-"8"
highest$CHROM[highest$CHROM== "NC_055965.1"] <-"9"
highest$CHROM[highest$CHROM== "NC_055966.1"] <-"10"
highest$CHROM[highest$CHROM== "NC_055967.1"] <-"11"
highest$CHROM[highest$CHROM== "NC_055968.1"] <-"12"
highest$CHROM[highest$CHROM== "NC_055969.1"] <-"13"
highest$CHROM[highest$CHROM== "NC_055970.1"] <-"14"
highest$CHROM[highest$CHROM== "NC_055971.1"] <-"15"
highest$CHROM[highest$CHROM== "NC_055972.1"] <-"16"
highest$CHROM[highest$CHROM== "NC_055973.1"] <-"17"
highest$CHROM[highest$CHROM== "NC_055974.1"] <-"18"
highest$CHROM[highest$CHROM== "NC_055975.1"] <-"19"
highest$CHROM[highest$CHROM== "NC_055976.1"] <-"20"
highest$CHROM[highest$CHROM== "NC_055977.1"] <-"21"
highest$CHROM[highest$CHROM== "NC_055978.1"] <-"22"
highest$CHROM[highest$CHROM== "NC_055979.1"] <-"23"
highest$CHROM[highest$CHROM== "NC_055980.1"] <-"24"
highest$CHROM[highest$CHROM== "NC_014690.1"] <-"MT"
highest$CHROM[is.na(highest$CHROM)] <-"Unplaced"
highest$CHROM[highest$CHROM== "NW_024582131.1"] <- "Unplaced"
highest$CHROM[highest$CHROM== "NW_024582138.1"]<- "Unplaced"
highest$CHROM[highest$CHROM== "NW_024582116.1"] <- "Unplaced"
highest$CHROM[highest$CHROM== "NW_024582117.1"] <- "Unplaced"
highest$CHROM[highest$CHROM== "NW_024582121.1"] <-"Unplaced"
highest$CHROM[highest$CHROM== "NW_024582123.1"]<- "Unplaced"
highest$CHROM[highest$CHROM== "NW_024582125.1"]<- "Unplaced"
highest$CHROM[highest$CHROM== "NW_024582127.1"]<- "Unplaced"
highest$CHROM[highest$CHROM== "NW_024582130.1"]<- "Unplaced"
highest$CHROM[highest$CHROM== "NW_024582134.1"]<- "Unplaced"
highest$CHROM[highest$CHROM== "NW_024582136.1"]<- "Unplaced"
highest$CHROM[highest$CHROM== "NW_024582139.1"]<- "Unplaced"
highest$CHROM[highest$CHROM== "NW_024582142.1"]<- "Unplaced"
highest$CHROM[highest$CHROM== "NW_024582143.1"]<- "Unplaced"
highest$CHROM[highest$CHROM== "NW_024582146.1"]<- "Unplaced"
highest$CHROM[highest$CHROM== "NW_024582147.1"]<- "Unplaced"
highest$CHROM[highest$CHROM== "NW_024582151.1"]<- "Unplaced"
highest$CHROM[highest$CHROM== "NW_024582153.1"]<- "Unplaced"
highest$CHROM[highest$CHROM== "NW_024582154.1"]<- "Unplaced"
highest$CHROM[highest$CHROM== "NW_024582155.1"]<- "Unplaced"
highest$CHROM[highest$CHROM== "NW_024582157.1"]<- "Unplaced"
#check NA unplaced values
check<- maybe1%>% filter(is.na(CHROM))
#import joined dataset from NCBI(see previous code)
y<- ncbisamples %>% select(Chromosome, GeneId)
#join both datasets
joined<- left_join(highest, y, by= "GeneId")
#replace NA values with any that are not NA in the other column
joined<- joined %>% mutate(CHROM = coalesce(Chromosome))
#check NA values and genes
#still NAs, but randomly chosen to check on NCBI. Genes were either from unplaced scaffold or otherwise unlocated
check1<- joined%>% filter(is.na(CHROM))
highest<- joined
top50<- highest %>% slice(1:50)
#make graph of top 50 most high impact variants
ggplot(top50, aes(x=GeneId, y=variants_impact_HIGH)) +
  geom_col(position= 'identity')+
  labs(x= "Gene Name", y= "Total Variants") +
  ggtitle( "50 Genes with Highest Average High-Impact Variants")+
  theme(axis.text.x = element_text(angle = 90, size= 9, vjust= .5, hjust= 0.95))
#graph top 50 by Chromosome
ggplot(top50, aes(x= factor(CHROM,levels =c("1", "2","4", "6", "7","8", "9", "10", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "MT", "Unplaced")), y=variants_impact_HIGH, label=GeneId)) +
  geom_point(size= 1)+
  labs(x= "Chromosome", y= "Number of High Impact Variants") +
  ggtitle( "50 Genes with Most High-Impact Variants") +
  geom_text_repel(aes(label = GeneId),
                  size= 1.5,
                  hjust= 0.0,
                  vjust= 0.0,
                  point.padding= 0,
                  force=0.005,
                  force_pull=5,
                  color= "black")+
  scale_y_continuous(breaks= seq(5, 25, by= 1)) +
  scale_color_discrete(breaks=c("1", "2","4", "6", "7","8", "9", "10", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "MT", "Unplaced"))

#graph all variants by chromosome, label those over a certain y value
ggplot(highest) +
  aes(x = factor(CHROM,levels =c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "MT", "Unplaced")), y = variants_impact_HIGH) +
  geom_jitter() +
  labs(x= "Chromosome", y= "Number of High Impact Variants") +
  ggtitle( "Number of High Impact Variants by Gene") +
  geom_label(aes(label = GeneId), data= subset(highest, variants_impact_HIGH > 5),
             size= 2,
             hjust= 0.0,
             vjust= -0.5)+
  scale_color_discrete(breaks=c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "MT", "Unplaced"))
#calculated mean 
mean(highest$variants_impact_HIGH)
#make new column with sum of moderate and high impact variants
highest <- highest%>% 
  mutate(total= rowSums(across(c(variants_impact_HIGH, variants_impact_MODERATE, variants_impact_LOW, variants_impact_MODIFIER))))
highest <- highest%>% 
  mutate(tothighmodvariants= rowSums(across(c(variants_impact_HIGH, variants_impact_MODERATE))))
#graph by chromsome and repeat with other data
ggplot(highest) +
  aes(x = factor(CHROM,levels =c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "MT", "Unplaced")), y = tothighmodvariants) +
  geom_jitter() +
  labs(x= "Chromosome", y= "Number of High and Mod Impact Variants") +
  ggtitle( "Number of High and Mod Impact Variants by Gene") +
  geom_label(aes(label = GeneId), data= subset(highest, tothighmodvariants > 20),
             size= 2,
             hjust= 0.0,
             vjust= -0.5)+
  scale_color_discrete(breaks=c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "MT", "Unplaced"))
ggplot(highest) +
  aes(x = factor(CHROM,levels =c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "MT", "Unplaced")), y = total) +
  geom_jitter() +
  labs(x= "Chromosome", y= "Number of Variants") +
  ggtitle( "Number of High, Moderate, Low, and Modifier Impact Variants by Gene") +
  geom_label(aes(label = GeneId), data= subset(highest, total > 1800),
             size= 2,
             hjust= 0.0,
             vjust= -0.5)+
  scale_color_discrete(breaks=c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "MT", "Unplaced"))
#repeat with new data
highest<- highest %>% arrange(desc(tothighmodvariants))
top50<- highest %>% slice(1:50)
ggplot(top50, aes(x=GeneId, y=tothighmodvariants)) +
  geom_col(position= 'identity')+
  labs(x= "Gene Name", y= "Total Variants") +
  ggtitle( "50 Genes with the Most High and Moderate Impact Variants")+
  theme(axis.text.x = element_text(angle = 90, size= 9, vjust= .5, hjust= 0.95))+
  coord_flip()
#repeat with new total
highest<- highest %>% arrange(desc(total))
top50<- highest %>% slice(1:50)
ggplot(top50, aes(x=GeneId, y=total)) +
  geom_col(position= 'identity')+
  labs(x= "Gene Name", y= "Total Variants") +
  ggtitle( "50 Genes with the Most Total Impact Variants")+
  theme(axis.text.x= element_text(angle = 90, size= 9, vjust= .5, hjust= 0.95))+
  coord_flip()




