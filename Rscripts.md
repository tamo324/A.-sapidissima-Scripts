**R code for data anylses** ##update w Tara's code 

R script for Variant Density Plot:
```
#SNP density by rolling window code

library(tidyverse)
#read in data
snpdens <- read_tsv("Asap001.tsv.snpden") #tsv from vcftools of 20000 bp windows
#check 
str(snpdens) #chromosomes "NC_055957.1" etc, bins look good, SNP counts per bin, variants/kb
view(snpdens)
unique(snpdens$CHROM) #nws are unplaced scaffolds
slice_max(snpdens, order_by = `VARIANTS/KB`) #12.1
slice_min(snpdens, order_by = `VARIANTS/KB`) #0


#make chromosomes numbered and take out unplaced scaffolds
snpdenchr <- snpdens %>% filter(!str_starts(CHROM, "NW")) %>% mutate(CHROM = case_when(
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
    TRUE ~ CHROM))

str(snpdenchr) #tibble looks good with the chromosomes numbered now.

#plot across genome in bins, y axis should be variants/kb
snpplot <- ggplot(snpdenchr, aes(x = BIN_START, y = `VARIANTS/KB`, color = CHROM)) +
  geom_line() +
  scale_color_viridis_d(name = "CHROM") + 
  labs(
    x = "Position", 
    y = "Variant Density",
    title = "Variant Density Across Genome") +
  theme_bw () +
  facet_wrap(~CHROM, scales = "free_x") + 
  theme(axis.text.x = element_text(size = 4))

print(snpplot)
```

#Tara insert data 

#read in snpEff_genes.txt
snpEff_genes <- read_delim("C:/Users/XXX/Downloads/snpEff_genes.txt", 
delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1)
genesamples <- snpEff_genes %>%
  group_by(GeneId) %>%
  sample_n(size=1)

#read in SNPeff annotated VCF
VCFann<- read.vcfR("C:/Users/XXX/Downloads/Asap001.filtered.ann.noMito.vcf/Asap001.filtered.ann.noMito.vcf")
#fix within vcf has the info needed(annotations)
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
#rename column so easier to join to gene sample dataset
colnames(final)[1]<- "GeneId"
#join datatset
maybe1<- left_join(genesamples, final, by= "GeneId")
# arrange by number of high impact variants
highest<- maybe1 %>% arrange(desc(variants_impact_HIGH))

#renaming to chromosome
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
highest$CHROM[highest$CHROM== "NC_014690.1"] <-"Unplaced"
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
highest$CHROM[highest$CHROM== "NW_024582157.1"]<- "Unplaced"
highest$CHROM[grepl("unassigned", highest$TranscriptId)] <- "Unplaced"

# define chromosome levels
chrom_levels <- c(as.character(1:24), "Unplaced")

# apply factor levels
highest$CHROM <- factor(highest$CHROM, levels = chrom_levels)
subset_high <- subset(highest, variants_impact_HIGH > 8)
subset_high$CHROM <- factor(subset_high$CHROM, levels = chrom_levels)
subset_low <- subset(highest, variants_impact_HIGH <= 8)
subset_low$CHROM <- factor(subset_low$CHROM, levels = chrom_levels)

#set seed for jitter
jitter_nudge <- position_jitternudge(
  seed = 143, height = 0.4, y=1.2, x=0,
  direction = "split",
  nudge.from = "jittered")
jitter<- position_jitter(height = 0.4, seed = 143)

# custom axis labels to hide "Unplaced"
chrom_labels <- chrom_levels
chrom_labels[chrom_labels == "Unplaced"] <- ""

# generate graph using grid so numbers and "unplaced" x-axis labels easily readable
p <- ggplot(highest, aes(x = CHROM, y = variants_impact_HIGH, label = GeneId)) +
  geom_jitter(data = subset_high, position = jitter, size = 1.5, shape = 21, fill = "red", color = "black") +
  geom_jitter(data = subset_low, position = jitter, size = 1.5) +
  ggrepel::geom_label_repel(data = subset_high, position = jitter_nudge,
                            size = 3, force_pull = 100, force = 0.7) +
  geom_hline(yintercept = 9, linetype = "dashed", color = "red", linewidth = 0.7)+
  scale_x_discrete(drop = FALSE, labels = chrom_labels) +
  scale_y_continuous(breaks = seq(1, 23, 2)) +
  theme(
    axis.text.x = element_text(size = 12),
    plot.margin = margin(t = 10, r = 10, b = 60, l = 10),
    axis.title.x = element_text(margin=margin(t=40), size= 14),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size= 14)
  ) +
  labs(x = "Chromosome", y = "Number of High Impact Variants")

# draw plot
grid.newpage()
vp <- viewport(layout = grid.layout(1, 1))
pushViewport(vp)
print(p, vp = vp)


grid.text("Unplaced", x = unit(0.965, "npc"), y = unit(0.225, "npc"),
          gp = gpar(fontsize = 10), rot = -90)


#Exported as 920 width x 600 height image



