############## Antibiotic resistance genes detected in lichens: insights from Cladonia stellaris 
############## Marta Alonso-García, Paul B. L George, Samantha Leclerc, Marc Veillette, Caroline Duchaine and Juan Carlos Villarreal

# Set working directory (adjust as needed)
#setwd("C:/Users/Marta/ARGs") 

# Load required packages
library("dada2")      # Processing 16S rRNA amplicon data
library("phyloseq")   # Microbial community analysis
library("vegan")      # Diversity and multivariate analyses
library("Hmisc")      # Correlation analyses
library("igraph")     # Network analyses
# Add other packages as needed



#####
# ------------------------------------------ 1. Analyses of 16S sequences
#####

## Set path to raw 16S FASTQ files (from Alonso-García & Villarreal 2022)
Path16S <- "Sequences_16S_42samples_ARG"  
list.files(Path16S)
 
##------------------------------------------ 1.1 Filter and trim
fnFs <- sort(list.files(Path16S, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(Path16S, pattern="_R2_001.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

## Inspect read quality 
plotQualityProfile(fnFs[1:4])
plotQualityProfile(fnRs[1:4])

## Define output filtered files
filtFs <- file.path(Path16S, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(Path16S, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Filter/trimming: retain high-quality bases only (Q>30), remove primers, trim low-quality ends
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(280,220),
                     maxN=0, maxEE=c(2,2), truncQ=2, trimLeft=c(10,10), rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
head(out)

## Re-check quality after trimming
plotQualityProfile(filtFs[1:4])
plotQualityProfile(filtRs[1:4])


##------------------------------------------ 1.2 Infer amplicon sequence variants (ASVs)
errFs <- learnErrors(filtFs, multithread = TRUE)
errRs <- learnErrors(filtRs, multithread = TRUE)
plotErrors(errFs, nominalQ=TRUE) 
plotErrors(errRs, nominalQ=TRUE)

## Sample inference using pseudo pooling (recommended for sparse microbial datasets)
dadaFsPseudo <- dada(filtFs, err=errFs, multithread = TRUE, pool = "pseudo")
dadaRsPseudo <- dada(filtRs, err=errRs, multithread = TRUE, pool = "pseudo")
dadaFsPseudo[[1]] ## 770 sequence variants were inferred from 30632 input unique sequences.
dadaRsPseudo[[1]] ## 667 sequence variants were inferred from 26132 input unique sequences..

## Merge paired reads and construct sequence table
mergers <- mergePairs(dadaFsPseudo, filtFs, dadaRsPseudo, filtRs, verbose=TRUE)
head(mergers[[1]])
seqtab <- makeSequenceTable(mergers)
dim(seqtab) ## 42 samples and 23923 unique ASVs identified across all the samples

## Inspect distribution of sequence lengths 
table(nchar(getSequences(seqtab))) ## The majority of your sequences are around 420-445 base pairs (bp) long, 

## Filter sequences by expected length (420-447 bp)
seqtab2<- seqtab[,nchar(colnames(seqtab)) %in% seq(420,447)] ## 42 samples and 23883 unique ASVs identified across all the samples

## Track reads through pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFsPseudo, getN), sapply(dadaRsPseudo, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
track
write.table(track,"Step_01_Processed_reads_16S.txt", sep= "\t",
            row.names = TRUE)


##------------------------------------------ 1.3 Assign taxonomy for all samples

set.seed(100) # Random number generator
taxa <- assignTaxonomy(seqtab.nochim, "01_Assign_taxonomy/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE, minBoot = 80)
dim(taxa) # 2013 ASVs and 6 taxonomic levels



#####
##------------------------------------------ 2. Handoff to phyloseq 
#####

## Check the row names of the OTU table (ensure they match the CSV file)
rownames(seqtab.nochim) 

## Load sample metadata CSV
samdf_ARG <- read.csv('./02_Inputs/samdf_ARG.csv', header = TRUE, sep=';')
rownames(samdf_ARG) <- samdf_ARG$File_name  
samdf_ARG

library(phyloseq); packageVersion("phyloseq")
library("DECIPHER")
library("phangorn")

## Build phyloseq object (excluding the phylogenetic tree)
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), sample_data(samdf_ARG),
               tax_table(taxa))

## Check structure
ps #  2013 taxa and 42 samples
head(sample_data(ps), 20)

## Use short names for ASVs instead of full DNA sequence
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

## Remove features with ambiguous phylum, chloroplasts, or mitochondria
ps <- subset_taxa(ps, !is.na(Phylum) & !Order %in% c("Chloroplast") 
                  & !Family %in% c("Mitochondria"))

## Table of features per phylum
table(tax_table(ps)[, "Phylum"], exclude = NULL)

## Remove phyla with fewer than 10 taxa (manual selection based on previous table)
ps <- subset_taxa(ps, !Phylum %in% c("Firmicutes", "Chloroflexi"))
table(tax_table(ps)[, "Phylum"], exclude = NULL) 
ps #  1659 taxa and 42 samples

## All samples have >5000 reads, no need to remove samples



#####
##------------------------------------------ 3. Prevalence filtering
#####

## Load required packages
library(ggplot2); packageVersion("ggplot2")
theme_set(theme_bw())

## Reads per sample
sample_sums(ps)
hist(sample_sums(ps), main="Histogram: Read Counts", 
     xlab="Total Reads", border="black", col="darkgreen", las=1, breaks=15)

## Reads per taxa
head(taxa_sums(ps), 10)

## Compute prevalence of each feature
prevdf = apply(X = otu_table(ps),
               MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
			   
## Add taxonomy and total abundance
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps),
                    tax_table(ps))
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),
                                                  sum(df1$Prevalence))})


## Subset remaining phyla
prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(ps, "Phylum"))
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ps),color=Phylum)) +
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")

## Keep ASVs present in at least 2 samples
prevalenceThreshold = 0.047 * nsamples(ps) ## # 2/42 = 0.047
keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
ps2 = prune_taxa(keepTaxa, ps)

table(tax_table(ps2)[, "Phylum"], exclude = NULL) 
ps2 ## 481 taxa and 42 samples

## Compute prevalence again for filtered dataset
prevdf = apply(X = otu_table(ps2),
               MARGIN = ifelse(taxa_are_rows(ps2), 
                              yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps2),
                    tax_table(ps2))
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),
                                                  sum(df1$Prevalence))})
prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(ps2, "Phylum"))

ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ps2),
                    color=Phylum)) +
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  
  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + 
  ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")

## Export OTU table
write.csv2(otu_table(ps2), file = "Step_03_otu_table_ps2.csv")



#####
##------------------------------------------ 4. Phylogenetic tree
#####

## Load required packages
library(phangorn)
library(Biostrings)
library(DECIPHER)

## Align sequences and calculate distances
dna<- Biostrings::DNAStringSet(refseq(ps2.norm))
names(dna) <- taxa_names(ps2.norm)
aligned_seqs <- AlignSeqs(dna)
distances<- dist.ml(as.phyDat(as.DNAbin(aligned_seqs)))

## Construct Neighbor-Joining tree
NJ_tree <- NJ(distances)

## Maximum likelihood tree
fit = pml(NJ_tree, data=as.phyDat(as.DNAbin(aligned_seqs)))
fitGTR <- optim.pml(fit, model="GTR", optInv=TRUE, optGamma=TRUE, rearrangement = "stochastic", control = pml.control(trace = 0))

## Root the tree
rooted_tree <- midpoint(fitGTR$tree)

## Visualize tree
library("ggtree")
ggtree(rooted_tree) + 
  geom_tiplab(size = 2) + 
  theme_tree2() +
  labs(title = "Phylogenetic Tree of Bacterial Communities in Cladonia stellaris",
       subtitle = "Constructed using GTR+I model with midpoint rooting")

## Merge tree with phyloseq object
ps2.tree <- merge_phyloseq(ps2, rooted_tree)
ps2.tree



#####
##------------------------------------------ 5. Transform to relative abundance (ps2 = données brutes)
#####

## Load required packages
library(ggplot2)

## Convert raw counts to relative abundances to account for differences in library size
ps2.ra <- transform_sample_counts(ps2, function(x){x / sum(x)})
all(round(sample_sums(ps2.ra), 10) == 1) # verify sums = 1

## Check OTU table: 481 ASVs across 42 lichen samples
otu_table(ps2.ra) ## Values represent the relative contribution of each ASV to the sample's total sequences


## ------------------------------------------ 5.1 Phylum-level relative abundance (Figure S1)
theme_set(theme_bw())

## Exclude non-bacterial taxa (already filtered in ps2)
plot_bar(ps2.ra, x="Voucher", fill="Phylum") + 
  geom_bar(aes(color = Phylum, fill = Phylum), stat="identity", position="stack") +
  labs(
    x = "Samples", 
    y = "Relative Abundance\n",
    title = "Relative abundance of bacterial phyla in Cladonia stellaris\nfrom northern and southern lichen woodlands") +
  theme(panel.background = element_blank()) +
  facet_wrap(~LW, scales="free_x",nrow=1)

## Save figure
ggsave("FigureS1A_RelativeAbundancePhyla.svg", width = 12, height = 6)Figure: Step_05_1_RelativeAbundance.svg

## Genus-level relative abundance (Figure S1B)
## Remove unassigned genera
ps2.ra.genus <- subset_taxa(ps2.ra, !is.na(Genus) & Genus != "")

plot_bar(ps2.ra.genus, x="Voucher", fill="Genus") + 
  geom_bar(aes(color = Genus, fill = Genus), stat="identity", position="stack") +
  labs(
    x = "Samples", 
    y = "Relative Abundance\n",
    title = "Relative abundance of assigned genera in Cladonia stellaris\nfrom northern and southern Lichen woodlands") +
  theme(panel.background = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  facet_wrap(~LW, scales="free_x", nrow=1)

## Save figure
ggsave("FigureS1B_RelativeAbundanceAssignedGenus.svg", width = 12, height = 8)



#####
##------------------------------------------ 6. Normalize ASV counts using DESeq2
#####

## Load required packages
library(DESeq2)
library(phyloseq)

## Convert phyloseq object to DESeq2 dataset and estimate size factors using poscounts
dds <- phyloseq_to_deseq2(ps2.tree, ~ 1)
dds <- estimateSizeFactors(dds, type="poscounts")

## Get normalized counts and update phyloseq object
normalized_counts <- counts(dds, normalized=TRUE)
normalized_counts_t <- t(normalized_counts)

## Create a new phyloseq object with DESeq2-normalized counts
ps2.deseq <- ps2.tree
otu_table(ps2.deseq) <- otu_table(normalized_counts_t, taxa_are_rows = FALSE)



#####
##------------------------------------------ 7. Subset samples by lichen woodland 
##### 

## Load required packages
library(phyloseq)

## ps2 subsets
ps2.NLW <- subset_samples(ps2, LW == "Kuujuarapik")
ps2.NLW <- prune_taxa(taxa_sums(ps2.NLW) > 0, ps2.NLW)
table(tax_table(ps2.NLW)[, "Phylum"], exclude = NULL)

ps2.SLW <- subset_samples(ps2, LW == "PNGJ")
ps2.SLW <- prune_taxa(taxa_sums(ps2.SLW) > 0, ps2.SLW)
table(tax_table(ps2.SLW)[, "Phylum"], exclude = NULL)

## ps2.deseq subsets
ps2.deseq.NLW <- subset_samples(ps2.deseq, LW == "Kuujuarapik")
ps2.deseq.NLW <- prune_taxa(taxa_sums(ps2.deseq.NLW) > 0, ps2.deseq.NLW)
table(tax_table(ps2.deseq.NLW)[, "Phylum"], exclude = NULL)

ps2.deseq.SLW <- subset_samples(ps2.deseq, LW == "PNGJ")
ps2.deseq.SLW <- prune_taxa(taxa_sums(ps2.deseq.SLW) > 0, ps2.deseq.SLW)
table(tax_table(ps2.deseq.SLW)[, "Phylum"], exclude = NULL)



#####
##------------------------------------------ 8. Alpha diversity in Cladonia stellaris 
#####

library(phyloseq)
library(vegan)
library(car)
library(ggplot2)

## Calculate Shannon diversity using DESeq2-normalized counts
alpha_div <- data.frame(
  Shannon = diversity(otu_table(ps2.deseq), index = "shannon")
)
alpha_div <- cbind(sample_data(ps2.deseq), alpha_div)


##------------------------------------------ 8.1 Statistical comparison between LWs  

## Shapiro-Wilk test for normality
shapiro_NLW <- shapiro.test(alpha_div$Shannon[alpha_div$LW == "Kuujuarapik"])
print(shapiro_NLW)
shapiro_SLW <- shapiro.test(alpha_div$Shannon[alpha_div$LW == "PNGJ"])
print(shapiro_SLW)

## Levene's test for homogeneity of variances
leveneTest(Shannon ~ LW, data = alpha_div) 

## Unpaired t-test between LWs
t.test(Shannon ~ LW, data = alpha_div, paired=FALSE, alternative="two.sided")


##------------------------------------------ 8.2 Plot Shannon alpha diversity (Figura 1A)
ggplot(alpha_div, aes(x = LW, y = Shannon, fill = LW)) +
  geom_boxplot(outlier.shape = NA, color = "black") +  # Ajoute une bordure noire
  geom_jitter(width = 0.2, alpha = 0.5) +
  labs(title = "Bacteria alpha diversity in Cladonia stellaris",
       x = "Lichen woodlands",
       y = "Alpha diversity") +
  theme_minimal() +
  theme(
    #panel.grid = element_blank(),  # Supprime toutes les lignes de grille
    panel.grid.major = element_line(color = "grey90", size = 1.2),  # Grille majeure grise
    panel.grid.minor = element_line(color = "grey90", size = 0.8),  # Grille mineure plus claire
    axis.line = element_line(color = "black"),  # Ajoute des lignes d'axe noires
    legend.position = "none",  # Supprime la légende si elle n'est pas nécessaire
    panel.border = element_rect(color = "black", fill = NA, size = 0.5)  # Ajoute un cadre noir autour du graphique
  )

## Save figure
ggsave("Figure1A_ShannonAlphaDiversity.svg", plot = p_alpha, width = 6, height = 5)



#####
##------------------------------------------ 9. Beta diversity comparison between northern and southern lichen woodlands
#####

library("vegan")
library("phyloseq")
library("ape")
library("ggplot2")

##------------------------------------------ 9.1. Bray-Curtis distance
BrayC_dist <- phyloseq::distance(ps2.ra, method="bray")

## Homogeneity of dispersions
sample_data_df <- data.frame(sample_data(ps2.ra))
betadisperBC  <- betadisper(BrayC_dist, sample_data_df$LW)
anova(betadisperBC) # Check if dispersions are homogeneous

## PERMANOVA
permanovaBC <- adonis2(BrayC_dist ~ LW, data = sample_data_df, permutations = 999)
print(permanovaBC) # Test for differences in bacterial community composition


##------------------------------------------ 9.2. Ordination - PCoA
set.seed(420)

pcoa.LW.BC <- ordinate(ps2.ra, method = "PCoA", distance = "bray", normalized=FALSE)
evals <- pcoa.LW.BC$values[,1]
plot_ordination(ps2.ra, pcoa.LW.BC, color = "LW") + # label = "Collection_number"
  labs(col = "Lichen woodland")+
  coord_fixed(sqrt(evals[2] / evals[1]))+
  theme_bw() + theme(text = element_text(size=12))+
  geom_point(size = 3)

## Save figure for manuscript (Figure 1B)
ggsave("Figure1B_PCoA_BrayCurtis_LW.svg", width = 8, height = 6)



#####
##------------------------------------------ 10. Variance explained by LW
#####

library(vegan)

##  Prepare explanatory variable (LW)
LW_df <- data.frame(LW = sample_data_df$LW)
LW_df$LW <- as.factor(LW_df$LW)

## db-RDA: effect of LW on bacterial community composition
dbRDA_16S_LW <- capscale(BrayC_dist ~ LW, data = LW)
anova(dbRDA_16S_LW) # Test global significance
anova(dbRDA_16S, by="terms") # Significance of LW
summary(dbRDA_16S_LW)  # Check proportion of variance explained by LW



##### 
##------------------------------------------ 11. Visualiser les ARGs présentes dans Cladonia stellaris
#####

library(tidyr)
library(dplyr)
library(ggplot2)

## Load ARG dataset
samdf_ARG_heatmap <- read.csv("02_Inputs/samdf_ARG_heatmap.csv", sep = ";", header = TRUE)

## Data cleaning
samdf_ARG_heatmap <- samdf_ARG_heatmap %>%
  mutate(across(c(Latitude, Longitude, Sul1, Sul2, blaCTXM1, blaMOXCMY, blaTEM, blaVIM, aac6Ib, aac3, qnrB, qepA, int1amarko), 
                ~as.numeric(gsub(",", ".", .))))%>%
  mutate(across(c(Sul1, Sul2, blaCTXM1, blaMOXCMY, blaTEM, blaVIM, aac6Ib, aac3, qnrB, qepA, int1amarko), 
                ~replace_na(., 0)))

## Convert to long format for plotting
samdf_ARG_heatmap_long <- samdf_ARG_heatmap %>%
  pivot_longer(cols = c(Sul1, Sul2, blaCTXM1, blaMOXCMY, blaTEM, blaVIM,
                        aac6Ib, aac3, qnrB, qepA, int1amarko),
               names_to = "ARG",
               values_to = "Abundance") %>%
  mutate(Log_Abundance = ifelse(Abundance > 0, log10(Abundance), NA))

## Heatmap of ARG abundance
p <- ggplot(samdf_ARG_heatmap_long, aes(x = ARG, y = Voucher, fill = Log_Abundance)) +
  geom_tile(color = "grey80", size = 0.3) +  # Ligne grise autour de chaque tuile
  scale_fill_viridis_c(option = "plasma", 
                       na.value = "white", 
                       breaks = log10(c(0.0001, 0.001, 0.01, 0.1, 1)),
                       labels = c("0.0001", "0.001", "0.01", "0.1", "1"),
                       limits = c(min(samdf_ARG_heatmap_long$Log_Abundance[is.finite(samdf_ARG_heatmap_long$Log_Abundance)]), 
                                  max(samdf_ARG_heatmap_long$Log_Abundance, na.rm = TRUE))) +
  facet_grid(LW ~ ., scales = "free_y", space = "free_y") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(title = "Relative abundance of ARGs present in Cladonia stellaris",
       x = "Antibiotic resistance genes (ARGs) et mobile genetic elements (MGE)",
       y = "Samples of Cladonia stellaris",
       fill = "Relative abundance")

## Save figure for manuscript (Figure 2)
ggsave("Figure2_ARGs_Heatmap.svg", plot = p, width = 8, height = 6, units = "in")



#####
##------------------------------------------ 12. Comparison of ARG relative abundances between LWs
#####

library(tidyr)
library(dplyr)
library(ggplot2)

##------------------------------------------ 12.1 Load and prepare ARG dataset
samdf_11ARG <- read.csv("02_Inputs/samdf_ARG_LW_clades.csv", sep=";", header=TRUE)

samdf_11ARG <- samdf_11ARG %>%
  pivot_longer(cols = c(Sul1, Sul2, blaCTXM1, blaMOXCMY, blaTEM, blaVIM, aac6Ib, aac3, qnrB, qepA, int1amarko),
               names_to = "ARG",
               values_to = "Abundance")

samdf_11ARG$Abundance[is.na(samdf_11ARG$Abundance)] <- 0

samdf_11ARG <- samdf_11ARG %>%
  mutate(Abundance = as.numeric(gsub(",", ".", Abundance)))
print(n=25, samdf_11ARG)


##------------------------------------------ 12.2 Identify prevalent ARGs (> 50% of samples)

prevalence_11ARG <- samdf_11ARG %>%
  group_by(ARG) %>% 
  summarise(prevalence = mean(Abundance > 0))

samdf_3ARG <- prevalence_11ARG %>%
  filter(prevalence > 0.5) %>%
  pull(ARG)
  
samdf_3ARG <- samdf_11ARG %>%
  filter(ARG %in% samdf_3ARG)


##------------------------------------------ 12.3 Test normality (Shapiro-Wilk) -> distributions non-normal
normality_results <- samdf_3ARG %>%
  group_by(ARG) %>%
  summarise(
    p_value = ifelse(length(unique(Abundance)) > 2,
                     shapiro.test(Abundance)$p.value, NA)
  )

print(normality_results)

##------------------------------------------ 12.4 Wilcoxon rank-sum test for prevalent ARGs (filtering zeros)
library(dplyr)

samdf_3ARG_filtered <- samdf_3ARG %>%
  filter(Abundance > 0)

apply_wilcox_test_filtered <- function(data) {
  wilcox.test(Abundance ~ LW, data = data, exact = FALSE)$p.value
}

wilcox_results_filtered <- samdf_3ARG_filtered %>%
  group_by(ARG) %>%
  summarise(p_value = apply_wilcox_test_filtered(cur_data_all()))

wilcox_results_filtered <- wilcox_results_filtered %>%
  mutate(adjusted_p = p.adjust(p_value, method = "BH"))
  
print(wilcox_results_filtered) ## qnrB shows significant difference between LWs (adjusted p < 0.05)


##------------------------------------------ 12.5 Visualize relative abundance of prevalent ARGs (Figure 3)
library(ggplot2)

ggplot(samdf_3ARG_filtered, aes(x = LW, y = Abundance, fill = LW)) +
  geom_boxplot() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") + # Ajout de la ligne de référence
  facet_wrap(~ ARG, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Relative abundance of the three prevalent ARGs by Lichen woodland",
       x = "Lichen woodland",
       y = "Relative abundance")

## Save figure for manuscript (Figure 3)
ggsave("Figure3_ARGs_Prevalent.svg", plot = p3, width = 8, height = 6, units = "in")


##------------------------------------------ 12.6 Total relative abundance of all 11 ARGs per sample
total_RelAbu_11ARG <- samdf_11ARG %>%
  group_by(Voucher, LW) %>%
  summarise(total_RelAbu_11ARG = sum(Abundance, na.rm = TRUE))

## Test normality (Shapiro-Wilk)
shapiro_total_RelAbu_11ARG <- total_RelAbu_11ARG %>%
  ungroup() %>%
  summarise(p_value = shapiro.test(total_RelAbu_11ARG)$p.value)
print(shapiro_total_RelAbu_11ARG)

## Wilcoxon rank-sum test for total ARG abundance
wilcox_results_RelAbu_11ARG <- wilcox.test(total_RelAbu_11ARG ~ LW, data = total_RelAbu_11ARG, exact = FALSE)
print(wilcox_results_RelAbu_11ARG) # No significant difference between LWs (p > 0.05)


##------------------------------------------ 12.7 Visualize total relative abundance of 11 ARGs (Figure S2)
library(ggplot2)

ggplot(total_RelAbu_11ARG, aes(x = LW, y = total_RelAbu_11ARG, fill = LW)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +  # Ajoute des points pour chaque échantillon
  #geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Total relative abundance of 11 ARGs by Lichen woodland",
       x = "Lichen woodland",
       y = "Total relative abundance") +
  scale_y_log10() +  # Échelle logarithmique pour mieux visualiser les différences
  annotation_logticks(sides = "l")  # Ajoute des marques de graduation logarithmiques

## Save figure for manuscript (Figure S2)
ggsave("FigureS2_Total_ARGs.svg", plot = pS2, width = 8, height = 6, units = "in")



#####
##------------------------------------------ 13. Variance of relative abundance of ARGs explained by LW and Clade
#####

##------------------------------------------ 13.1 Data preparation
# Load required packages
library(vegan)
library(phyloseq)

## Load ARGs relative abundance file
samdf_ARG_RDA <- read.csv("02_Inputs/samdf_ARG_RDA.csv", sep=";", header=TRUE)

## Convert numeric columns
numeric_cols <- c("Latitude", "Longitude", "Sul1", "Sul2", "blaCTXM1", "blaMOXCMY", "blaTEM", "blaVIM", "aac6Ib", "aac3", "qnrB", "qepA", "int1amarko")
for(col in numeric_cols) {
  samdf_ARG_RDA[[col]] <- as.numeric(gsub(",", ".", samdf_ARG_RDA[[col]]))
}
samdf_ARG_RDA[is.na(samdf_ARG_RDA)] <- 0
options(scipen = 999) # avoid scientific notation

## Select only ARGs columns
ARG_abundance <- samdf_ARG_RDA[, c("Sul1", "Sul2", "blaCTXM1", "blaMOXCMY", "blaTEM", "blaVIM", "aac6Ib", "aac3", "qnrB", "qepA", "int1amarko")]
rownames(ARG_abundance) <- samdf_ARG_RDA$File_name

## Compute Bray-Curtis distance matrix
BrayC_dist_ARG <- vegdist(ARG_abundance, method = "bray")
BrayC_dist_ARG <- as.matrix(BrayC_dist_ARG)


##------------------------------------------ 13.2 Variance explained by LW
dbRDA_ARG_LW <- capscale(BrayC_dist_ARG ~ LW, data = LW)
anova(dbRDA_ARG_LW) # Test the global effect of LW
summary(dbRDA_ARG_LW) # Check proportion of variance explained by LW



#####
##------------------------------------------ 14. Prepare matrices to compare ARGs and 16S
#####

## Load required packages
library(vegan)
library(phyloseq)

## ARG relative abundance data frame (already prepared)
ARG_abundance

## Bray-Curtis distance matrix for ARGs (already prepared from ARG_abundance)
BrayC_dist_ARG

## Extract taxonomy table
genusASV <- tax_table[, "Genus"]
familyASV <- tax_table[, "Family"]
orderASV <-tax_table[, "Order"]
classASV <-tax_table[, "Class"]
phylaASV <- tax_table[, "Phylum"]

## Replace missing genus with "Family_Unclassified" or "Unclassified_Genus"
genusModified <- genusASV
for (i in 1:length(genusASV)) {
  if (is.na(genusASV[i])) {
    if (!is.na(familyASV[i])) {
      genusModified[i] <- paste0(familyASV[i], "_Unclassified")
    } else {
      genusModified[i] <- "Unclassified_Genus"
    }
  }
}
genusASV_Fam <- genusModified

## Convert phyloseq OTU table to data frame
otu_table_ps2.ra <- as.data.frame(otu_table(ps2.ra))

## Calculate relative abundances per taxonomic rank
genusAbun <- aggregate(t(otu_table_ps2.ra), by = list(Genus = genusASV_Fam), FUN = sum)
familyAbun <- aggregate(t(otu_table_ps2.ra), by = list(Family = familyASV), FUN = sum)

## Transpose and convert to numeric, keeping taxon names
genusAbun_t <- t(genusAbun)
genusAbun_t <- as.data.frame(genusAbun_t)
genus_names <- genusAbun_t[1,] 
genusAbun_t_num <- data.frame(lapply(genusAbun_t[-1,], as.numeric))
colnames(genusAbun_t_num) <- genus_names  
rownames(genusAbun_t_num) <- rownames(genusAbun_t)[-1]  

familyAbun_t <- t(familyAbun)
familyAbun_t <- as.data.frame(familyAbun_t)
family_names <- familyAbun_t[1,] 
familyAbun_t_num <- data.frame(lapply(familyAbun_t[-1,], as.numeric))
colnames(familyAbun_t_num) <- family_names
rownames(familyAbun_t_num) <- rownames(familyAbun_t)[-1]

# Verification
str(genusAbun_t_num)
str(familyAbun_t_num)

# Prepare OTU table transposed for downstream analyses
otu_table_ps2.ra_t <- t(otu_table_ps2.ra)



#####
##------------------------------------------ 15. Variance of ARG relative abundance explained by bacterial family abundances
#####

## Load required packages
library(vegan)
library(phyloseq)
library(ggplot2)
library(reshape2)

## Response variables: ARG relative abundances
ARG_abundance

## Explanatory variables: bacterial family relative abundances
str(familyAbun_t_num)


##------------------------------------------ 15.1 Multicollinearity check
## Spearman correlation among explanatory variables
cor_matrix_rda <- cor(familyAbun_t_num, method = "spearman") 
write.csv2(cor_matrix_rda, file="cor_matrix_rda.csv")

## Remove highly correlated bacterial families (|ρ| > 0.4)
familyAbun_t_num_reduced <- familyAbun_t_num[, !colnames(familyAbun_t_num) %in% c("Bdellovibrionaceae", "Sphingomonadaceae", 
                                                                                  "Isosphaeraceae", "SM2D12", "Caulobacteraceae", "Microbacteriaceae")]

## Check correlations again
cor_matrix_rda_reduced <- cor(familyAbun_t_num_reduced, method = "spearman")
write.csv2(cor_matrix_rda_reduced, file="cor_matrix_rda_reduced.csv")

## Heatmap of correlations between reduced explanatory variables
cor_melted <- melt(cor_matrix_rda_reduced)
ggplot(cor_melted, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1, 1), space = "Lab", 
                       name="Spearman\nCorrelation") +
  theme_minimal() +
  labs(x = "ARGs", y = "Families") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


##------------------------------------------ 15.2 Linear relationship check
## Pearson correlations between ARGs and reduced bacterial families
cor_matrix_linear <- cor(ARG_abundance, familyAbun_t_num_reduced, method = "pearson")
write.csv2(cor_matrix_linear, file="cor_matrix_linear.csv")

##Heatmap of Pearson correlations
cor_melted <- melt(cor_matrix_linear)
ggplot(cor_melted, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1, 1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal() +
  labs(x = "ARGs", y = "Families") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


##------------------------------------------ 15.3 Canonical Correspondence Analysis (CCA)
## Execute CCA with ARGs as response and reduced bacterial families as explanatory variables
cca_ARGFam <- cca(ARG_abundance ~ ., data = familyAbun_t_num_reduced)
anova(cca_ARGFam, permutations = 999) # Test significance of the CCA model using permutation ANOVA



#####
##------------------------------------------ 16. Data preparation for network analysis: ARGs relative abundance vs. bacterial abundances
#####

##------------------------------------------ 16.1 Absolute abundances transformed with CLR for bacteria grouped by genus and family

library(compositions)
library(phyloseq)


## Relative abundance of ARGs is the recommended metric for correlation analyses
str(ARG3_abundance) # Dataset generated in Step 13.4

## CLR transformation of ASV counts for absolute abundances
otu_table_ps2 <- as.data.frame(otu_table(ps2))
otu_table_ps2_nozero <- otu_table_ps2 + 1  # Ajoute 1 à toutes les valeurs
otu_table_ps2_clr <- as.data.frame(clr(otu_table_ps2_nozero))


### GENUS-level abundances
genusASV_Fam 

## Sum abundances per genus
genusAbun_clr <- aggregate(t(otu_table_ps2_clr), by = list(Genus = genusASV_Fam), FUN = sum)

## Transpose and convert to data.frame
genusAbun_clr_t <- t(genusAbun_clr)
genusAbun_clr_t <- as.data.frame(genusAbun_clr_t)
genus_names <- genusAbun_clr_t[1,] 
genusAbun_clr_t_num <- data.frame(lapply(genusAbun_clr_t[-1,], as.numeric))
colnames(genusAbun_clr_t_num) <- genus_names  
rownames(genusAbun_clr_t_num) <- rownames(genusAbun_clr_t)[-1] 

## Verification
str(genusAbun_clr_t_num)

## Combine ARG and genus CLR datasets
comb_ARG3_genusAbunclr <- cbind(ARG3_abundance, genusAbun_clr_t_num)
dim(comb_ARG3_genusAbunclr)
str(comb_ARG3_genusAbunclr)


### FAMILY-level abundances
familyASV 

## Sum abundances per family
famAbun_clr <- aggregate(t(otu_table_ps2_clr), by = list(Family = familyASV), FUN = sum)

## Transpose and convert to data.frame
famAbun_clr_t <- t(famAbun_clr)
famAbun_clr_t <- as.data.frame(famAbun_clr_t)
family_names <- famAbun_clr_t[1,] 
famAbun_clr_t_num <- data.frame(lapply(famAbun_clr_t[-1,], as.numeric))
colnames(famAbun_clr_t_num) <- family_names  
rownames(famAbun_clr_t_num) <- rownames(famAbun_clr_t)[-1]

## Verification
str(famAbun_clr_t_num)

## Combine ARG and family CLR datasets
comb_ARG3_famAbunclr <- cbind(ARG3_abundance, famAbun_clr_t_num)
dim(comb_ARG3_famAbunclr)
str(comb_ARG3_famAbunclr)


##------------------------------------------ 16.2 NLW: CLR-transformed bacterial abundances by genus
library(compositions)
library(phyloseq)

## Subset ARGs and OTU table for northern lichen woodland (NLW)
str(ARG3_abundance) 
nlw_samples <- rownames(sample_data_df[sample_data_df$LW == "Kuujuarapik", ])
ARG3_abundance.NLW <- ARG3_abundance[rownames(ARG3_abundance) %in% nlw_samples, ]
otu_table_ps2.NLW <- as.data.frame(otu_table(ps2.NLW))
otu_table_ps2.NLW_nozero <- otu_table_ps2.NLW + 1  # Ajoute 1 à toutes les valeurs
otu_table_ps2.NLW_clr <- as.data.frame(clr(otu_table_ps2.NLW_nozero))

## Subset genus assignment for NLW ASVs
asvs_in_NLW <- colnames(otu_table_ps2.NLW_clr)
genusASV_Fam_NLW <- genusASV_Fam[rownames(genusASV_Fam) %in% asvs_in_NLW, , drop=FALSE]

## Sum CLR-transformed abundances per genus
genusAbun_clr.NLW <- aggregate(t(otu_table_ps2.NLW_clr), by = list(Genus = genusASV_Fam_NLW), FUN = sum)
genusAbun_clr.NLW_t <- t(genusAbun_clr.NLW)
genusAbun_clr.NLW_t <- as.data.frame(genusAbun_clr.NLW_t)
genus_names <- genusAbun_clr.NLW_t[1,] ## Extraire la première ligne avec la première colonne
genusAbun_clr.NLW_t_num <- data.frame(lapply(genusAbun_clr.NLW_t[-1,], as.numeric))
colnames(genusAbun_clr.NLW_t_num) <- genus_names  # Attribuer les noms des genus comme noms de colonnes
rownames(genusAbun_clr.NLW_t_num) <- rownames(genusAbun_clr.NLW_t)[-1]  # Conserver les noms de lignes

## Verification
str(genusAbun_clr.NLW_t_num)

## Combine ARGs and NLW genus CLR abundances
comb_ARG3_genusAbunclr.NLW <- cbind(ARG3_abundance.NLW, genusAbun_clr.NLW_t_num)
dim(comb_ARG3_genusAbunclr.NLW)
str(comb_ARG3_genusAbunclr.NLW)


##------------------------------------------ 16.3 SLW: CLR-transformed bacterial abundances by genus

library(compositions)
library(phyloseq)

## Subset ARGs and OTU table for southern lichen woodland (SLW)
str(ARG3_abundance) 
slw_samples <- rownames(sample_data_df[sample_data_df$LW == "PNGJ", ])
ARG3_abundance.SLW <- ARG3_abundance[rownames(ARG3_abundance) %in% slw_samples, ]
otu_table_ps2.SLW <- as.data.frame(otu_table(ps2.SLW))
otu_table_ps2.SLW_nozero <- otu_table_ps2.SLW + 1  # Ajoute 1 à toutes les valeurs
otu_table_ps2.SLW_clr <- as.data.frame(clr(otu_table_ps2.SLW_nozero))


## Subset genus assignment for SLW ASVs
asvs_in_SLW <- colnames(otu_table_ps2.SLW_clr)
genusASV_Fam_SLW <- genusASV_Fam[rownames(genusASV_Fam) %in% asvs_in_SLW, , drop=FALSE]

## Sum CLR-transformed abundances per genus
genusAbun_clr.SLW <- aggregate(t(otu_table_ps2.SLW_clr), by = list(Genus = genusASV_Fam_SLW), FUN = sum)
genusAbun_clr.SLW_t <- t(genusAbun_clr.SLW)
genusAbun_clr.SLW_t <- as.data.frame(genusAbun_clr.SLW_t)
genus_names <- genusAbun_clr.SLW_t[1,] ## Extraire la première ligne avec la première colonne
genusAbun_clr.SLW_t_num <- data.frame(lapply(genusAbun_clr.SLW_t[-1,], as.numeric))
colnames(genusAbun_clr.SLW_t_num) <- genus_names  # Attribuer les noms des genus comme noms de colonnes
rownames(genusAbun_clr.SLW_t_num) <- rownames(genusAbun_clr.SLW_t)[-1]  # Conserver les noms de lignes

## Verification
str(genusAbun_clr.SLW_t_num)

## Combine ARGs and SLW genus CLR abundances
comb_ARG3_genusAbunclr.SLW <- cbind(ARG3_abundance.SLW, genusAbun_clr.SLW_t_num)
dim(comb_ARG3_genusAbunclr.SLW)
str(comb_ARG3_genusAbunclr.SLW)



#####
##------------------------------------------ 17. NETWORK analyses
#####

library(Hmisc)
library(igraph)
library(ggplot2)
library(ggnetwork)

### Construct correlation matrix

## We calculate Spearman correlations between ARG relative abundances and bacterial genera/families
## Spearman correlation is suitable because it detects monotonic relationships, is robust to outliers, and does not assume normality

## Use combined dataset: comb_ARG3_genusAbunclr (ARGs + genera) or comb_ARG3_famAbunclr (ARGs + families)
cor_results <- rcorr(as.matrix(comb_ARG3_genusAbunclr), type = "spearman")
cor_matrix <- cor_results$r
p_matrix <- cor_results$P

## Adjust p-values for multiple testing using BH method
p_adj <- p.adjust(p_matrix, method = "BH")
dim(p_adj) <- dim(p_matrix)

## Apply thresholds for correlation and significance
cor_threshold <- 0.2
p_threshold <- 0.05
significant_cor <- (abs(cor_matrix) > cor_threshold) & (p_adj < p_threshold)
cor_matrix[!significant_cor] <- 0

## Construct undirected network with weighted edges
network <- graph_from_adjacency_matrix(cor_matrix, mode = "undirected", weighted = TRUE, diag = FALSE) 

## Create edge dataframe
edges <- as_data_frame(network, what = "edges")

## Add adjusted p-values
edges$p_value <- p_adj[cbind(match(edges$from, colnames(cor_matrix)), 
                             match(edges$to, colnames(cor_matrix)))]

## Annotate direction and strength of correlation
edges$direction <- ifelse(edges$weight > 0, "positive", "negative")
edges$strength <- cut(abs(edges$weight), 
                      breaks = c(0.2, 0.4, 0.6, 0.8, 1), # force de la correlation de Spearman
                      labels = c("weak", "moderate", "strong", "very strong"))

## Significant edges for Cytoscape filtering
edges$significant <- edges$p_value < 0.05

## Export edges
write.csv(edges, "network_edges.csv", row.names = FALSE)

## Create node dataframe
nodes <- data.frame(id = V(network)$name, 
                    label = V(network)$name,
                    type = ifelse(V(network)$name %in% colnames(ARG3_abundance), "ARG", "ASV"))

## Add ARG family information
info_3ARG <- data.frame(
  id = c("qnrB", "qepA", "blaCTXM1"),
  ARGFamily = c("Quinolone resistance", "Quinolone resistance", "Beta-lactam resistance")
)
nodes <- merge(nodes, info_3ARG, by = "id", all.x = TRUE)
nodes$ARGFamily[is.na(nodes$ARGFamily)] <- "Not Applicable"

## Exporter les nodes
write.csv(nodes, "network_nodes.csv", row.names = FALSE)


## we continued the analyses and visualization directly in Cytoscape.

