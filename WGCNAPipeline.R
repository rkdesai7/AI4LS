#Load WGCNA
if (!("DESeq2" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("DESeq2", update = FALSE)
}

if (!("impute" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("impute")
}

if (!("WGCNA" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  install.packages("WGCNA")
}

if (!("ggforce" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  install.packages("ggforce")
}

if (!("ComplexHeatmap" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  install.packages("ComplexHeatmap")
}

# Attach the DESeq2 library
library(DESeq2)

# We will need this so we can use the pipe: %>%
library(magrittr)

# We'll need this for finding gene modules
library(WGCNA)

# We'll be making some plots
library(ggplot2)

library(WGCNA)
options(stringsAsFactors = FALSE)

#Filename
data = 'GLDS105_transcriptomic_normalized_VST.csv'
GLDS105_transcriptomic_normalized_VST_1_ <- read.csv(data,header=TRUE,row.names=1,check.names=F)
df = t(GLDS105_transcriptomic_normalized_VST_1_)

#Find soft threshold power
sft = pickSoftThreshold(df, dataIsExpr = TRUE, corFnc = cor, networkType = "signed")
sft_df <- data.frame(sft$fitIndices)%>%
  dplyr::mutate(model_fit = -sign(slope)*SFT.R.sq)

ggplot(sft_df, aes(x = Power, y = model_fit, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.80, col = "red") +
  ylim(c(min(sft_df$model_fit), 1.05)) +
  xlab("Soft Threshold (power)") +
  ylab("Scale Free Topology Model Fit, signed R^2") +
  ggtitle("Scale independence") +
  theme_classic()


bwnet <- blockwiseModules(datExpr = df,   maxBlockSize = 5000,
                          TOMType = "signed", 
                          power = 20, 
                          numericLabels = TRUE, 
                          randomSeed = 1234)

#Write main WGNCA results object to file
readr::write_rds(bwnet, file = file.path("results", "GLDS105Transcriptomics_WGCNA_results.RDS"))

metadata <- df %>%
  dplyr::mutate(
    time_point = dplyr::case_when(
      # Create our new variable based on refinebio_title containing AV/CV
      stringr::str_detect(df[,0], "_FLT_") ~ "flight",
      stringr::str_detect(df[,0], "_GC_") ~ "ground control"
      
    ),
    time_point = as.factor(time_point)
  )
    
module_df <- data.frame(gene_id = names(bwnet$colors), colors = labels2colors(bwnet$colors))
mergedColors = labels2colors(bwnet$colors)
#getmodule eigengenes per cluster
MEs0 <- moduleEigengenes(df, mergedColors)$eigengenes
#reorder modules so similar modules are next to each other
MEs0 <- orderMEs(MEs0)
module_order = names(MEs0) %>% gsub("ME", "", .)
#Add treatment names
MEs0$treatment = row.names(MEs0)
#tidy and plot data
mME = MEs0 %>% pivot_longer(-treatment) %>% mutate(name = gsub("ME", "", name), 
                                                   name = factor(name, levels = module_order))
#Finding modules of interest
mME %>% ggplot(., aes(x=treatment, y=name, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-trait Relationships", y = "Modules", fill="corr")
head(MEs0)


modules_of_interest = c("orangered4", "firebrick4", "maroon", "blue", "skyblue1","ivory")

genes_of_interest = module_df %>% subset(colors %in% modules_of_interest)


