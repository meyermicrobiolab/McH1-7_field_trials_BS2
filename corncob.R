### Use corncob to look for differentially abundant taxa by the continuous variable of total tissue loss

# use an ASV table that has been filtered to remove low-abundance ASVs, no blanks, replaced sequences with ASV numbers in both otu table and taxa table, added a column named "Species" to the taxa table to fill in later with ASV number
otu <- read.table("silva_nochloronomito_otu_table_ps5_BS2_jan corncob.txt",sep="\t",header=TRUE, row.names=1)
taxon <- read.table("silva_nochloronomito_taxa_table_ps5_BS2_jan_corncob.txt",sep="\t",header=TRUE,row.names=1)
samples<-read.table("metadata_ps5_BS2_jan.txt",sep="\t",header=T,row.names=1)
OTU = otu_table(otu, taxa_are_rows=FALSE)
taxon<-as.matrix(taxon)
TAX = tax_table(taxon)
sampledata = sample_data(samples)
ps5_BS2_jan <- phyloseq(otu_table(otu, taxa_are_rows=FALSE), 
                        sample_data(samples), 
                        tax_table(taxon))
ps5_BS2_jan  #756 taxa and 63 samples

ps_corn = subset_samples(ps5_BS2_jan, Proportion_TissueLost_Jan21 != "nd")
ps_corn #756 taxa and 54 samples
ps_corn <- prune_taxa(taxa_sums(ps_corn) > 1, ps_corn) #remove taxa not found in remaining samples
ps_corn #651 taxa and 54 samples, removed 105 taxa only found in samples without tissue loss data



# What microorganisms are significantly changing with changes in tissue loss? Use Jan samples only; total tissue loss in January

set.seed(1)
tissueJan.da <- differentialTest(formula = ~ Proportion_TissueLost_Jan21, 
                                 phi.formula = ~ Proportion_TissueLost_Jan21,
                                 formula_null = ~ 1,
                                 phi.formula_null = ~ Proportion_TissueLost_Jan21,
                                 test = "LRT", boot = FALSE,
                                 data = ps_corn,
                                 fdr_cutoff = 0.05)

# model overspecified - too many taxa compared to samples


```


### Use corncob to look for differentially abundant taxa by category: treatment type


```{r, echo=FALSE}
# use an ASV table that has been filtered to remove low-abundance ASVs, no blanks, replaced sequences with ASV numbers in both otu table and taxa table, added a column named "Species" to the taxa table to fill in later with ASV number
otu <- read.table("silva_nochloronomito_otu_table_ps5_BS2_jan corncob.txt",sep="\t",header=TRUE, row.names=1)
taxon <- read.table("silva_nochloronomito_taxa_table_ps5_BS2_jan_corncob.txt",sep="\t",header=TRUE,row.names=1)
samples<-read.table("metadata_ps5_BS2_jan.txt",sep="\t",header=T,row.names=1)
OTU = otu_table(otu, taxa_are_rows=FALSE)
taxon<-as.matrix(taxon)
TAX = tax_table(taxon)
sampledata = sample_data(samples)
ps5_BS2_jan <- phyloseq(otu_table(otu, taxa_are_rows=FALSE), 
                        sample_data(samples), 
                        tax_table(taxon))
ps5_BS2_jan  #756 taxa and 63 samples


# What microorganisms are significantly different by treatment type? Use Jan samples only first since the files are already formatted

# DA determined relative to first reference level, so set this manually to the "none" treatment
sample_data(ps5_BS2_jan)$Treatment<-factor(sample_data(ps5_BS2_jan)$Treatment,levels=c("none","probiotic bag","probiotic paste","control bag","control paste","resistant"))


set.seed(1)
treatment.da <- differentialTest(formula = ~ Treatment, 
                                 phi.formula = ~ Treatment,
                                 formula_null = ~ 1,
                                 phi.formula_null = ~ Treatment,
                                 test = "Wald", boot = FALSE,
                                 data = ps5_BS2_jan,
                                 fdr_cutoff = 0.05)


#quick look at results with corncob's plotting feature
plot(treatment.da)
plot(treatment.da, level="Genus")

# create a function that will return the values of interest
# This one is for when there is only one comparison
extractmods <- function(model) {
  result <- data.frame("Estimate" = model$coefficients[2, 1], 
                       "Std.Error" = model$coefficients[2, 2], 
                       "p" = model$coefficients[2, 4])
  return(result)
}

# save the treatment.da output and add asv, p.adj values to it, and ultimately taxonomny info
treatment.da.models <- lapply(treatment.da$significant_models, extractmods)
names(treatment.da.models) <- treatment.da$significant_taxa

# Add ASVs to the taxonomy table and save the significant asvs
tax_table(ps5_BS2_jan)[,7] <- rownames(tax_table(ps5_BS2_jan))
sig.taxonomy.treatment <- as.data.frame(tax_table(ps5_BS2_jan)[treatment.da$significant_taxa,]) 

# Move the data from a list to a dataframe and add taxonomy info
treatment.da.models.df <- ldply(treatment.da.models, data.frame) %>% 
  left_join(sig.taxonomy.treatment, by = c(".id" = "Species")) %>%
  mutate(genusasv = paste0(Genus, "_(", .id,")"))

### Save these data so you don't have to re-run the model ###
write.table(treatment.da.models.df, "DA_treatment.txt", sep="\t",row.names = FALSE)

