### Use corncob to look for differentially abundant taxa by the continuous variable of total tissue loss

library(phyloseq)
library(corncob)
library(plyr)
library(dplyr)
library(reshape2)
library(tibble)
library(ggplot2)


# use an ASV table that has been filtered to remove low-abundance ASVs, no blanks, replaced sequences with ASV numbers in both otu table and taxa table, added a column named "Species" to the taxa table to fill in later with ASV number
otu <- read.table("silva_nochloronomito_otu_table_ps5_BS2_jan_corncob.txt",sep="\t",header=TRUE, row.names=1)
taxon <- read.table("silva_nochloronomito_taxa_table_ps5_BS2_jan_corncob.txt",sep="\t",header=TRUE,row.names=1)
samples<-read.table("metadata_ps5_BS2_jan.txt",sep="\t",header=T,row.names=1)
OTU = otu_table(otu, taxa_are_rows=FALSE)
taxon<-as.matrix(taxon)
TAX = tax_table(taxon)
sampledata = sample_data(samples)
ps5_BS2_jan <- phyloseq(otu_table(otu, taxa_are_rows=FALSE), 
                        sample_data(samples), 
                        tax_table(taxon))
ps5_BS2_jan  #756 taxa and 59 samples

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

# All p-values are NA!




### Use corncob to look for differentially abundant taxa by category: treatment type

####### JANUARY
# use an ASV table that has been filtered to remove low-abundance ASVs, no blanks, replaced sequences with ASV numbers in both otu table and taxa table, added a column named "Species" to the taxa table to fill in later with ASV number
otu <- read.table("silva_nochloronomito_otu_table_ps5_BS2_jan_corncob.txt",sep="\t",header=TRUE, row.names=1)
taxon <- read.table("silva_nochloronomito_taxa_table_ps5_BS2_jan_corncob.txt",sep="\t",header=TRUE,row.names=1)
samples<-read.table("metadata_ps5_BS2_jan.txt",sep="\t",header=T,row.names=1)
OTU = otu_table(otu, taxa_are_rows=FALSE)
taxon<-as.matrix(taxon)
TAX = tax_table(taxon)
sampledata = sample_data(samples)
ps5_BS2_jan <- phyloseq(otu_table(otu, taxa_are_rows=FALSE), 
                        sample_data(samples), 
                        tax_table(taxon))
ps5_BS2_jan  #756 taxa and 59 samples


# What microorganisms are significantly different by treatment type? 
# DA determined relative to first reference level, so set this manually to the "none" treatment
sample_data(ps5_BS2_jan)$Treatment<-factor(sample_data(ps5_BS2_jan)$Treatment,levels=c("none","probiotic bag","probiotic paste","control bag","control paste"))

set.seed(1)
treatment.da <- differentialTest(formula = ~ Treatment, 
                                 phi.formula = ~ Treatment,
                                 formula_null = ~ 1,
                                 phi.formula_null = ~ Treatment,
                                 test = "Wald", boot = FALSE,
                                 data = ps5_BS2_jan,
                                 fdr_cutoff = 0.05)

summary(treatment.da)
treatment.da
treatment.da$significant_taxa #differentially abundant taxa
#ASV329

#quick look at results with corncob's plotting feature
plot(treatment.da, c("Family", "Genus"))
# Positive values mean that taxon is enriched in the treatment in the header relative to no treatment
# Only Fusibacter ASV329 was detected as DA: it was enriched in treatments relative to no treatment

# create a function that will return the values of interest
sigtaxto_df <- function(treatment.da, ps5_BS2_jan) {
  df.out <- c()
  for (i in 1:length(treatment.da$significant_models)) { #from 1 to the number of significant taxa
    df.out = rbind(df.out, treatment.da$significant_models[[i]]$coefficients[2, 1:4]) #bind the sig.mcav row with coefficient, std. error, t value, Pr. coefficients[2,] is the DA model output.
  }
  asv = as.vector(treatment.da$significant_taxa) #get significant ASVs
  Class = as.vector(otu_to_taxonomy(data = ps5_BS2_jan, treatment.da$significant_taxa, level = c("Class")))
  Order = as.vector(otu_to_taxonomy(data = ps5_BS2_jan, treatment.da$significant_taxa, level = c("Order")))
  Family = as.vector(otu_to_taxonomy(data = ps5_BS2_jan, treatment.da$significant_taxa, level = c("Family")))
  Genus = as.vector(otu_to_taxonomy(data = ps5_BS2_jan, treatment.da$significant_taxa, level = c("Genus")))
  
  df.bind <- as.data.frame(cbind(df.out, asv, Class, Order, Family, Genus))
  df.bind$asvgenus <- paste0(df.bind$Genus, " (", df.bind$asv, ")")
  colnames(df.bind) <- c("Estimate", "StdError", "tvalue", "Pr", "asv", "Class", "Order", "Family", "Genus", "asvtaxa")
  df.bind$Estimate = as.numeric(as.character(df.bind$Estimate)) #estimate column should be numeric 
  df.bind$StdError = as.numeric(as.character(df.bind$StdError)) #stderror column should be numeric 
  
  return(df.bind)
}

df <- sigtaxto_df(treatment.da, ps5_BS2_jan) #make a df of the corncob output and the taxonomy

### Save these data so you don't have to re-run the model ###
write.table(df, "DA_treatment-Jan.txt", sep="\t",row.names = FALSE)

# retrieve relative abundances of significant taxa for plotting
ps5_BS2_jan #756 taxa and 59 samples
ps5_BS2_jan_ra<-transform_sample_counts(ps5_BS2_jan, function(OTU) OTU/sum(OTU))
ps5_BS2_jan_ra #756 families and 59 samples
otu_da = as(otu_table(ps5_BS2_jan_ra), "matrix")
taxon_da = as(tax_table(ps5_BS2_jan_ra), "matrix")
da<-t(otu_da)
da <- as.data.frame(da)
da <- tibble::rownames_to_column(da, "ASV")
taxon_da <- as.data.frame(taxon_da)
taxon_da <- tibble::rownames_to_column(taxon_da, "ASV")
da_tax <- merge(da, taxon_da, by="ASV", all=FALSE)
da_tax$Species<-NULL
# now to retrieve only the taxa that I want to plot
DA1 <- da_tax %>% filter(ASV =='ASV329')


# Merge with metadata for plotting aesthetics
samples <- as.data.frame(samples)
samples <- tibble::rownames_to_column(samples, "sample")
DA_long<-melt(DA1,value.name="proportion",variable.name="sample",id.vars=c("Kingdom","Phylum","Class","Order","Family"))
DA_long <- DA_long[-1,]
DAplot <- merge(DA_long,samples,by="sample", all=TRUE)



####### OCTOBER
# use an ASV table that has been filtered to remove low-abundance ASVs, no blanks, replaced sequences with ASV numbers in both otu table and taxa table, added a column named "Species" to the taxa table to fill in later with ASV number
otu <- read.table("silva_nochloronomito_otu_table_ps5_BS2_oct_corncob.txt",sep="\t",header=TRUE, row.names=1)
taxon <- read.table("silva_nochloronomito_taxa_table_ps5_BS2_oct_corncob.txt",sep="\t",header=TRUE,row.names=1)
samples<-read.table("metadata_ps5_BS2_oct.txt",sep="\t",header=T,row.names=1)
OTU = otu_table(otu, taxa_are_rows=FALSE)
taxon<-as.matrix(taxon)
TAX = tax_table(taxon)
sampledata = sample_data(samples)
ps5_BS2_oct <- phyloseq(otu_table(otu, taxa_are_rows=FALSE), 
                        sample_data(samples), 
                        tax_table(taxon))
ps5_BS2_oct  #756 taxa and 66 samples


# What microorganisms are significantly different by treatment type?

# DA determined relative to first reference level, so set this manually to the "none" treatment
sample_data(ps5_BS2_oct)$Treatment<-factor(sample_data(ps5_BS2_oct)$Treatment,levels=c("none","probiotic bag","probiotic paste","control bag","control paste"))

set.seed(1)
treatment.da <- differentialTest(formula = ~ Treatment, 
                                 phi.formula = ~ Treatment,
                                 formula_null = ~ 1,
                                 phi.formula_null = ~ Treatment,
                                 test = "Wald", boot = FALSE,
                                 data = ps5_BS2_oct,
                                 fdr_cutoff = 0.05)

summary(treatment.da)
treatment.da
treatment.da$significant_taxa #differentially abundant taxa
#ASV67

#quick look at results with corncob's plotting feature
plot(treatment.da, c("Family", "Genus"))
# Positive values mean that taxon is enriched in the treatment in the header relative to no treatment
# Only Ottowia ASV67 was detected as DA: it was enriched in treatments compared to no treatment

# create a function that will return the values of interest

sigtaxto_df <- function(treatment.da, ps5_BS2_oct) {
  df.out <- c()
  for (i in 1:length(treatment.da$significant_models)) { #from 1 to the number of significant taxa
    df.out = rbind(df.out, treatment.da$significant_models[[i]]$coefficients[2, 1:4]) #bind the significant taxa with coefficient, std. error, t value, Pr. coefficients[2,] is the DA model output.
  }
  asv = as.vector(treatment.da$significant_taxa) #get significant ASVs
  Class = as.vector(otu_to_taxonomy(data = ps5_BS2_oct, treatment.da$significant_taxa, level = c("Class")))
  Order = as.vector(otu_to_taxonomy(data = ps5_BS2_oct, treatment.da$significant_taxa, level = c("Order")))
  Family = as.vector(otu_to_taxonomy(data = ps5_BS2_oct, treatment.da$significant_taxa, level = c("Family")))
  Genus = as.vector(otu_to_taxonomy(data = ps5_BS2_oct, treatment.da$significant_taxa, level = c("Genus")))
  
  df.bind <- as.data.frame(cbind(df.out, asv, Class, Order, Family, Genus))
  df.bind$asvgenus <- paste0(df.bind$Genus, " (", df.bind$asv, ")")
  colnames(df.bind) <- c("Estimate", "StdError", "tvalue", "Pr", "asv", "Class", "Order", "Family", "Genus", "asvtaxa")
  df.bind$Estimate = as.numeric(as.character(df.bind$Estimate)) #estimate column should be numeric 
  df.bind$StdError = as.numeric(as.character(df.bind$StdError)) #stderror column should be numeric 
  
  return(df.bind)
}

df <- sigtaxto_df(treatment.da, ps5_BS2_oct) #make a df of the corncob output and the taxonomy

### Save these data so you don't have to re-run the model ###
write.table(df, "DA_treatment-Oct.txt", sep="\t",row.names = FALSE)






