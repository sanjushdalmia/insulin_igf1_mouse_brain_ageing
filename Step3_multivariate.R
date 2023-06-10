
# Load packages ------
library(tidyverse) # you're familiar with this fromt the past two lectures
library(DT) # for making interactive tables
library(plotly) # for making interactive plots
library(gt) # A layered 'grammar of tables' - think ggplot, but for tables
library(ggeasy)

# Identify variables of interest in study design file ----

batch <- studydesign$Batch
batch.factor <- factor(batch)

cohort <- studydesign$Cohort
cohort.factor <- factor(cohort)

age <- studydesign$Age
age.factor <- factor(age)

genotype <- studydesign$Genotype
genotype.factor <- factor(genotype)

age_genotype <- studydesign$Age_Genotype
age_genotype.factor <- factor(age_genotype)


# Prepare your data -------
# for this part of the class you'll use your normalized and filtered data in log2 cpm
# make sure you have this object already in your work environment
# if you don't, go back to the Step2 script and generate it
log2.cpm.filtered.norm.df
sampleLabels <- paste(studydesign$Sample)
# Hierarchical clustering ---------------
#hierarchical clustering can only work on a data matrix, not a data frame, so we're using the data matrix version of our data rather than the data frame version
#try using filtered and unfiltered data...how does this change the result?
#try other distance methods (e.g. switch from 'maximum' to 'euclidean')...how does this change the result?
distance <- dist(t(log2.cpm.filtered.norm), method = "maximum") #other distance methods are "euclidean", maximum", "manhattan", "canberra", "binary" or "minkowski"
clusters <- hclust(distance, method = "average") #other agglomeration methods are "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", or "centroid"
plot(clusters, labels=sampleLabels)

# Principal component analysis (PCA) -------------
pca.res <- prcomp(t(log2.cpm.filtered.norm), scale.=F, retx=T)
#look at the PCA result (pca.res) that you just created
ls(pca.res)
summary(pca.res) # Prints variance summary for all principal components.
pca.res$rotation #$rotation shows you how much each gene influenced each PC (called 'scores')
pca.res$x # 'x' shows you how much each sample influenced each PC (called 'loadings')
#note that these have a magnitude and a direction (this is the basis for making a PCA plot)
screeplot(pca.res) # A screeplot is a standard way to view eigenvalues for each PCA
pc.var<-pca.res$sdev^2 # sdev^2 captures these eigenvalues from the PCA result
pc.per<-round(pc.var/sum(pc.var)*100, 1) # we can then use these eigenvalues to calculate the percentage variance explained by each PC
pc.per

pc.pertibble<-as_tibble(pc.per)

top10pc.per<-top_n(pc.pertibble, 10)

gt(top10pc.per)


# Visualize your PCA result ------------------
#lets first plot any two PCs against each other
#We know how much each sample contributes to each PC (loadings), so let's plot
pca.res.df <- as_tibble(pca.res$x)

# PC1 vs PC2 plots --------------------------------------------------------

# Plain
PC12plain <- ggplot(pca.res.df) +
  aes(x=PC1, y=PC2, label=sampleLabels) +
  geom_point(size=2) +
  xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  labs(title="Scatterplot: PC1 and PC2",
       caption=paste0("produced on ", Sys.time())) +
  coord_fixed() +
  theme_bw() +
  theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        )

# Age
PC12age <- ggplot(pca.res.df) +
  aes(x=PC1, y=PC2, label=sampleLabels, color = age.factor) +
  geom_point(size=2) +
  xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  labs(title="Scatterplot: PC1 and PC2",
       color = "Age",
       caption=paste0("produced on ", Sys.time())) +
  coord_fixed() +
  theme_bw() +
  theme(
        plot.title = element_text(face = "bold", hjust = 0.5)
        )

# Genotype
PC12genotype <- ggplot(pca.res.df) +
  aes(x=PC1, y=PC2, label=sampleLabels, color = genotype.factor) +
  geom_point(size=2) +
  xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  labs(title="Scatterplot: PC1 and PC2",
       color = "Genotype",
       caption=paste0("produced on ", Sys.time())) +
  coord_fixed() +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        )

# Age_Genotype
PC12age_genotype <- ggplot(pca.res.df) +
  aes(x=PC1, y=PC2, label=sampleLabels, color = age_genotype.factor) +
  geom_point(size=2) +
  xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  labs(title="Scatterplot: PC1 and PC2",
       color = "Age and Genotype",
       caption=paste0("produced on ", Sys.time())) +
  coord_fixed() +
  theme_bw() +
  theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        )

# Batch
PC12batch <- ggplot(pca.res.df) +
  aes(x=PC1, y=PC2, label=sampleLabels, color = batch.factor) +
  geom_point(size=2) +
  xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  labs(title="Scatterplot: PC1 and PC2",
       color = "Batch",
       caption=paste0("produced on ", Sys.time())) +
  coord_fixed() +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
  )

PC12 <- plot_grid(PC12age, PC12genotype, PC12age_genotype, PC12batch, labels = c('A', 'B', 'C', 'D'), label_size = 12)

# PC2 vs PC3 plots --------------------------------------------------------

# Plain
PC23plain <- ggplot(pca.res.df) +
  aes(x=PC2, y=PC3, label=sampleLabels) +
  geom_point(size=2) +
  xlab(paste0("PC2 (",pc.per[2],"%",")")) + 
  ylab(paste0("PC3 (",pc.per[3],"%",")")) +
  labs(title="Scatterplot: PC2 and PC3",
       caption=paste0("produced on ", Sys.time())) +
  coord_fixed() +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
  )

# Age
PC23age <- ggplot(pca.res.df) +
  aes(x=PC2, y=PC3, label=sampleLabels, color = age.factor) +
  geom_point(size=2) +
  xlab(paste0("PC2 (",pc.per[2],"%",")")) + 
  ylab(paste0("PC3 (",pc.per[3],"%",")")) +
  labs(title="Scatterplot: PC2 and PC3",
       color = "Age",
       caption=paste0("produced on ", Sys.time())
      ) +
  coord_fixed() +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

# Genotype
PC23genotype <- ggplot(pca.res.df) +
  aes(x=PC2, y=PC3, label=sampleLabels, color = genotype.factor) +
  geom_point(size=2) +
  xlab(paste0("PC2 (",pc.per[2],"%",")")) + 
  ylab(paste0("PC3 (",pc.per[3],"%",")")) +
  labs(title="Scatterplot: PC2 and PC3",
       color = "Genotype",
       caption=paste0("produced on ", Sys.time())
       ) +
  coord_fixed() +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
  )

# Age_Genotype
PC23age_genotype <- ggplot(pca.res.df) +
  aes(x=PC2, y=PC3, label=sampleLabels, color = age_genotype.factor) +
  geom_point(size=2) +
  xlab(paste0("PC2 (",pc.per[2],"%",")")) + 
  ylab(paste0("PC3 (",pc.per[3],"%",")")) +
  labs(title="Scatterplot: PC2 and PC3",
       color = "Age and Genotype",
       caption=paste0("produced on ", Sys.time())
      ) +
  coord_fixed() +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
  )

# Batch
PC23batch <- ggplot(pca.res.df) +
  aes(x=PC2, y=PC3, label=sampleLabels, color = batch.factor) +
  geom_point(size=2) +
  xlab(paste0("PC2 (",pc.per[2],"%",")")) + 
  ylab(paste0("PC3 (",pc.per[3],"%",")")) +
  labs(title="Scatterplot: PC2 and PC3",
       color = "Batch",
       caption=paste0("produced on ", Sys.time())
      ) +
  coord_fixed() +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
  )


PC23 <- plot_grid(PC23age, PC23genotype, PC23age_genotype, PC23batch, labels = c('A', 'B', 'C', 'D'), label_size = 12)
# Let's discuss and iteratively refine the PCA code and plot from above
# First, take note of the fact that we can use information from our PCA analysis to label our axes
# Remember that PCA is unsupervised, so knows nothing about group assignment (healthy vs disease)
# But *we* know, and so we can use this knowledge to enhance the plot.  Add a 'color=group' mapping to the aes of the plot above
# Can we figure out the identity of the outlier?  We have already provided samplelabel mapping in aes, so just uncomment the 'geom_label()'
# Uncomment 'coord_fixed()' to apply the correct aspect ratio
# Uncomment 'stat_ellipse()' to see how you can circle clusters on the PCA
# How would this PCA look if you used raw counts (myCounts) instead of log2 CPM?
# What are the disadvantages of looking at a PCA result using such a simple XY plot?

# Use dplyr 'verbs' to modify our dataframe ----
# use dplyr 'mutate' function to add new columns based on existing data

FilterWToldminusDKOold.df <- log2.cpm.filtered.norm.df %>% 
  mutate(WTold.AVG.column = (S10 + S11 + S12)/3,
         DKOold.AVG.column = (S01 + S02 + S03)/3,
         #now make columns comparing each of the averages above that you're interested in
         FilterWToldminusDKOold.column = (WTold.AVG.column - DKOold.AVG.column)) %>% 
  mutate_if(is.numeric, round, 2)

FilterWTyoungminusDKOyoung.df <- log2.cpm.filtered.norm.df %>% 
  mutate(WTyoung.AVG.column = (S22 + S23 + S24)/3,
         DKOyoung.AVG.column = (S13 + S14 + S15)/3,
         #now make columns comparing each of the averages above that you're interested in
         FilterWTyoungminusDKOyoung.column = (WTyoung.AVG.column - DKOyoung.AVG.column)) %>% 
  mutate_if(is.numeric, round, 2)

Filterinteraction.df <- log2.cpm.filtered.norm.df %>% 
  mutate(WTold.AVG.column = (S10 + S11 + S12)/3,
         DKOold.AVG.column = (S01 + S02 + S03)/3,
         WTyoung.AVG.column = (S22 + S23 + S24)/3,
         DKOyoung.AVG.column = (S13 + S14 + S15)/3,
         #now make columns comparing each of the averages above that you're interested in
         Filterinteraction.column = (WTold.AVG.column - WTyoung.AVG.column) - (DKOold.AVG.column -  DKOyoung.AVG.column))
         

#now make columns comparing each of the averages above that you're interested in
        

#now look at this modified data table


# Use dplyr 'arrange' and 'select' to sort your dataframe based on any variable
# first, we'll use dplyr "arrange" function to sort rows based on the values in a column of interest
# then we'll display 'select' only the columns we're interested in seeing



# Use dplyr "filter" and "select" functions to pick out genes of interest 
# ways to tweak the 'select' function:
# use ':' between two column names to select all columns between
# use 'contains', 'starts_with' or 'ends_with' to modify how you select
# can refer to columns using exact name or numerical indicator
# use boolean operators such as '&' (and), '|' (or), '==' (equal to), '!' (not)



