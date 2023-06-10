
# Load packages -----
library(tidyverse) # already know about this from Step 1 script
library(edgeR) # well known package for differential expression analysis, but we only use for the DGEList object and for normalization methods
library(matrixStats) # let's us easily calculate stats on rows or columns of a data matrix
library(cowplot) # allows you to combine multiple plots in one figure
library (openxlsx)
library (gt)

# Examine your data up to this point ----
myTPM <- Txi_gene$abundance
myCounts <- Txi_gene$counts
colSums(myTPM)
colSums(myCounts)

# capture sample labels from the study design file that you worked with and saved as 'targets' in step 1

sampleLabels <- studydesign$Sample 

# Make a DGElist from your counts, and plot 
myDGEList <- DGEList(myCounts)
# take a look at the DGEList object 
myDGEList
#DEGList objects are a good R data file to consider saving to you working directory
save(myDGEList, file = "myDGEList")
#Saved DGEList objects can be easily shared and loaded into an R environment
load(file = "myDGEList")

# use the 'cpm' function from EdgeR to get counts per million
cpm <- cpm(myDGEList) 
colSums(cpm)
log2.cpm <- cpm(myDGEList, log=TRUE)

# 'coerce' your data matrix to a dataframe so that you can use tidyverse tools on it
log2.cpm.df <- as_tibble(log2.cpm, rownames = "geneID")
log2.cpm.df
# add your sample names to this dataframe (we lost these when we read our data in with tximport)
colnames(log2.cpm.df) <- c("geneID", sampleLabels)
# use the tidy package to 'pivot' your dataframe (from wide to long)
log2.cpm.df.pivot <- pivot_longer(log2.cpm.df, # dataframe to be pivoted
                                  cols = S01:S24, # column names to be stored as a SINGLE variable
                                  names_to = "samples", # name of that new variable (column)
                                  values_to = "expression") # name of new variable (column) storing all the values (data)

# let's look at the impact of pivoting the data
log2.cpm.df.pivot

# note it is easy to plot this pivoted data
p1 <- ggplot(log2.cpm.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="Log2 Expression", x = "Sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="unfiltered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        plot.subtitle = element_text(face = "italic"),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14)) 


# Filter your data ----
#first, take a look at how many genes or transcripts have no read counts at all
table(rowSums(myDGEList$counts==0)==24)
# breaking down the line above is a little tricky.  Let's try:
# 1st - 'myDGEList$counts==0' returns a new 'logical matrix' where each observation (gene) is evaluated (TRUE/FALSE) for each variable (sample) as to whether it has zero counts
# 2nd - passing this logical matrix to 'rowsums' allows you to sum the total number of times an observation was 'TRUE' across all samples
# 3rd - adding the '==10' is a simple way of asking how many of the rowsums equaled 10. In other words, how many genes had 0 counts (TRUE) for all samples in our dataset
# 4th - passing all this to the 'table' function just provides a handy way to summarize the large logical produced in the previous step

# now set some cut-off to get rid of genes/transcripts with low counts
# again using rowSums to tally up the 'TRUE' results of a simple evaluation
# how many genes had more than 1 CPM (TRUE) in at least 3 samples

# The line below is important! This is where the filtering starts
# Be sure to adjust this cutoff for the number of samples in the smallest group of comparison.
keepers <- rowSums(cpm>1)>=3
# now use base R's simple subsetting method to filter your DGEList based on the logical produced above
myDGEList.filtered <- myDGEList[keepers,]
dim(myDGEList.filtered)

log2.cpm.filtered <- cpm(myDGEList.filtered, log=TRUE)
log2.cpm.filtered.df <- as_tibble(log2.cpm.filtered, rownames = "geneID")
colnames(log2.cpm.filtered.df) <- c("geneID", sampleLabels)
# pivot this FILTERED data, just as you did earlier
log2.cpm.filtered.df.pivot <- pivot_longer(log2.cpm.filtered.df, # dataframe to be pivoted
                                           cols = S01:S24, # column names to be stored as a SINGLE variable
                                           names_to = "samples", # name of that new variable (column)
                                           values_to = "expression") # name of new variable (column) storing all the values (data)


p2 <- ggplot(log2.cpm.filtered.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="Log2 Expression", x = "Sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        plot.subtitle = element_text(face = "italic"),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14)) 

# Normalize your data ----
myDGEList.filtered.norm <- calcNormFactors(myDGEList.filtered, method = "TMM")
# take a look at this new DGEList object...how has it changed?

# use the 'cpm' function from EdgeR to get counts per million from your normalized data
log2.cpm.filtered.norm <- cpm(myDGEList.filtered.norm, log=TRUE)
log2.cpm.filtered.norm.df <- as_tibble(log2.cpm.filtered.norm, rownames = "geneID")
colnames(log2.cpm.filtered.norm.df) <- c("geneID", sampleLabels)
# pivot this NORMALIZED data, just as you did earlier
log2.cpm.filtered.norm.df.pivot <- pivot_longer(log2.cpm.filtered.norm.df, # dataframe to be pivoted
                                                cols = S01:S24, # column names to be stored as a SINGLE variable
                                                names_to = "samples", # name of that new variable (column)
                                                values_to = "expression") # name of new variable (column) storing all the values (data)


p3 <- ggplot(log2.cpm.filtered.norm.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="Log2 Expression", x = "Sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        plot.subtitle = element_text(face = "italic"),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14)) 


# we'll use the 'plot_grid' function from the cowplot package to put these together in a figure
p4 <- plot_grid(p1, p2, p3, labels = c('A', 'B', 'C'), label_size = 12)


# Filter Counts and TMP data to verify IGF-1 and Insr KO ------------------------------------------------------

#Using Counts, WT Old vs DKO Old

myCounts.tibble <- as_tibble(myCounts, rownames = "geneID")

Filter.myCounts.tibble <- myCounts.tibble %>% 
  mutate(WTold.AVG.column = (V10 + V11 + V12)/3,
         IGF1RKOold.AVG.column = (V4 + V5 + V6) / 3,
         IRKOold.AVG.column = (V7 + V8 + V9) / 3,
         DKOold.AVG.column = (V1 + V2 + V3)/3,
         WTyoung.AVG.column = (V22 + V23 + V24)/3,
         IGF1RKOyoung.AVG.column = (V16 + V17 + V18)/3,
         IRKOyoung.AVG.column = (V19 + V20 + V21)/3,
         DKOyoung.AVG.column = (V13 + V14 + V15)/3) %>% 
  mutate_if(is.numeric, round, 2)

Filter.myCounts.filter.old <- Filter.myCounts.tibble %>%
  dplyr::filter(geneID=="Insr" | geneID=="Igf1r" ) %>%
  dplyr::select(geneID, 
                WTold.AVG.column, 
                IGF1RKOold.AVG.column,  
                IRKOold.AVG.column, 
                DKOold.AVG.column)

gt(Filter.myCounts.filter.old)

Filter.myCounts.filter.young <- Filter.myCounts.tibble %>%
  dplyr::filter(geneID=="Insr" | geneID=="Igf1r" ) %>%
  dplyr::select(geneID, 
                WTyoung.AVG.column,
                IGF1RKOyoung.AVG.column,
                IRKOyoung.AVG.column,
                DKOyoung.AVG.column)

gt(Filter.myCounts.filter.young)


#Again using TPM

myTPM.tibble <- as_tibble(myTPM, rownames = "geneID")

Filter.myTPM.tibble <- myTPM.tibble %>% 
  mutate(WTold.AVG.column = (V10 + V11 + V12)/3,
         IGF1RKOold.AVG.column = (V4 + V5 + V6) / 3,
         IRKOold.AVG.column = (V7 + V8 + V9) / 3,
         DKOold.AVG.column = (V1 + V2 + V3)/3,
         WTyoung.AVG.column = (V22 + V23 + V24)/3,
         IGF1RKOyoung.AVG.column = (V16 + V17 + V18)/3,
         IRKOyoung.AVG.column = (V19 + V20 + V21)/3,
         DKOyoung.AVG.column = (V13 + V14 + V15)/3) %>% 
  mutate_if(is.numeric, round, 2)

Filter.myTPM.filter.old <- Filter.myTPM.tibble %>%
  dplyr::filter(geneID=="Insr" | geneID=="Igf1r" ) %>%
  dplyr::select(geneID, 
                WTold.AVG.column, 
                IGF1RKOold.AVG.column,  
                IRKOold.AVG.column, 
                DKOold.AVG.column, 
                ) 
  
gt(Filter.myTPM.filter.old)

Filter.myTPM.filter.young <- Filter.myTPM.tibble %>%
  dplyr::filter(geneID=="Insr" | geneID=="Igf1r" ) %>%
  dplyr::select(geneID, 
                WTyoung.AVG.column,
                IGF1RKOyoung.AVG.column,
                IRKOyoung.AVG.column,
                DKOyoung.AVG.column) 

gt(Filter.myTPM.filter.young)




# TPM data for DGEs -------------------------------------------------------

myTPM.tibble <- as_tibble(myTPM, rownames = "geneID")

Filter.myTPM.tibble <- myTPM.tibble %>% 
  mutate(WTold.AVG.column = (V10 + V11 + V12)/3,
         DKOold.AVG.column = (V1 + V2 + V3)/3,
         WTyoung.AVG.column = (V22 + V23 + V24)/3,
         DKOyoung.AVG.column = (V13 + V14 + V15)/3,
         WTaging.AVG.column = WTold.AVG.column - WTyoung.AVG.column,
         DKOaging.AVG.column = DKOold.AVG.column - DKOyoung.AVG.column) %>% 
  mutate_if(is.numeric, round, 2)

Filter.myTPM.DGEs <- Filter.myTPM.tibble %>%
  dplyr::filter(geneID=="Cops7b" | geneID=="Fam171a2" | 
                geneID=="Bai2" | geneID=="Parp14" |
                  geneID=="Ensa" | geneID=="Hspa1a" |
                  geneID=="Rasl11b" | geneID=="Rasal1" |
                  geneID=="Junb" | geneID=="Ctdsp1") %>%
  dplyr::select(geneID, V10, V11, V12, V1, V2, V3, V22, V23, V24, V13, V14, V15) 

Filter.myTPM.DGEs %>%
gt() 

openxlsx::write.xlsx(Filter.myTPM.DGEs, "TopDGEs.xlsx")
