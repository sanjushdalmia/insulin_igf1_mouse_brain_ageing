# Introduction to this script -----------
# the goal of this script is to identify differentially expressed genes (DEGs) and differential transcript usage (DTU)
# you should already know which pairwise comparisons are most important to you
# whether you look for differential expression at the gene or transcript level depends on how you read the Kallisto output into R using TxImport back in Step 1
# if you have no biological replicates, you will NOT be able to leverage statistical tools for differential expression analysis
# instead, you will ONLY rely on fold changes, and can use the dplyr 'verbs' we discussed in Step 3 and 4 to identify genes based on log fold-change

# Load packages -----
library(tidyverse) # you know it well by now!
library(limma) # venerable package for differential gene expression using linear modeling
library(edgeR)
library(gt)
library(DT)
library(plotly)
library(gtExtras)

# Set up your design matrix ----
age_genotype.factor <- factor(age_genotype)
batch.factor <- factor(batch)

# adding batch as a blocking factor was messing up the matrix for some reason, come back to this

designinteraction <- model.matrix(~0+age_genotype.factor)
colnames(designinteraction) <- levels(age_genotype.factor)


# NOTE: if you need a paired analysis (a.k.a.'blocking' design) or have a batch effect, the following design is useful
# design <- model.matrix(~block + treatment)
# this is just an example. 'block' and 'treatment' would need to be objects in your environment

# Model mean-variance trend and fit linear model to data ----
# Use VOOM function from Limma package to model the mean-variance relationship and then correct it / weight genes by the relationship

v.DEGList.filtered.norm <- voom(myDGEList.filtered.norm, designinteraction, plot = TRUE)

fit <- lmFit(v.DEGList.filtered.norm, designinteraction)


# Contrast matrix ----
interaction.contrast.matrix <- makeContrasts(
                                DKOoldminusWTold = Old.DKO - Old.WT,
                                DKOYoungminusWTYoung = Young.DKO - Young.WT,
                              Interaction = (Old.DKO - Young.DKO) - (Old.WT - Young.WT),
                                 levels=designinteraction)

# make.Contrasts means make comparisons
# extract the linear model fit -----
interaction.fits <- contrasts.fit(fit, interaction.contrast.matrix)


#get bayesian stats for your linear model fit
ebFit.interaction <- eBayes(interaction.fits)
write.fit(ebFit.interaction, file = "interactionlmfit_results")



# TopTable to view DEGs -----
DKOoldminusWTold_allDEGs <- topTable(ebFit.interaction, adjust ="BH", coef=1, number=40000, sort.by="p")

DKOoldminusWTold_top10DEGs <- topTable(ebFit.interaction, adjust ="BH", coef=1, number=10, sort.by="p")

DKOYoungminusWTYoung_allDEGs <- topTable(ebFit.interaction, adjust ="BH", coef=2, number=40000, sort.by="p")

DKOYoungminusWTYoung_top10DEGs <- topTable(ebFit.interaction, adjust ="BH", coef=2, number=10, sort.by="p")

Interaction_allDEGs <- topTable(ebFit.interaction, adjust ="BH", coef=3, number=40000, sort.by="p")

Interaction_top10DEGs <- topTable(ebFit.interaction, adjust ="BH", coef=3, number=10, sort.by="p")

Interactions_DEGsforGProfiler <- topTable(ebFit.interaction, adjust ="BH", coef=3, number=40000, sort.by="p", "p.value" = 0.05)

# convert to a tibble
DKOoldminusWTold_top10DEGs.tibble <- DKOoldminusWTold_top10DEGs %>%
  as_tibble(rownames = "geneID")

DKOYoungminusWTYoung_top10DEGs.tibble <- DKOYoungminusWTYoung_top10DEGs %>%
  as_tibble(rownames = "geneID")

Interaction_top10DEGs.tibble <- Interaction_top10DEGs %>%
  as_tibble(rownames = "geneID")

DKOoldminusWTold_allDEGs.tibble <- DKOoldminusWTold_allDEGs %>%
  as_tibble(rownames = "geneID")

DKOYoungminusWTYoung_allDEGs.tibble <- DKOYoungminusWTYoung_allDEGs %>%
  as_tibble(rownames = "geneID")

Interaction_allDEGs.tibble <- Interaction_allDEGs %>%
  as_tibble(rownames = "geneID")

Interactions_DEGsforGProfiler.tibble <- Interactions_DEGsforGProfiler %>%
  as_tibble(rownames = "geneID")

Interactions_Upregulated_DEGsforGProfiler <- Interactions_DEGsforGProfiler.tibble %>%
  dplyr::filter (logFC > 0) %>%
  dplyr::select(geneID)

write.csv(Interactions_Upregulated_DEGsforGProfiler, "DEGs_increaseinDKOageing_decreaseinWTageing.txt", row.names = FALSE, quote = FALSE)

Interactions_Downregulated_DEGsforGProfiler <- Interactions_DEGsforGProfiler.tibble %>%
  dplyr::filter (logFC < 0) %>%
  dplyr::select(geneID) 

write.csv(Interactions_Downregulated_DEGsforGProfiler, "DEGs_decreaseinDKOageing_increaseinWTageing.txt", row.names = FALSE, quote = FALSE)

# for old WT minus Old DKO
DKOoldminusWTold_top10DEGs.tibble %>%
  gt() %>%
  cols_hide(columns = c(AveExpr, t, B)) %>%
  cols_label( geneID = "Gene",
              logFC = "Log2 (fold-change)", 
              P.Value = "p",
             adj.P.Val = "Adjusted p") %>%
  tab_header(
    title = md("**Differential gene expression Old DKO vs Old WT**"),
    subtitle = md("*Log2 (fold-change) shows Old DKO minus Old WT*")) %>%
  fmt_number(columns = logFC,n_sigfig = 3) %>%
  fmt_number(columns = P.Value, decimals = 4) %>%
  fmt_number(columns = adj.P.Val, decimals = 3)  %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()) %>%
  tab_style(
    style = cell_text(align = "center"),
    locations = cells_column_labels(P.Value))%>%
  tab_style(
    style = cell_text(align = "center"),
    locations = cells_column_labels(adj.P.Val)) %>%
  tab_style(
    style = cell_text(align = "center"),
    locations = cells_body(
      columns = adj.P.Val)) 

# for young WT minus young DKO

DKOYoungminusWTYoung_top10DEGs.tibble %>%
  gt() %>%
  cols_hide(columns = c(AveExpr, t, B)) %>%
  cols_label( geneID = "Gene",
              logFC = "Log2 (fold-change)", 
              P.Value = "p",
              adj.P.Val = "Adjusted p") %>%
  tab_header(
    title = md("**Differential gene expression: Young DKO vs Young WT**"),
    subtitle = md("*Log2 (fold-change) shows Young DKO minus Young WT*")) %>%
  fmt_number(columns = logFC,n_sigfig = 3) %>%
  fmt_number(columns = P.Value, decimals = 5) %>%
  fmt_number(columns = adj.P.Val, decimals = 3)  %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()) %>%
  tab_style(
    style = cell_text(align = "center"),
    locations = cells_column_labels(P.Value))%>%
  tab_style(
    style = cell_text(align = "center"),
    locations = cells_column_labels(adj.P.Val)) %>%
  tab_style(
    style = cell_text(align = "center"),
    locations = cells_body(
      columns = adj.P.Val))

# for interaction

interactions_DEGtable_sortedbyP <- Interaction_top10DEGs.tibble%>%
  gt() %>%
  cols_hide(columns = c(AveExpr, t, B)) %>%
  cols_label( geneID = "Gene",
              logFC = "Log2 (fold-change)", 
              P.Value = "p",
              adj.P.Val = "Adjusted p") %>%
  tab_header(
    title = md("**Differential gene expression: Interaction between Genotype and Age**"),
    subtitle = md("*Log2 (fold-change) shows difference in differences between gene expression for old and young DKO mice, and old and young WT mice*")) %>%
  fmt_number(columns = logFC,n_sigfig = 3) %>%
  fmt_number(columns = P.Value, decimals = 5) %>%
  fmt_number(columns = adj.P.Val, decimals = 3)  %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()) %>%
  tab_style(
    style = cell_text(align = "center"),
    locations = cells_column_labels(P.Value))%>%
  tab_style(
    style = cell_text(align = "center"),
    locations = cells_column_labels(adj.P.Val)) %>%
  tab_style(
    style = cell_text(align = "center"),
    locations = cells_body(
      columns = adj.P.Val)) %>%
tab_footnote(
  footnote = "statistically signficant at alpha level 0.05",
  locations = cells_body(
    columns = adj.P.Val,
    rows = c(1:10)))

# Volcano Plots ----

# plot for Old WT vs Old DKO
DKOoldminusWTold_volcanoplot <- ggplot(DKOoldminusWTold_allDEGs.tibble) +
  aes(y=-log10(adj.P.Val), x=logFC, text = paste("Symbol:", geneID)) +
  geom_point(size=0.5) +
  labs(title="Differential Gene Expression Analysis: Old DKO vs Old WT",
       subtitle = "Volcano Plot",
       caption=paste0("produced on ", Sys.time())) +
  xlab("Log2 (fold-change)") +
  ylab("-Log10 (adjusted P value)") +
  theme_bw() +
  theme(axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12),
        plot.title = element_text(size=16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(face = "italic", hjust = 0.5)
        )
  
# now plot for Young WT vs Young DKO
  
DKOYoungminusWTYoung_volcanoplot <- ggplot(DKOYoungminusWTYoung_allDEGs.tibble) +
  aes(y=-log10(adj.P.Val), x=logFC, text = paste("Symbol:", geneID)) +
  geom_point(size=0.5) +
  labs(title="Differential Gene Expression Analysis: Young DKO vs Young WT",
       subtitle = "Volcano Plot",
       caption=paste0("produced on ", Sys.time())) +
  xlab("Log2 (fold-change)") +
  ylab("-Log10 (adjusted P value)") +
  theme_bw() +
  theme(axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12),
        plot.title = element_text(size=16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(face = "italic", hjust = 0.5)
  )

# now plot for interaction

Interaction_volcanoplot <- ggplot(Interaction_allDEGs.tibble) +
  aes(y=-log10(adj.P.Val), x=logFC, text = paste("Symbol:", geneID)) +
  geom_point(size=0.5) +
  annotate("rect", xmin = 0, xmax = 12, ymin = -log10(0.05), ymax = 1.56, alpha=.2, fill="#BE684D") +
  annotate("rect", xmin = 0, xmax = -14, ymin = -log10(0.05), ymax = 1.56, alpha=.2, fill="#2C467A") +
  labs(title="Differential Gene Expression Analysis: DKO ageing vs WT ageing",
       subtitle = "Volcano Plot",
       caption=paste0("produced on ", Sys.time())) +
  xlab("Log2 (fold-change)") +
  ylab("-Log10 (adjusted P value)") +
  theme_bw() +
  theme(axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12),
        plot.title = element_text(size=16, face = "bold"),
        plot.subtitle = element_text(face = "italic")
  )



# decideTests to pull out the DEGs and make Venn Diagram ----
results.interaction <- decideTests(ebFit.interaction, method="global", adjust.method="BH", p = 0.1
                                   )

# take a look at what the results of decideTests looks like
head(results.interaction)
summary(results.interaction)
vennDiagram(results.interaction, include="both")

# retrieve expression data for your DEGs ----
head(v.DEGList.filtered.norm$E)
colnames(v.DEGList.filtered.norm$E) <- sampleLabels

interaction.diffGenes <- v.DEGList.filtered.norm$E[results.interaction[,3] !=0,]
head(interaction.diffGenes)
dim(interaction.diffGenes)

WToldminusDKOold.diffGenes <- v.DEGList.filtered.norm$E[results.interaction[,1] !=0,]
head(WToldminusDKOold.diffGenes)
dim(WToldminusDKOold.diffGenes)

WTyoungminusDKOyoung.diffGenes <- v.DEGList.filtered.norm$E[results.interaction[,2] !=0,]
head(WTyoungminusDKOyoung.diffGenes)
dim(WTyoungminusDKOyoung.diffGenes)

#convert your DEGs to a dataframe using as_tibble
interaction.diffGenes.tibble <- as_tibble(interaction.diffGenes, rownames = "geneID")

WToldminusDKOold.diffGenes.tibble <- as_tibble(WToldminusDKOold.diffGenes, rownames = "geneID")

WTyoungminusDKOyoung.diffGenes.tibble <- as_tibble(WTyoungminusDKOyoung.diffGenes, rownames = "geneID")


gt(interaction.diffGenes.tibble)


# Display DGEs with highest log-fold change for interaction -------------------------------

interaction.diffGenes.withaverages.tibble <- interaction.diffGenes.tibble %>% 
  mutate(WTold.AVG.column = (S10 + S11 + S12)/3,
         WTyoung.AVG.column = (S22 + S23 + S24)/3,
         DKOold.AVG.column = (S01 + S02 + S03)/3,
         DKOyoung.AVG.column = (S13 + S14 + S15)/3,
         Interaction.column = (WTold.AVG.column - WTyoung.AVG.column) - (DKOold.AVG.column - DKOyoung.AVG.column),
         WTaging.column = WTold.AVG.column - WTyoung.AVG.column,
         DKOaging.column = DKOold.AVG.column - DKOyoung.AVG.column,
         Redorblue = case_when(
           Interaction.column > 0 ~ "red",
           Interaction.column < 0 ~ "blue")
         ) %>% 
  mutate_if(is.numeric, round, 2)

# descending order
interaction.higherinWTageinggenes.bylogFC.tibble <- interaction.diffGenes.withaverages.tibble %>% 
  dplyr::filter(Interaction.column > 0) %>% 
  dplyr::arrange(desc(Interaction.column)) %>% 
   dplyr::select(geneID,
                 Interaction.column,
                 DKOaging.column,
                 WTaging.column,
                 Redorblue
                 )

# ascending order
interaction.higherinDKOageinggenes.bylogFC.tibble <- interaction.diffGenes.withaverages.tibble %>% 
  dplyr::filter(Interaction.column < 0) %>%    
  dplyr::arrange(Interaction.column) %>% 
                      dplyr::select(geneID,
                                    Interaction.column,
                                    DKOaging.column,
                                    WTaging.column,
                                    Redorblue)

# absolute values                   
interaction.Diffgenes.absolutevalues.bylogFC.tibble <- interaction.diffGenes.withaverages.tibble %>% 
  dplyr::arrange(desc(abs(Interaction.column))) %>% 
  dplyr::select(geneID,
                Interaction.column)  

write_tsv(interaction.higherinWTageinggenes.bylogFC.tibble,"Higher in WT ageing genes.txt")

write_tsv(interaction.higherinDKOageinggenes.bylogFC.tibble,"Higher in DKO ageing genes.txt")

write_tsv(interaction.Diffgenes.absolutevalues.bylogFC.tibble,"Interaction diff genes, absolute values.txt")

# scatterplots for signficant genes and log2FC ----------------------------

WTaging_vs_DKOaging_scatterplot <- ggplot(interaction.diffGenes.withaverages.tibble) +
  aes(y=DKOaging.column, x=WTaging.column, text = paste("Symbol:", geneID)) +
  geom_point(size=0.5, aes(colour = Redorblue)) +
  scale_color_hue(direction = -1) +
  labs(title="Differential Gene Expression Analysis: DKO ageing vs WT ageing",
       subtitle = "Scatterplot",
       caption=paste0("produced on ", Sys.time())) +
  xlab("Log2 (fold-change) during WT aging") +
  ylab("Log2 (fold-change) during DKO aging") +
  theme_bw() +
  theme(axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12),
        plot.title = element_text(size=16, face = "bold"),
        plot.subtitle = element_text(face = "italic")
  ) +
  guides(colour="none")


# now plot scatterplot of signficant genes of logFC in WT vs DKO aging, only inlcuding higher in WT aging

WTaging_vs_DKOaging_scatterplot_higherinWTageingonly <- ggplot(interaction.higherinWTageinggenes.bylogFC.tibble) +
  aes(y=DKOaging.column, x=WTaging.column, text = paste("Symbol:", geneID)) +
  geom_point(size=0.5, aes(colour = Redorblue)) +
  scale_color_hue(direction = -1) +
  labs(title="Differential Gene Expression Analysis: DKO ageing vs WT ageing",
       subtitle = "Scatterplot",
       caption=paste0("produced on ", Sys.time())) +
  xlab("Log2 (fold-change) during WT aging") +
  ylab("Log2 (fold-change) during DKO aging") +
  theme_bw() +
  theme(axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12),
        plot.title = element_text(size=16, face = "bold"),
        plot.subtitle = element_text(face = "italic")
  ) +
  guides(colour="none") +
  scale_x_continuous(position = "top") 

# now plot scatterplot of signficant genes of logFC in WT vs DKO aging, only including higher in DKO aging

WTaging_vs_DKOaging_scatterplot_higherinDKOageingonly <- ggplot(interaction.higherinDKOageinggenes.bylogFC.tibble) +
  aes(y=DKOaging.column, x=WTaging.column, text = paste("Symbol:", geneID)) +
  geom_point(size=0.5, color = "#619CFF") +
  labs(title="Differential Gene Expression Analysis: DKO ageing vs WT ageing",
       subtitle = "Scatterplot",
       caption=paste0("produced on ", Sys.time())) +
  xlab("Log2 (fold-change) during WT aging") +
  ylab("Log2 (fold-change) during DKO aging") +
  theme_bw() +
  theme(axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12),
        plot.title = element_text(size=16, face = "bold"),
        plot.subtitle = element_text(face = "italic")
  ) +
  guides(colour="none") +
  scale_y_continuous(position = "right") +
  scale_color_hue(direction = 1)

# Extract data for bar charts for top 10 DEGs -----------------------------
barchartdata <- interaction.diffGenes.tibble %>%
  dplyr::filter(geneID=="Cops7b" | geneID=="Fam171a2" | 
                  geneID=="Bai2" | geneID=="Parp14" |
                  geneID=="Ensa" | geneID=="Hspa1a" |
                  geneID=="Rasl11b" | geneID=="Rasal1" |
                  geneID=="Junb" | geneID=="Ctdsp1") %>%
  dplyr::select(geneID, S01, S02, S03, S04, S05, S06, S07, S08, S09, S10, S11, S12, S13, S14, S15, S16, S17, S18, S19, S20, S21, S22, S23, S24)

library(openxlsx)

openxlsx::write.xlsx(barchartdata, "barchartdata.xlsx")