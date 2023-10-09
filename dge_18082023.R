#DGE analysis for SB-positive and negative-cell lines
BiocManager::install('tximport')
library("tximport")
library("RColorBrewer")
library("DESeq2")
library("readr")
library("gplots")
library("ggrepel")
library("ggplot2")
library("circlize")
library("ComplexHeatmap")

setwd("E:/MSR BLY/lab work/202190 Autophagy")
#read model information
model<- read.csv("Model.csv")
profiles <- read.csv("OmicsProfiles.csv")
profiles_rna <-profiles[which(profiles$Datatype=="rna"),]

#make a cell line and results table
sampleDat <- data.frame(sample=c("A549", "AGS", "CACO2", "DU145", "HCT116",
                                 "HEK293T", "HT29", "HUH7", "MCF10A", 
                                 "MCF7", "PC3", "SHSY5Y", "SW480", "WM1617",
                                 "WM793", "NCIH460", "MDAMB231", "U87MG", 
                                 "HACAT", "LN229"),
                        status=as.factor(c("pos", "pos", "pos", "pos", "pos", "neg",
                                 "pos", "neg", "pos", "neg", "neg", "neg", "pos",
                                 "neg", "neg", "pos", "pos", "pos", "neg", "pos")))

#read expected counts
expected_counts <- read_csv("OmicsExpressionGenesExpectedCountProfile.csv")
rsem<- as.data.frame(expected_counts)

#matching cell line names
model_idx <- profiles$ModelID[match(rsem[, 1], profiles$ProfileID)]
rsem$cell_line<- model$StrippedCellLineName[match(model_idx, model$ModelID)]

#match cell lines
sampleDat$idxs <- match(sampleDat$sample, rsem$cell_line)
sampleDat1 <- na.omit(sampleDat)

rsem_sb90<- t.data.frame(rsem[sampleDat1$idxs,])
colnames(rsem_sb90)<-rsem_sb90[54345,]

#remove the profileID and cell line name and then round off the data
#removing genes if they have zero counts in 5 or more samples
rsem_sb90<- rsem_sb90[-c(1, 54345),]
rsem_sb90R <- apply(rsem_sb90, 2, function(x){
  round(as.numeric(x), digits=0)
})

genes_keep<- sapply(1:nrow(rsem_sb90R), function(i){
  if (length(which(rsem_sb90R[i,]==0))<10){
    return(i)
  } else
  return(NA)
})

rsem_sb90R<- cbind.data.frame(rownames(rsem_sb90), rsem_sb90R)
rsem_sb90R<-rsem_sb90R[na.omit(genes_keep),]

#DESeq2
dds <- DESeqDataSetFromMatrix(countData=rsem_sb90R, 
                              colData=sampleDat1[,1:2], 
                              design=~status, tidy = TRUE)
dds<- DESeq(dds)
#results
res <- results(dds)
head(results(dds, tidy=TRUE))

#summary
summary(res)
res <- res[order(res$padj),]
head(res)
write.table(res, file="18082023_dge_sb90_filter=10.txt", sep="\t")

#vst
vsdata <- vst(dds, blind=FALSE)
plotPCA(vsdata, intgroup="status")

genes<- sapply(res@rownames, function(x){strsplit(x, split=" ")[[1]][1]})
write.table(genes, file="18072023_genes.txt", row.names = FALSE)

vst_matrix<- vsdata@assays@data@listData[[1]]

genes_pos<- vector()
all_genes<- c(res@rownames[1:50])

for (gene in all_genes){
  pos<-which(rownames(vst_matrix)==gene)
  genes_pos <- c(genes_pos, pos)
}
genes_matrix<- vst_matrix[genes_pos,]

heatmap(genes_matrix[1:50,])
#colors for heatmap
col_fun = colorRamp2(c(2.5, 10, 15), c("#00549f", "#edd1c5", "#ac0e2b"))
col_vac<- sapply(sampleDat1$status, function(x){
  if(x=="pos"){status="#990000"}
  else if(x=="neg"){status="#406931"}
})

col_heatmap <- c("#00549f", "#edd1c5", "#ac0e2b") 
pal <- colorRampPalette(col_heatmap)(100)
rownames(genes_matrix) <- sapply(rownames(genes_matrix), function(x){unlist(strsplit(x, split=" \\("))[1]})
heatmap(genes_matrix[1:50,], col=pal, ColSideColors=col_vac, cexRow = 0.9, cexCol = 1.15)

targets <- c("MAPK11", "MAPK14", "MAPK13", "MAPK12", "WNK1", "WNK2", "WNK3", "WNK4",
             "RIPK2", "GAK", "BRAF", "RAF1", "GSK3B", "LCK", "PDK1", "TGFBR2", "CSNK1A1",
             "RIPK3", "AKT1", "BECN1")

rownames(vst_matrix)<- sapply(rownames(vst_matrix), function(x){unlist(strsplit(x, split=" \\("))[1]})
heatmap(vst_matrix[match(targets, rownames(vst_matrix)),], col=pal, ColSideColors=col_vac, cexRow = 0.9, cexCol = 1.15)
legend(x="topleft", legend=c("min", "ave", "max"), fill=colorRampPalette(col_heatmap)(3), trace=TRUE, title="Expression")

#a better volcano plot!
de <- read.delim("18082023_dge_sb90_filter=10.txt", row.names = 1)
de$gene<-sapply(rownames(de), function(x){unlist(strsplit(x, split=" \\("))[1]})
de$diffexpressed <- "NO"
# if log2Foldchange > 2 and pvalue < 0.05, set as "UP" 
de$diffexpressed[de$log2FoldChange > 2 & de$padj < 0.01] <- "UP"
# if log2Foldchange < -2 and pvalue < 0.01, set as "DOWN"
de$diffexpressed[de$log2FoldChange < -2 & de$padj < 0.01] <- "DOWN"
de$delabel <- NA
de$delabel[de$diffexpressed != "NO"] <- de$gene[which(de$diffexpressed != "NO")]
de_labels <- match(c("AEBP1", "RELN", "ANPEP", "LAMP5", "DLK1", "UBE2QL1",
                     "IGFL2", "GRIN2B", "MYOM3", "NLRP2"), de$delabel)

#adding colors for DEGs
mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")
de_lab<- de[de_labels,]
ggplot(data=de, aes(x=log2FoldChange, y=-log10(padj), 
                     col=diffexpressed)) +
  geom_point(alpha=0.5, size=5)+
  geom_vline(xintercept=c(-2, 2), col="red", linetype=2) +
  geom_hline(yintercept=-log10(0.01), col="red", linetype=2)+
  theme(plot.background = element_blank(),
        axis.title = element_text(face="bold"),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        legend.position = "bottom",
        panel.grid.major = element_line(colour = "azure3"))+
  geom_label_repel(data=de_lab, label=de_lab$delabel, show.legend = FALSE)+
  labs(x="log2 Fold Change",
       y="-log10(adjusted p-value)",
       title = "Vacuole (+) versus Vacuole (-) cell lines")+
  scale_colour_manual(name="Differential Expression", values = mycolors)

#boxplot
status_vac<- sapply(sampleDat1$status, function(x){
  if(x=="pos"){status="positive"}
  else if(x=="neg"){status="negative"}
})
genes<- c("NLRP2", "GRIN2B", "IGFL2", "MYOM3", "RELN", "LAMP5", "AEBP1", "UBE2QL1")
genes_plot<- genes_matrix[match(genes, row.names(genes_matrix)),]
colnames(genes_plot)<- status_vac
cols=c("positive" = "red3", "negative" = "forestgreen")


plot_dat<- reshape2::melt(genes_plot)
plot_exp <- ggplot(data=plot_dat, aes(x=Var2, y=value, fill=Var2, colour=Var2))+
  geom_boxplot(alpha=0.7)+geom_jitter(show.legend = FALSE)+
  theme(plot.background = element_blank(),
        axis.title = element_text(face="bold"),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        legend.position="bottom", 
        legend.title=element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.x = element_text(face="bold"),
        panel.grid.minor = element_line(colour = "azure3")
        )+
  labs(y="Gene Expression",
       x="SB202190-induced Vacuolation Status")+
  scale_colour_manual(values=cols)+
  scale_fill_manual(values=cols)+facet_wrap(~ Var1, ncol=4, nrow=2)
plot(plot_exp)

print(sessionInfo())

