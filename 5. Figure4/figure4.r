library(dplyr)
library(ggsci)
library(ggplot2)
library(ggpubr)

group <- read.csv("00.rawdata/BCH_MP_metadata_88.csv",check.names = F,row.names = 1)
samplesNormal <- group %>% dplyr::filter(Group == 'Normal') %>% rownames()
samplesMucES <- group %>% dplyr::filter(Group == 'MucES-HS') %>% rownames()
samplesMucHS <- group %>% dplyr::filter(Group == 'MucHS') %>% rownames()

############################## diff KOs ####################
data <- read.table("00.rawdata/KO.screen.txt",sep="\t",header=T,row.names=1,check.names = F)
wilcox.01.p <- apply(data,1,function(x) {if(sum(x[c(samplesNormal,samplesMucES)]) == 0){return('NA')}else{test <- wilcox.test(x[samplesNormal],x[samplesMucES],conf.int = T);p <- test$p.value;return(p);}}) %>% as.numeric()
MucES.HS_Mean <- apply(data,1,function(x) {mean(x[samplesMucES])}) %>% as.numeric()
Normal_Mean <-  apply(data,1,function(x) {mean(x[samplesNormal])}) %>% as.numeric()
MucHS_Mean <- apply(data,1,function(x) {mean(x[samplesMucHS])}) %>% as.numeric()
MucES.HS_Median <- apply(data,1,function(x) {median(x[samplesMucES])}) %>% as.numeric()
Normal_Median <-  apply(data,1,function(x) {median(x[samplesNormal])}) %>% as.numeric()
result <- data.frame(pValue=wilcox.01.p,MucES.HS_Mean=MucES.HS_Mean,Normal_Mean,MucES.HS_Median=MucES.HS_Median,Normal_Median=Normal_Median,MucHS_Mean=MucHS_Mean)
rownames(result) <- rownames(data)

KOsign1 <- result %>% dplyr::filter(pValue < 0.05) %>% dplyr::mutate(log2FC=ifelse(MucES.HS_Mean > Normal_Mean,log2(MucES.HS_Mean/Normal_Mean),-log2(Normal_Mean/MucES.HS_Mean))) %>% dplyr::filter(abs(log2FC) > 1) %>% dplyr::filter(MucES.HS_Median > 0 | Normal_Median > 0 ) %>% rownames()
KO2pathway <- read.table("00.rawdata/signKO_pathway.txt",sep = "\t",header = T,row.names = 1)
rownames(KO2pathway) <- gsub(" .*","",rownames(KO2pathway))
heatmapIN <- result[KOsign1,c("MucES.HS_Mean","MucHS_Mean","Normal_Mean")]
rownames(heatmapIN) <- gsub(":.*","",rownames(heatmapIN))
pathCol <- c(pal_bmj()(8),pal_d3()(8)[c(3,5,6,7)],"#C6DBEF","#EE0000","#7F7F7F")
names(pathCol) <- unique(KO2pathway$Pathway)
pheatmap::pheatmap(heatmapIN[rownames(KO2pathway),],annotation_row=KO2pathway,color = colorRampPalette(colors = c("navy","white","firebrick3"))(100),annotation_colors = list(Pathway=pathCol),border=NA,show_colnames=T,show_rownames = F,width=4,heigh=4,fontsize=6,cluster_cols=F,cluster_rows = F,filename = "KOs.heatmap.pdf",scale = "row",angle_col = "45")

############################## diff VFs ####################
VFtable <- read.table("00.rawdata/VF.screen.txt",sep="\t",header=T,row.names=1,check.names = F)  %>% t() %>% as.data.frame()
VFtable$Group <- as.vector(group[rownames(VFtable),'Group'])

VFtable %>% reshape2::melt(by='Group') %>% ggplot(aes(Group,value,fill=Group)) + geom_boxplot() + facet_wrap('variable',nrow = 3,scales = 'free') + geom_signif(comparisons = list(c('Normal','MucES-HS'),c('Normal','MucHS'),c('MucES-HS','MucHS')), color="black",step_increase = 0.08, map_signif_level=function(p) if(p < 0.01){'+'}else{if(p < 0.05){'*'}else{sprintf("%.2f",p)}},test = wilcox.test) 

VFtable %>% reshape2::melt(by='Group') %>% 
  dplyr::filter(variable %in% c('VFG002035(gb|WP_010874803)','VFG002029(gb|WP_010874666)','VFG002034(gb|WP_010874808)','VFG002031(gb|WP_010874665)','VFG002036(gb|WP_010874809)','VFG046465(gb|WP_003028672)')) %>% 
  ggplot(aes(variable,value,fill=Group)) + geom_boxplot(outlier.size = 0.5)  + coord_cartesian(ylim=c(0,0.2))  + xlab("") + ylab("Abundance") + theme_bw() + theme(axis.text.x = element_text(angle = 45,hjust = 1))

