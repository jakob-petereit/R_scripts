# libraries ----
library(DESeq2)
library(data.table)
library(tidyverse)
library(tximport)
library(AnnotationHub)
library(ggrepel)
library(ensembldb)
library(Biostrings)




# load  db and quant from salmon ----
#arabidopsis reference transcripts (tx2gene)
## Load the annotation resource.

ah <- AnnotationHub()
ahDb <- query(ah, pattern = c("thaliana"))
ahEdb <- ahDb[[4]]
k <- keys(ahEdb, keytype = "TXNAME")
tx2gene <- select(ahEdb, k, "GENEID", "TXNAME")


##import reads after running Salmon
# dir is path/to/dir
#make a data/quant folder in  working directory

dir <- ('data/quants/')
samples <- list.files(dir)
files <- file.path(dir, samples, "quant.sf")
names(files) <- c('clp1_1','clp1_2','clp1_3',
                  'clp2_1','clp2_2','clp2_3',
                  'col0_1','col0_2','col0_3')
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene)
head(txi.salmon$counts)

# make deseq dataset ----

coldata <- tibble(rowname=names(files)) %>% 
  mutate(genotype=sapply(strsplit(rowname,'_'),'[',1)) %>% 
  column_to_rownames()

dds <- DESeqDataSetFromTximport(txi.salmon,
                                   colData = coldata,
                                   design = ~ genotype)




# Prefiltering----

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

##factor leveling, setting col0 as default

dds$genotype <- factor(dds$genotype, levels = c('col0','clp1','clp2'))

# Differential expression analysis----
dds <- DESeq(dds)

res <- results(dds)
res


#stand alone sections
# DESEQ workflow from guide ----
#Log fold change shrinkage for visualization and ranking
resultsNames(dds)
resLFC <- lfcShrink(dds, coef="genotype_clp2_vs_col0", type="apeglm")
resLFC
#p-values and adjusted p-values
resOrdered <- res[order(res$pvalue),]

summary(res)
#How many adjusted p-values were less than 0.1?
sum(res$padj < 0.1, na.rm=TRUE)

#same for 0.05
res05 <- results(dds, alpha=0.05)
summary(res05)
sum(res$padj < 0.05, na.rm=TRUE)

##Exploring and exporting results
#MA-plot
#raw
plotMA(res, ylim=c(-2,2))
#shrunken
plotMA(resLFC, ylim=c(-2,2))

#Alternative shrinkage estimators
resNorm <- lfcShrink(dds, coef=2, type="normal")
resAsh <- lfcShrink(dds, coef=2, type="ashr")
par(mfrow=c(1,3), mar=c(4,4,2,1))
xlim <- c(1,1e5); ylim <- c(-3,3)
plotMA(resLFC, xlim=xlim, ylim=ylim, main="apeglm")
plotMA(resNorm, xlim=xlim, ylim=ylim, main="normal")
plotMA(resAsh, xlim=xlim, ylim=ylim, main="ashr")

#Plot counts

d <- plotCounts(dds, gene=toupper('AT5G63510'), intgroup="genotype", 
                returnData=TRUE)
# ggplot(d, aes(x=genotype, y=count)) + 
#   geom_point(position=position_jitter(w=0.1,h=0)) + 
#   scale_y_log10()


ggplot(d, aes(genotype, count, color=genotype)) + 
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", size = 0.5)+
  geom_point(size=3,color='black')+
  expand_limits(y=0)+
  labs(title='Raw Count comparison')
  
#Exporting results to CSV files

#clp2 vs col-0
write.csv(as.data.frame(resOrdered), 
          file="condition_treated_results.csv")




#QC
#works on count transformation
#shifted log transformation (Log2(n+1))
ntd <- normTransform(dds)
#variance stabiliszing transformation
vsd <- vst(dds, blind=FALSE)
#regularized log transforamtion
rld <- rlog(dds, blind=FALSE)

#meanSD plots
#NTD
library("vsn")
meanSdPlot(assay(ntd))
#VSD
meanSdPlot(assay(vsd))
#rld
meanSdPlot(assay(rld))

#Heatmap of the count matrix
library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("genotype")])
rownames(df) <- colnames(dds)

#ntd
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=T,
         cluster_cols=FALSE, annotation_col=df)
#vsd
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=T,
         cluster_cols=FALSE, annotation_col=df)
#rld
pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=T,
         cluster_cols=FALSE, annotation_col=df)

#Heatmap of the sample-to-sample distances
sampleDists <- dist(t(assay(vsd)))

library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- vsd$genotype
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

#Principal component plot of the samples
pcaData <- plotPCA(vsd, intgroup=c("genotype"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pcaData$genotype <- factor(pcaData$genotype, levels=c('col0','clp1','clp2'))
ggplot(pcaData, aes(PC1, PC2, color=genotype)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()+
  scale_colour_manual(values=c('#339900','#0066cc','#990033'))+
  labs(title='PCA RNA Seq Col-0 vs clp1/cl2')+
  theme(title = element_text(size=18,face='bold'),axis.text.x = element_text(size=16,face='bold'),axis.title = element_text(size=16,face='bold'),axis.text.y = element_text(size=12),
        strip.text = element_text(size=12,face='bold'),legend.text=element_text(size=14,face='bold'),legend.key.size = unit(2,'cm'))
  





# data Visualisation----

#Overlap volcano====

## custom mutant overlap volcano

#work with result files for clp1 and clp2 by using different resultsnames

#clp1
resLFC_clp1 <- lfcShrink(dds, coef="genotype_clp1_vs_col0", type="apeglm") %>% 
  as.data.frame() %>% 
  rownames_to_column() %>%  
  as_tibble() %>% 
  mutate(genotype='clp1')
#clp2
resLFC_clp2 <- lfcShrink(dds, coef="genotype_clp2_vs_col0", type="apeglm") %>% 
  as.data.frame() %>% 
  rownames_to_column() %>%  
  as_tibble() %>% 
  mutate(genotype='clp1')


#combine and make long format

reslfc <- bind_rows(resLFC_clp1,resLFC_clp2)

#read from file
clpvswt <- fread('clp1_2vs_wt_results.csv')
suba <- fread('data/1-7-19_SUBA4_display.csv') %>% 
  dplyr::select(1,3,40) %>% 
  mutate(locus=substr(locus,1,9)) %>% 
  distinct(locus, .keep_all = T)

#filter for both mutants under 0.1pvalue, and use lowest FC and highest pvalue
data_volcano <-  
  reslfc %>% 
  group_by(rowname) %>% 
  mutate(colour_p=ifelse(max(padj) <= 0.1,'red','blue')) %>% 
  mutate(min_ratio=min(abs(log2FoldChange)),
         max_p=max(padj)) %>% 
  dplyr::filter(abs(log2FoldChange)==min_ratio) %>% 
  dplyr::select(-min_ratio) %>% 
  mutate(colour_r=ifelse(log2FoldChange <=-0.4 | log2FoldChange >= 0.4,'red','blue')) %>% 
  mutate(sig=ifelse(colour_p=='blue'|colour_r=='blue','non_sig','sig')) %>% 
  distinct(rowname, .keep_all = T)

#add good annotation
#write significant genes to clipboard
#writeClipboard(dplyr::filter(data_volcano, sig=='sig')$rowname)
#go to tair > bulkd data retrieval > gene descriptions > download as text
#read sig_gene_descriptions.txt 
gene_desc <- fread('data/sig_gene_desc.txt') %>% 
  mutate(desc=ifelse(`Gene Model Type`!='',`Gene Model Type`,`Gene Model Name`)) %>% 
  dplyr::select(-c(2:7))
#adjust annotation by hand in text editor
#write(commented our for safety)
#write.table(gene_desc,'data/sig_gene_desc2.txt')
#read
gene_desc <- fread('data/sig_gene_desc2.txt') %>% dplyr::select(-1)



#merge annnotation into volcano data

data_volcano <- data_volcano %>% 
  left_join(gene_desc, by=c('rowname'='gene')) 


#levels
data_volcano$sig <- factor(data_volcano$sig, levels=c('sig','non_sig'))

## format data for plot
#change-log10 pvalue  high numbers to a maximum of 7.5
#change symbol of  >7.5 to op arrow

data_volcano <- data_volcano %>% 
  mutate(max_p_adj=ifelse(rowname=='AT5G23140',1e-65,max_p),
         pch=ifelse(rowname=='AT5G23140',17,16),
         size=ifelse(pch==17,3,2),
         ratio_adj=ifelse(log2FoldChange > 5,5,ifelse(log2FoldChange < -5,-5,log2FoldChange))) 

#remove description for hypothetical proteins, transposable elements, retrotransposons and pseudogenes, change their symbol to an empty point
data_volcano <- data_volcano %>% 
  mutate(pch=as.numeric(pch),
         desc2=ifelse(desc =='h'|
                      desc=='TE'|
                      desc=='rTE'|
                      str_detect(desc,'pseudo')==T,'',desc),
         pch=ifelse(desc2=='',1,pch))

data_volcano <- data_volcano %>% 
  mutate(pch=ifelse(is.na(pch),16,pch))

#Write_data
write.csv(data_volcano,'data/data_voclano_rnaseq.csv')

p <- ggplot(data_volcano, aes(x=ratio_adj,y=-log(max_p_adj,10),col=sig))+
  geom_point(pch=data_volcano$pch,size=data_volcano$size,alpha=0.75)+
  geom_text_repel(data=dplyr::filter(data_volcano,sig=='sig'),aes(label=desc2),col='black',size=2.5, fontface='bold')+
  scale_colour_manual(values=c('#990000','#99ccff'))+
  theme(legend.position = 'none', axis.title = element_text(face='bold',size = 18),
        axis.text = element_text(face='bold',size = 16), strip.text = element_text(face='bold',size=18),
        title = element_text(face='bold',size=18))+
  labs(title='Clp vs col-O', x='log2 fold change', y='-log10 p-value')
p

ggsave(filename = paste0('volcano','.png'),path = 'images',device = 'png',dpi=1080,plot = p)




#complex I (incomplete) ####

#select Genes
cI_anno <- fread('data/complex I annotation.csv', header = T) %>% dplyr::select(1:3) %>% as_tibble()

cI_anno <- cI_anno %>% 
  mutate(ACC=gsub(pattern = ', ',replacement = '/',x = ACC)) %>% 
  mutate(ACC1=sapply(strsplit(ACC,'/'),'[',1),
         ACC2=sapply(strsplit(ACC,'/'),'[',2),
         ACC3=sapply(strsplit(ACC,'/'),'[',3)) %>% 
  dplyr::select(-ACC) %>% 
  gather(nr,ACC,ACC1:ACC3) %>% 
  dplyr::select(-nr) %>% 
  na.omit() %>% 
  mutate(ACC=toupper(substr(ACC,1,9))) %>% 
  dplyr::select(-group)
#Plot counts

d <- plotCounts(dds, gene=toupper('AT1G49630'), intgroup="genotype", 
                returnData=TRUE)
# ggplot(d, aes(x=genotype, y=count)) + 
#   geom_point(position=position_jitter(w=0.1,h=0)) + 
#   scale_y_log10()


ggplot(d, aes(genotype, count, color=genotype)) + 
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", size = 0.5)+
  geom_point(size=3,color='black')+
  expand_limits(y=0)+
  labs(title='Raw Count comparison')
#CLPp and control Genes for paper figure 1 ####
#Gene selection
#clpp
clpp <- 'AT5G23140'
#most stabel house keeper overall
#https://link.springer.com/article/10.1007/s13562-017-0403-0
tub6 <- 'AT5G12250'
#mito housekeeper vdac1
vdac1 <- 'AT3G01280'

#selction

sel <- c(clpp,tub6,vdac1)

#Plot counts
#custom selection, Plotcounts() doesn't work on multiple Genes
subset <- assay(dds[sel]) %>% 
  as.data.frame() %>% 
  rownames_to_column(var='gene') %>% 
  as_tibble() %>% 
  gather(sample,count,2:10) %>% 
  mutate(genotype=sapply(strsplit(sample,'_'),'[',1),
         genotype=ifelse(genotype=='col0','WT',genotype))
#description
desc_sel <- tibble(gene=c('AT5G23140','AT5G12250','AT3G01280'),desc=c('CLPP2','Tubulin6','ATVDAC1'))


#join desc onto pvalue

#join desc and pvalue to subset
subset <- subset %>% 
  left_join(desc_sel)

#make pvalue frame and get descriptions
#resultsNames(dds)
resclp1 <- results(dds[sel], name = 'genotype_clp1_vs_col0',tidy = T) %>% as_tibble() %>% 
  dplyr::select(row,padj) %>% 
  dplyr::rename(gene=row) %>% 
  mutate(genotype='clp1')
resclp2 <- results(dds[sel], name = 'genotype_clp2_vs_col0',tidy = T) %>% as_tibble() %>% 
  dplyr::select(row,padj) %>% 
  dplyr::rename(gene=row) %>% 
  mutate(genotype='clp2')

res <- bind_rows(resclp1,resclp2) %>% 
  left_join(desc_sel) %>% 
  dplyr::select(-gene) %>% 
  mutate(sig_level=ifelse(padj > 0.01,'',
                                 ifelse(padj <= 0.01 & padj > 0.001,'*',
                                        ifelse(padj <= 0.001,'*\n*','failsave'))))
#transform into thousands 'K'
subset <- subset %>% 
  mutate(count=count/1000)


#add median as ref y axis point 
y_ref <- subset %>% group_by(genotype,desc) %>% 
  summarise(median=median(count),max=max(count))
res <- res %>% 
  left_join(y_ref)

#levels
subset$genotype <- factor(subset$genotype, levels = c('WT','clp1','clp2'))
subset$desc <- factor(subset$desc, levels = c('CLPP2','Tubulin6','ATVDAC1'))

res$genotype <- factor(res$genotype, levels = c('WT','clp1','clp2'))
res$desc <- factor(res$desc, levels = c('CLPP2','Tubulin6','ATVDAC1'))



#plot
p <- ggplot(subset, aes(genotype, count,bg=genotype)) + 
  facet_wrap(~desc,scales = 'free_y')+
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar",col='black', size = 0.3)+
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "bar",col='black', size = 0.15,alpha= 0.6)+
  geom_point(pch = 21,size=2,color='black',alpha=0.5)+
  geom_text(data=res,aes(genotype,c(0.4,1,1,0.3,1,1),label=sig_level),size=6, lineheight = 0.25)+
  expand_limits(y=0)+
  scale_colour_manual(values=c('#339900','#3399cc','#3366cc'))+
  scale_fill_manual(values=c('#339900','#3399cc','#3366cc'))+
  labs(title='Total transcript abundance',y='Normalised counts [k]')+
  theme(axis.title.x = element_blank(),legend.position = 'none',
        axis.text.x = element_text(face=c('plain','italic','italic'),size=8, angle = 30),
        axis.title.y = element_text(face='bold',size='8'),
        axis.text.y = element_text(face='bold',size=8),
        strip.text = element_text(face='bold',size=8),
        title=element_text(size=10))

#save
ggsave('RNA_KO_figure1.pdf',device = 'pdf',dpi=1080,plot = p,height = 6.52,width = 6,units = 'cm')




#proteases

##GO term analysis


#mito encoded ####
#Gene selection


#selction

sel <- str_detect(rownames(dds),'ATM')==T

#Plot counts
#custom selection, Plotcounts() doesn't work on multiple Genes
subset <- assay(dds[sel]) %>% 
  as.data.frame() %>% 
  rownames_to_column(var='gene') %>% 
  as_tibble() %>% 
  gather(sample,count,2:10) %>% 
  mutate(genotype=sapply(strsplit(sample,'_'),'[',1),
         genotype=ifelse(genotype=='col0','WT',genotype))
#description
desc_sel <- fread('data/desc_atm.csv')


#join desc onto pvalue

#join desc and pvalue to subset
subset <- subset %>% 
  left_join(desc_sel, by=c('gene'='AGI'))

#make pvalue frame and get descriptions
#resultsNames(dds)
resclp1 <- results(dds[sel], name = 'genotype_clp1_vs_col0',tidy = T) %>% as_tibble() %>% 
  dplyr::select(row,padj) %>% 
  dplyr::rename(gene=row) %>% 
  mutate(genotype='clp1')
resclp2 <- results(dds[sel], name = 'genotype_clp2_vs_col0',tidy = T) %>% as_tibble() %>% 
  dplyr::select(row,padj) %>% 
  dplyr::rename(gene=row) %>% 
  mutate(genotype='clp2')

res <- bind_rows(resclp1,resclp2) %>% 
  left_join(desc_sel,by=c('gene'='AGI')) %>% 
  dplyr::select(-gene) %>% 
  mutate(sig_level=ifelse(padj > 0.01,'',
                          ifelse(padj <= 0.01 & padj > 0.001,'*',
                                 ifelse(padj <= 0.001,'*\n*','failsave'))))
#transform into thousands 'K'
subset <- subset %>% 
  mutate(count=count/1000)


#add median as ref y axis point 
y_ref <- subset %>% group_by(genotype,desc) %>% 
  summarise(median=median(count),max=max(count))
res <- res %>% 
  left_join(y_ref)


#atg on facet
subset <- subset %>% 
  mutate(strip=paste0(gene,'\n',desc))

#levels
subset$genotype <- factor(subset$genotype, levels = c('WT','clp1','clp2'))

res$genotype <- factor(res$genotype, levels = c('WT','clp1','clp2'))



#plot
p <- ggplot(dplyr::filter(subset,str_detect(desc,'ORF',negate = T)), aes(genotype, count,bg=genotype)) + 
  facet_wrap(~strip,scales = 'free_y')+
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar",col='black', size = 0.3)+
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "bar",col='black', size = 0.15,alpha= 0.6)+
  geom_point(pch = 21,size=2,color='black',alpha=0.5)+
  expand_limits(y=0)+
  scale_colour_manual(values=c('#339900','#3399cc','#3366cc'))+
  scale_fill_manual(values=c('#339900','#3399cc','#3366cc'))+
  labs(title='Total transcript abundance',y='Normalised counts [k]')+
  theme(axis.title.x = element_blank(),legend.position = 'none',
        axis.text.x = element_text(face=c('plain','italic','italic'),size=8, angle = 30),
        axis.title.y = element_text(face='bold',size='8'),
        axis.text.y = element_text(face='bold',size=8),
        strip.text = element_text(face='bold',size=12),
        title=element_text(size=10))

#save
ggsave('atm_no_ORF.pdf',device = 'pdf',dpi=1080,plot = p,height = 10,width = 16,units = 'in')




# results table with annotation, location and chromosome position
#tables
#clp1
resLFC_clp1 <- lfcShrink(dds, coef="genotype_clp1_vs_col0", type="apeglm") %>% 
  as.data.frame() %>% 
  rownames_to_column() %>%  
  as_tibble() %>% 
  mutate(genotype='clp1')
#clp2
resLFC_clp2 <- lfcShrink(dds, coef="genotype_clp2_vs_col0", type="apeglm") %>% 
  as.data.frame() %>% 
  rownames_to_column() %>%  
  as_tibble() %>% 
  mutate(genotype='clp2')


#combine and make long format

reslfc <- bind_rows(resLFC_clp1,resLFC_clp2)

# Table -------------------------------------------------------------------

#work with result files for clp1 and clp2 by using different resultsnames

#clp1
resLFC_clp1 <- lfcShrink(dds, coef="genotype_clp1_vs_col0", type="apeglm") %>% 
  as.data.frame() %>% 
  rownames_to_column() %>%  
  as_tibble() %>% 
  mutate(genotype='clp1')
#clp2
resLFC_clp2 <- lfcShrink(dds, coef="genotype_clp2_vs_col0", type="apeglm") %>% 
  as.data.frame() %>% 
  rownames_to_column() %>%  
  as_tibble() %>% 
  mutate(genotype='clp2')


#combine and wide format

reslfc <- bind_rows(resLFC_clp1,resLFC_clp2) %>% 
  dplyr::select(rowname,log2FoldChange,padj,genotype) %>% 
  mutate(comb=paste(log2FoldChange,padj,sep='_')) %>% 
  dplyr::select(-log2FoldChange,-padj) %>% 
  spread(genotype,comb) %>% 
  mutate(clp1_lognfc=sapply(strsplit(clp1,'_'),'[',1),
         clp1_padj=sapply(strsplit(clp1,'_'),'[',2),
         clp2_lognfc=sapply(strsplit(clp2,'_'),'[',1),
         clp2_padj=sapply(strsplit(clp2,'_'),'[',2)) %>% 
  dplyr::select(-clp1,-clp2)

#add suba
#extra descriptions for the ones that got away (to suba description)
extra <- fread('data/extra desc for table.txt') %>% 
  dplyr::select(1,3)


suba <- fread('data/1-7-19_SUBA4_display.csv') %>% 
  dplyr::select(1,3,40) %>% 
  mutate(locus = substr(locus,1,9)) %>% 
  left_join(extra, by=c('locus'='V1')) %>% 
  distinct(.keep_all = T) %>% 
  mutate(description=ifelse(description==';',`Gene Model Name`,description)) %>% 
  dplyr::select(-`Gene Model Name`) %>% 
  group_by(locus) %>% 
  mutate(location_consensus=paste0(location_consensus,collapse = ',')) %>% 
  ungroup() %>% 
  distinct(locus,.keep_all = T)%>%
  mutate(location_consensus=gsub('plasma membrane','plasma_membrane',location_consensus),
         location_consensus=gsub('endoplasmic reticulum','endoplasmic_reticulum',location_consensus)) %>% 
  separate_rows(location_consensus) %>% 
  distinct (locus,location_consensus,.keep_all = T) %>% 
  group_by(locus,description) %>% 
  summarise(location_consensus = paste(location_consensus, collapse=","))



reslfc <- reslfc %>% 
  left_join(suba, by=c('rowname'='locus'))

#extract position from tair db fasta
pos <- readAAStringSet('data/Arabidopsis_thaliana.TAIR10.pep.all.fa') %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::select(1) %>% 
  mutate(AGI=sapply(strsplit(rowname,'[.]'),'[',1),
         left=sapply(strsplit(rowname,':'),'[',5),
         right=sapply(strsplit(rowname,':'),'[',6),
         direction=substr(sapply(strsplit(rowname,':'),'[',7),1,2),
         direction=ifelse(direction=='1 ','forward','reverse'),
         encoded=ifelse(str_detect(AGI,'ATM')==T,'Mitochondrion',
                                   ifelse(str_detect(AGI,'ATC')==T,'Plastid','Nucleus'))) %>% 
  dplyr::select(-1) %>% 
  distinct(AGI,.keep_all = T)


reslfc <- reslfc %>% 
  left_join(pos, by=c('rowname'='AGI'))

reslfc <- reslfc %>% 
  dplyr::rename(AGI=rowname)


#write table

write.csv(reslfc,'data/rnaseq_table_v1.csv')










# count plots for organelles====
#Start with table, make long, facet for organelles and coding

results <- fread('data/rnaseq_table_v1.csv') %>% 
  gather(type,value,clp1_lognfc:clp2_padj) %>% 
  mutate(genotype=sapply(strsplit(type,'_'),'[',1),
         type=sapply(strsplit(type,'_'),'[',2)) %>% 
  spread(type,value) %>% 
  mutate(lognfc=as.numeric(lognfc),
         padj=as.numeric(padj)) %>% 
  na.omit() %>% 
  mutate(location_consensus=gsub('endoplasmic,reticulum','endoplasmic reticulum',location_consensus),
         location_consensus=gsub('plasma,membrane','plasma membrane',location_consensus)) %>% 
  rowwise() %>% 
  mutate(location_consensus=paste(sort(unlist(strsplit(location_consensus,','))),collapse = ',')) %>% 
  dplyr::filter(str_detect(location_consensus,',')==F) %>% 
  mutate(regulation=ifelse(lognfc > 0.5 & padj <= 0.05,'up',
                           ifelse(lognfc < -0.5 & padj <= 0.05,'down',
                                  ifelse(padj > 0.05 | lognfc == 0 ,'none','failsave'))))


total <- results %>% 
  xtabs(formula = ~ genotype + encoded) %>%
  as.data.frame()
regulation <- results %>% 
  xtabs(formula = ~ genotype + encoded + regulation) %>%
  as.data.frame()


#add median as ref y axis point 
y_ref <- subset %>% group_by(genotype,desc) %>% 
  summarise(median=median(count),max=max(count))
res <- res %>% 
  left_join(y_ref)

#levels
subset$genotype <- factor(subset$genotype, levels = c('WT','clp1','clp2'))
subset$desc <- factor(subset$desc, levels = c('CLPP2','Tubulin6','ATVDAC1'))

res$genotype <- factor(res$genotype, levels = c('WT','clp1','clp2'))
res$desc <- factor(res$desc, levels = c('CLPP2','Tubulin6','ATVDAC1'))



#plot
p <- ggplot(results, aes(encoded, lognfc,fill=genotype,col=genotype)) + 
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar",col='black', size = 0.3,position = position_dodge(width=1))+
  geom_jitter(pch = 21,size=2,color='black',alpha=0.05,position = position_jitterdodge(jitter.width = 0.1,
                                                                            dodge.width = 1))+
  geom_text(data=dplyr::filter(regulation,regulation=='up'),aes(encoded,5,label=Freq),size=4,col='black', lineheight = 0.25,position = position_dodge(width=1))+
  geom_text(data=dplyr::filter(regulation,regulation=='down'),aes(encoded,-5,label=Freq),size=4,col='black', lineheight = 0.25,position = position_dodge(width=1))+
  geom_text(data=dplyr::filter(regulation,regulation=='none'),aes(encoded,0.8,label=Freq),size=4, col='black',lineheight = 0.25,position = position_dodge(width=1))+
  expand_limits(y=0)+
  scale_colour_manual(values=c('#3399cc','#3366cc'))+
  scale_fill_manual(values=c('#3399cc','#3366cc'))+
  labs(title='Total transcript abundance',y='lognFC',x='encoding origin')+
  theme(axis.title.x = element_text(face='bold',size='8'),
        axis.text.x = element_text(face=c('plain','italic','italic'),size=8,angle=30),
        axis.title.y = element_text(face='bold',size='8'),
        axis.text.y = element_text(face='bold',size=8),
        strip.text = element_text(face='bold',size=8),
        title=element_text(size=10))

#save
ggsave('transcript per origin.pdf',device = 'pdf',dpi=1080,plot = p,height = 10,width = 16,units = 'cm')
