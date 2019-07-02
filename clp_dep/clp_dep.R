

# libraries ---------------------------------------------------------------
library(DEP)
library(tidyverse)
library(data.table)
library(fitdistrplus)
library(Biostrings)
library(SummarizedExperiment)

###### read data ---------------------------------------------------------------
data <- fread('data/proteinGroups.txt')
data <- filter(data, Reverse != "+", `Potential contaminant` != "+")
colnames(data)

##### Data preparation #######



#extract gene names from fasta header and sort
data2 <- data %>%
  rowwise() %>% 
  mutate(Gene.names=paste(filter(as.tibble(unlist(strsplit(`Fasta headers`,'|',fixed = T))),str_detect(value,'AT.G')==F,str_detect(value,'Chr')==F,is.na(as.numeric(value)))$value, collapse = ';'))
data2 <- data2[,c(1,2,ncol(data2),3:(ncol(data2)-1))]

## remove pool and blank
data2 <- data2[,which(str_detect(colnames(data2),'blan')==F)]
data2 <- data2[,which(str_detect(colnames(data2),'pool')==F)]
data2 <- data2[,which(str_detect(colnames(data2),'clp2p_4')==F)]

##add better annotations of complex I and remainder
cI_anno <- fread('data/complex I annotation.csv', header = T) %>% dplyr::select(1:3) %>% as.tibble()

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

all_anno <- readAAStringSet('data/Arabidopsis_thaliana.TAIR10.pep.all.fa') %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  mutate(ACC=sapply(strsplit(rowname,'[.]'),'[',1),
         desc=sapply(strsplit(rowname,'symbol:',fixed=T),'[',2),
         desc=ifelse(is.na(desc),sapply(strsplit(rowname,' description:',fixed=T),'[',2),sapply(strsplit(desc,' description:',fixed=T),'[',1)),
         desc=sapply(strsplit(desc,';',fixed=T),'[',1),
         desc=sapply(strsplit(desc,' [Source',fixed=T),'[',1)) %>% 
  dplyr::select(ACC,desc)

all_anno <- all_anno %>% left_join(cI_anno) %>% 
  mutate(desc=ifelse(is.na(Name),desc,Name)) %>%
  dplyr::select(-Name) %>% 
  as.tibble() %>% 
  distinct(ACC, .keep_all = T)

#extract gene names from fasta header and sort

data3 <- separate_rows(data2,`Majority protein IDs`, convert = T) %>% 
  mutate(`Majority protein IDs`=toupper(substr(`Majority protein IDs`,1,9))) %>% 
  left_join(all_anno, by=c('Majority protein IDs'='ACC'))

data3 <- data3[,c(1:3,221,4:220)] 
data3 <- data3 %>% as.tibble() %>%  
  dplyr::distinct(`Majority protein IDs`,desc,.keep_all = T)

##take best names(remove unknown)

data3 <- data3 %>% 
  mutate(Gene.names=ifelse(str_detect(desc,'unknown')==T,Gene.names,desc))

#check for duplicated gene names and remove isoform duplications in gene names          
data3$Gene.names %>% duplicated() %>% any()

data2 <- data3 %>% 
  rowwise() %>% 
  mutate(Gene.names=paste(unique(unlist(strsplit(Gene.names, ';'))),collapse=';'))




# Make unique names using the annotation in the "Gene.names" column as primary names and the annotation in "Protein.IDs" as name for those that do not have an gene name.
data_unique <- make_unique(data2, "Gene.names", "Protein IDs", delim = ";")




##### Generate a SummarizedExperiment object and run DEP #######
# Generate a SummarizedExperiment object by parsing condition information from the column names
LFQ_columns <- grep("LFQ.", colnames(data_unique)) # get LFQ column numbers
data_se <- make_se_parse(data_unique, LFQ_columns)
data_se


# #normal distribution ?
# #before log2
# data_dist <- data_unique[,which(str_detect(colnames(data_unique),'LFQ')==T)] %>% 
#   gather(column,value)
# descdist(data_dist$value)
# hist(data_dist$value,col = 'royalblue',breaks = 500,xlim = c(0,50000000))
# #after
# data_dist_after <- assay(data_se) %>% as.data.frame() %>% 
#   gather(column,value) %>% 
#   na.omit()
# descdist(data_dist_after$value)
# hist(data_dist_after$value, col='royalblue')

# Less stringent filtering:
# Filter for proteins that are identified in 2 out of 3 replicates of at least one condition
data_filt2 <- filter_missval(data_se, thr = 1)
data_filt <- data_filt2

# Normalize the data
data_norm <- normalize_vsn(data_filt)
# Impute missing data using random draws from a Gaussian distribution centered around a minimal value (for MNAR)
data_imp <- impute(data_norm, fun = "MinProb", q = 0.01)
data_diff_manual <- test_diff(data_imp, type = "manual", 
                              test = c("clp1s__vs_wts_", "clp2s__vs_wts_",'clp1p__vs_wtp_','clp2p__vs_wtp_'))
data_diff <- data_diff_manual

# Denote significant proteins based on user defined cutoffs ---------------
dep <- add_rejections(data_diff, alpha = 0.1, lfc = log2(1))
dep <- add_rejections(data_diff, alpha = 0.05, lfc = log2(1))
data_results <- get_results(dep)
data_results %>% filter(significant) %>% nrow()

##### PCA PLOT #######

# Plot the first and second principal components
plot_pca(dep, x = 1, n=nrow(dep),y = 2, point_size = 4)


##### Cor Matrix #######


# Plot the Pearson correlation matrix
plot_cor(dep, significant = TRUE, lower = 0, upper = 1, pal = "Reds")


##### Heatmap #######


# Plot a heatmap of all significant proteins with the data centered per protein
plot_heatmap(dep, type = "centered", kmeans = TRUE, 
             k = 6, col_limit = 4, show_row_names = F,
             indicate = c("condition", "replicate"))

# Plot a heatmap of all significant proteins (rows) and the tested contrasts (columns)
plot_heatmap(dep, type = "contrast", kmeans = TRUE, 
             k = 6, col_limit = 10, show_row_names = FALSE)


##### Volcano Plot #######

# Plot a volcano plot for the contrast "Ubi6 vs Ctrl""
plot_volcano(dep, contrast = "clp1s__vs_wts_", label_size = 3, add_names = TRUE)
plot_volcano(dep, contrast = "clp2s__vs_wts_", label_size = 3, add_names = TRUE)

##### Barplots of a protein of interest #######

#change levels
dep_save <- dep
dep@colData$condition <- gsub('_', '',dep@colData$condition)
dep@colData$condition <- factor(dep@colData$condition, 
                                levels = c('wtp','clp1p','clp2p','wts','clp1s','clp2s',4))
# Plot a barplot for USP15 and IKBKG
plot_single(dep, proteins = ' presequence protease 1 ')
# Plot a barplot for the protein USP15 with the data centered
plot_single(dep, proteins = ' lon protease 1 ', type = "centered")



##### Results table #######


# Generate a results table
dep <- add_rejections(data_diff, alpha = 0.1, lfc = log2(1))

data_results <- get_results(dep)

# Number of significant proteins
data_results %>% filter(significant) %>% nrow()

# Column names of the results table
colnames(data_results)


#custom single plot
#old

#remove 'significant' colnames

data_results <- data_results[,which(str_detect(colnames(data_results),'significant')==F)]

#long format(gather)
data_long <- data_results %>% gather(sample,value,3:ncol(.)) %>% 
  mutate(type=sapply(strsplit(sample,'__'),'[',3),
         type=ifelse(is.na(type)==T,sapply(strsplit(sample,'__'),'[',2),type),
         sample=sapply(strsplit(sample,'__'),'[',1),
         fraction=ifelse(sample %in% c("clp1p","clp2p", "wtp"),'membrane','soluble'),
         genotype=substr(sample,1,nchar(sample)-1)) %>% 
  dplyr::select(-sample) %>%
  rowwise() %>% 
  spread(key = type,value = value) %>% 
  mutate(ID=substr(ID,1,9))



#levels for genotype and fraction

data_long$genotype <- factor(data_long$genotype, levels = c('wt','clp1','clp2'))
data_long$fraction <- factor(data_long$fraction, levels = c('membrane','soluble'))


#plot gene of interest
dep <- add_rejections(data_diff, alpha = 0.1, lfc = log2(1))

# rownames(dep) <- ifelse(substr(rownames(dep),nchar(rownames(dep)),
#                                nchar(rownames(dep)))==' ',
#                         substr(rownames(dep),2,nchar(rownames(dep))-1),
#                         substr(rownames(dep),2,nchar(rownames(dep))))

prot_of_interest <- 'AT1G48030'

subset <- dep[filter(data_long,ID==prot_of_interest)$name[1]]
means <- rowMeans(assay(subset), na.rm = TRUE)
df_reps <- data.frame(assay(subset)) %>% rownames_to_column() %>% 
  gather(ID, val, -rowname) %>% left_join(., data.frame(colData(subset)),by = "ID")

df_reps$replicate <- as.factor(df_reps$replicate) 
df_reps <- df_reps %>% 
  mutate(condition=sapply(strsplit(condition,'_'),'[',1),
         fraction=ifelse(condition %in% c("clp1p","clp2p", "wtp"),'membrane','soluble'),
         genotype=substr(condition,1,nchar(condition)-1))
  

df <- df_reps %>% 
  group_by(condition, rowname) %>% 
  summarize(mean = mean(val,na.rm = TRUE), sd = sd(val, na.rm = TRUE), n = n()) %>% 
  mutate(error = qnorm(0.975) * sd/sqrt(n), CI.L = mean -  error, CI.R = mean + error) %>% 
  as.data.frame() %>% 
  mutate(condition=sapply(strsplit(condition,'_'),'[',1),
         fraction=ifelse(condition %in% c("clp1p","clp2p", "wtp"),'membrane','soluble'),
         genotype=substr(condition,1,nchar(condition)-1))


rowname <- df$rowname[1]
#levels
df$genotype <- factor(df$genotype, levels = c('wt','clp1','clp2'))
df$fraction <- factor(df$fraction, levels = c('membrane','soluble'))

ggplot(df, aes(genotype, mean)) + geom_hline(yintercept = 0) + 
  geom_col(colour = "black", fill = "grey") + 
  geom_errorbar(aes(ymin = CI.L, ymax = CI.R), width = 0.3,size=1.2) +
  geom_jitter(data = df_reps, aes(genotype, val, col = replicate),alpha=0.5,
             size = 5, position = position_dodge(width = 0.3)) + 
  labs(title=paste(prot_of_interest,rowname, sep=' - '),x = "sample", y = expression(log[2] ~ "Centered intensity" ~ "(?95% CI)"), col = "Rep") + 
  facet_wrap(~fraction) + 
  theme_DEP2()



###no log


#remove 'significant' colnames

data_results <- data_results[,which(str_detect(colnames(data_results),'significant')==F)]

#long format(gather)
data_long <- data_results %>% gather(sample,value,3:ncol(.)) %>% 
  mutate(type=sapply(strsplit(sample,'__'),'[',3),
         type=ifelse(is.na(type)==T,sapply(strsplit(sample,'__'),'[',2),type),
         sample=sapply(strsplit(sample,'__'),'[',1),
         fraction=ifelse(sample %in% c("clp1p","clp2p", "wtp"),'membrane','soluble'),
         genotype=substr(sample,1,nchar(sample)-1)) %>% 
  dplyr::select(-sample) %>%
  rowwise() %>% 
  spread(key = type,value = value) %>% 
  mutate(ID=sapply(strsplit(ID,'[.]'),'[',1))


####


#plot gene of interest
dep <- add_rejections(data_diff, alpha = 0.1, lfc = log2(1))

# rownames(dep) <- ifelse(substr(rownames(dep),nchar(rownames(dep)),
#                                nchar(rownames(dep)))==' ',
#                         substr(rownames(dep),2,nchar(rownames(dep))-1),
#                         substr(rownames(dep),2,nchar(rownames(dep))))

prot_of_interest <- filter(data3, `Majority protein IDs` %in% cI_anno$ACC) %>% 
  dplyr::select(2,3) %>% 
  dplyr::rename(desc=Gene.names,ID=`Majority protein IDs`)
 

subset <- dep[unique(filter(data_long,ID %in% prot_of_interest$ID)$name)]
means <- rowMeans(assay(subset), na.rm = TRUE)
df_reps <- data.frame(assay(subset)) %>% rownames_to_column() %>% 
  gather(ID, val, -rowname) %>% left_join(., data.frame(colData(subset)),by = "ID") 

desc <- dplyr::select(data_long,name,ID) %>% 
  left_join(prot_of_interest) %>% 
  na.omit() %>% 
  dplyr::rename(ACC=ID)

df_reps <- left_join(df_reps,desc,by=c('rowname'='name')) %>% distinct() 


df_reps$replicate <- as.factor(df_reps$replicate) 
df_reps <- df_reps %>% 
  mutate(condition=sapply(strsplit(condition,'_'),'[',1),
         fraction=ifelse(condition %in% c("clp1p","clp2p", "wtp"),'membrane','soluble'),
         genotype=substr(condition,1,nchar(condition)-1)) %>% 
  mutate(val=2^val,
         facet=paste(ACC,desc,sep = ' - '))


df <- df_reps %>% 
  group_by(condition, rowname) %>% 
  summarize(mean = mean(val,na.rm = TRUE), sd = sd(val, na.rm = TRUE), n = n()) %>% 
  mutate(error = qnorm(0.975) * sd/sqrt(n), CI.L = mean -  error, CI.R = mean + error) %>% 
  as.data.frame() %>% 
  mutate(condition=sapply(strsplit(condition,'_'),'[',1),
         fraction=ifelse(condition %in% c("clp1p","clp2p", "wtp"),'membrane','soluble'),
         genotype=substr(condition,1,nchar(condition)-1)) %>% 
  left_join(desc,by=c('rowname'='name')) %>% 
  mutate(facet=paste(ACC,desc,sep = ' - ')) %>%
  distinct()

  


rowname <- df$rowname[1]
#levels
df$genotype <- factor(df$genotype, levels = c('wt','clp1','clp2'))
df$fraction <- factor(df$fraction, levels = c('membrane','soluble'))
df_reps$genotype <- factor(df_reps$genotype, levels = c('wt','clp1','clp2'))
df_reps$fraction <- factor(df_reps$fraction, levels = c('membrane','soluble'))

facet=unique(df$facet)

#order of facets
ACC_order <- data.frame(ACC=toupper(c('At5g08530','At4g02580','At5g37510','At1g16700','At1g79010','At5g11770','AtMg00070','AtMg00510','AtMg00516','AtMg01120','AtMg01275','AtMg00285','AtMg01320','AtMg00990','AtMg00650',
  'AtMg00580','AtMg00513','AtMg0060','AtMg00665','AtMg00270','At5g67590','At3g03070','At3g08610','At5g47890','At3g03100','At1g67785','At5g52840','At5g08060','At2g20360',
  'At1g72040','At2g46540','At3g12260','At5g18800','At3g06310','At2g42210','At1g04630','At2g33220','At4g16450','At1g76200','At2g02510','At1g14450','At2g31490','At2g02050','At5g47570',
  'At4g34700','At1g49140','At3g18410','At3g57785','At2g42310','At4g00585','At4g20150','At3g62790','At2g47690','At1g19580','At1g47260','At5g66510','At5g63510','At3g48680',
  'At1g67350','At2g27730','At1g18320','At3g07480','At5g14105','At1g68680','At3g10110','At1g72170','At2g28430','At1g72750')))

ACC_order <- ACC_order %>% left_join(dplyr::select(df,ACC,desc)) %>% distinct(ACC,.keep_all = T) %>% 
  mutate(facet=paste(ACC,desc,sep = ' - ')) %>% 
  distinct(facet)

df$facet <- factor(df$facet, levels=(ACC_order$facet))
df_reps$facet <- factor(df_reps$facet, levels=(ACC_order$facet))
#all single plots

#ggplot(df, aes(fraction, mean/1000000, fill=genotype)) + geom_hline(yintercept = 0) + 
 # geom_bar(stat='identity',position=position_dodge(width=1)) + 
  #geom_errorbar(aes(ymin = CI.L/1000000, ymax = CI.R/1000000),position=position_dodge(width=1), width = 0.3,size=1.2) +
  #geom_jitter(data = df_reps, aes(fraction, val/1000000,group=genotype, col = replicate),alpha=0.5,
  #            size = 1, position = position_dodge(width = 1)) + 
  #labs(title='LFQ Raw intensities CLP',x = "", y = ("LFQ Intensity" ~ "(?95% CI) [Million Area Count]"), col = "Rep") + 
  #facet_wrap(~facet,scales='free') + 
  #scale_fill_manual(values=c('#339900','#0066cc','#990033'))+
  #scale_y_continuous(label=unit_format(suffix=' M'))+
  #theme(axis.text.x = element_text(angle=0,hjust=0.5))

#stats
library(broom)
all_ttest <- data.frame()
for (prot in unique(df$ACC)){
  for(frac in unique(df$fraction)){
    stat=filter(df_reps,ACC==prot,fraction==frac) 
    stat=tidy(pairwise.t.test(stat$val,stat$genotype,p.adjust.method = 'BH')) %>%
      as.data.frame() %>% 
      mutate(fraction=frac,ACC=prot,
             p.value=round(p.value,digits=3)) %>% 
      filter(group2=='wt') %>% 
      dplyr::select(-group2) %>% 
      dplyr::rename(genotype=group1)
    all_ttest=bind_rows(all_ttest,stat)
  }
}
df <- left_join(df,all_ttest) %>% 
  mutate(p.value=ifelse(is.na(p.value),'1',p.value),sig=ifelse(p.value <= 0.05,'*',''))

df$genotype <- factor(df$genotype, levels = c('wt','clp1','clp2'))
df$fraction <- factor(df$fraction, levels = c('membrane','soluble'))
df_reps$genotype <- factor(df_reps$genotype, levels = c('wt','clp1','clp2'))
df_reps$fraction <- factor(df_reps$fraction, levels = c('membrane','soluble'))


#single all plots
for (chart in unique(df$facet)){
  filename=gsub('/','_',chart)
  p <-ggplot(filter(df, facet==chart), aes(fraction, mean/1000000, fill=genotype)) + geom_hline(yintercept = 0) + 
      geom_bar(stat='identity',position=position_dodge(width=1)) + 
      geom_errorbar(aes(ymin = CI.L/1000000, ymax = CI.R/1000000),position=position_dodge(width=1), width = 0.3,size=1.2) +
      geom_jitter(data = filter(df_reps, facet==chart), aes(fraction, val/1000000,group=genotype, col = replicate),alpha=0.5,
                  size = 3, position = position_dodge(width = 1)) + 
      geom_text(aes(y=CI.R/1000000+CI.R/1000000*0.1,label=sig),stat='identity',position=position_dodge(width=1))+
      labs(title=chart,x = "", y = ("LFQ Intensity" ~ "(?95% CI) [Million Area Count]"), col = "Rep") + 
      scale_fill_manual(values=c('#339900','#0066cc','#990033'))+
      scale_y_continuous(label=number_format(suffix=' M'))+
      theme_DEP2()+
      theme(axis.text.x = element_text(angle=0,hjust=0.5))
  ggsave(filename = paste0(filename,'.png'),path = 'images',device = 'png',dpi=1080,plot = p)
  message(paste0(filename,' done!'))
}




## custom mutant overlap volcano

#data

dep <- add_rejections(data_diff, alpha = 0.1, lfc = log2(1))

data_results <- get_results(dep)

#long format(gather)
data_long <- data_results %>% gather(sample,value,3:ncol(.)) %>% 
  mutate(type=sapply(strsplit(sample,'__'),'[',3),
         type=ifelse(is.na(type)==T,sapply(strsplit(sample,'__'),'[',2),type),
         sample=sapply(strsplit(sample,'__'),'[',1),
         fraction=ifelse(sample %in% c("clp1p","clp2p", "wtp"),'membrane','soluble'),
         genotype=substr(sample,1,nchar(sample)-1)) %>% 
  dplyr::select(-sample) %>%
  rowwise() %>% 
  filter(type != 'significant') %>% 
  spread(key = type,value = value) %>% 
  mutate(ID=sapply(strsplit(ID,'[.]'),'[',1))


#filter for both mutants under 0.1pvalue, and use lowest FC and highest pvalue
data_volcano <-  
  data_long %>% 
  filter(genotype !='wt') %>% 
  group_by(ID,fraction) %>% 
  mutate(colour_p=ifelse(max(p.adj) <= 0.1,'red','blue')) %>% 
  mutate(min_ratio=min(abs(ratio)),
         max_p=max(p.adj)) %>% 
  filter(abs(ratio)==min_ratio) %>% 
  dplyr::select(-min_ratio) %>% 
  mutate(colour_r=ifelse(ratio <=-0.4 | ratio >= 0.4,'red','blue')) %>% 
  mutate(sig=ifelse(colour_p=='blue'|colour_r=='blue','non_sig','sig')) %>% 
  distinct(ID,fraction, .keep_all = T)

#add good annotation
desc_volcano <- data2 %>% 
  dplyr::select(Gene.names,`Majority protein IDs`) %>% 
  dplyr::rename(desc=Gene.names,ACC=`Majority protein IDs`)

#hand curate some annotation

desc_volcano <- desc_volcano %>% 
  mutate(desc=ifelse(ACC=='AT5G08690','ATP synth beta2',
                     ifelse(ACC=='AT5G62270','mucin related AT5G62270',
                            ifelse(ACC=='AT1G26460','PPR AT1G26460',
                                   ifelse(ACC=='AT3G02650','PPR AT3G02650',
                                          ifelse(ACC=='AT5G64670','PPR AT5G64670',desc))))))

desc_volcano <- desc_volcano %>% 
  mutate(desc=ifelse(ACC=='AT3G62530','ARM repeat',
                     ifelse(ACC=='AT4G21020','Late embryo abundant',
                            ifelse(ACC=='AT5G55200','CoChaperone GrpE',
                                   ifelse(ACC=='AT3G19740','p-loop hydrolase',
                                          ifelse(ACC=='AT4G36680','PPR AT4G36680',
                                                 ifelse(ACC=='AT1G70190','Ribo L7/L12',desc)))))))






#merge annnotation into volcano data

data_volcano <- data_volcano %>% 
  left_join(desc_volcano, by=c('ID'='ACC')) 


#levels
data_volcano$sig <- factor(data_volcano$sig, levels=c('sig','non_sig'))

## format data for plot
  #change-log10 pvalue  high numbers to a maximum of 7.5
  #change symbol of  >7.5 to op arrow

data_volcano <- data_volcano %>% 
  mutate(max_p_adj=ifelse(-log10(max_p) > 7.5,0.0000000316215,max_p),
         pch=ifelse(-log10(max_p) > 7.5,17,16),
         ratio_adj=ifelse(ratio > 3.5,3.5,ifelse(ratio < -3.5,-3.5,ratio))) 
  


library(ggrepel)

p <- ggplot(data_volcano, aes(x=ratio_adj,y=-log10(max_p_adj),col=sig))+
  geom_point(pch=data_volcano$pch,alpha=0.75,size=2)+
  geom_text_repel(data=filter(data_volcano,sig=='sig'),aes(label=desc),col='black',size=2.5, fontface='bold')+
  facet_wrap(~fraction)+
  scale_colour_manual(values=c('#990000','#99ccff'))+
  theme(legend.position = 'none', axis.title = element_text(face='bold',size = 18),
        axis.text = element_text(face='bold',size = 16), strip.text = element_text(face='bold',size=18),
        title = element_text(face='bold',size=18))+
  labs(title='Clp vs col-O', x='log2 fold change', y='-log10 p-value')

ggsave(filename = paste0('volcano','.png'),path = 'images',device = 'png',dpi=1080,plot = p)



###custom table significant genes


dep <- add_rejections(data_diff, alpha = 0.1, lfc = log2(1))

data_results <- get_results(dep)
data_results <- data_results[,which(str_detect(colnames(data_results),'significant')==F)]
data_results <- data_results[,which(str_detect(colnames(data_results),'p.val')==F)]
data_results <- data_results[,which(str_detect(colnames(data_results),'center')==F)]

data_results <- data_results[,c(2,1,3,5,7,9,4,6,8,10)]

write.csv(data_results,'data/results_clp.csv')

#mapman annotation

data_results <- data_results %>% 
  mutate(ID=substr(ID,1,9))

mapman <- fread('data/X4_Araport11_R1.0.txt') %>% 
  mutate_all(.,funs(gsub("'","",.))) %>% 
  dplyr::rename(ID=IDENTIFIER) %>% 
  mutate(ID=toupper(substr(ID,1,9))) %>% 
  dplyr::select(-NAME)

data_results <- data_results %>% left_join(mapman)

##complex I subunits


cI_agi <- data_long %>% 
  filter(ID %in% cI_anno$ACC)
cI_agi <- cI_agi$ID


dep <- add_rejections(data_diff, alpha = 1, lfc = log2(1))



subset <- dep[unique(filter(data_long,ID %in% cI_agi)$name)]
means <- rowMeans(assay(subset), na.rm = TRUE)
df_reps <- data.frame(assay(subset)) %>% rownames_to_column() %>% 
  gather(ID, val, -rowname) %>% left_join(., data.frame(colData(subset)),by = "ID") 

desc <- dplyr::select(data_long,name,ID) %>% 
  filter(ID %in% cI_agi) %>% 
  na.omit() %>% 
  dplyr::rename(ACC=ID) %>% 
  distinct(ACC, .keep_all = T)

df_reps <- left_join(df_reps,desc,by=c('rowname'='name')) %>% distinct() %>% na.omit 


df_reps$replicate <- as.factor(df_reps$replicate) 
df_reps <- df_reps %>% 
  filter(ID !='clp2s__4') %>% 
  ungroup() %>% 
  mutate(condition=sapply(strsplit(condition,'_'),'[',1),
         genotype=substr(condition,1,nchar(condition)-1)) %>% 
  group_by(genotype,ACC,replicate) %>% 
  mutate(val=2^mean(val),
         facet=paste(ACC,rowname,sep = ' - ')) %>% 
  distinct(genotype,ACC,replicate,.keep_all = T)


df <- df_reps %>% 
  group_by(condition, rowname) %>% 
  summarize(mean = mean(val,na.rm = TRUE), sd = sd(val, na.rm = TRUE), n = n()) %>% 
  mutate(error = qnorm(0.975) * sd/sqrt(n), CI.L = mean -  error, CI.R = mean + error) %>% 
  as.data.frame() %>% 
  mutate(condition=sapply(strsplit(condition,'_'),'[',1),
         genotype=substr(condition,1,nchar(condition)-1)) %>% 
  left_join(desc,by=c('rowname'='name')) %>% 
  mutate(facet=paste(ACC,rowname,sep = ' - ')) %>%
  distinct()




rowname <- df$rowname[1]
#levels
df$genotype <- factor(df$genotype, levels = c('wt','clp1','clp2'))
df_reps$genotype <- factor(df_reps$genotype, levels = c('wt','clp1','clp2'))

facet=unique(df$facet)

df$facet <- factor(df$facet, levels=unique(df$facet[order(df$rowname)]))
df_reps$facet <- factor(df_reps$facet, levels=unique(df_reps$facet[order(df_reps$rowname)]))

#all single plots

p <- ggplot(df, aes(ACC, mean/1000000, fill=genotype)) + geom_hline(yintercept = 0) +
  geom_bar(stat='identity',position=position_dodge(width=1)) +
  geom_errorbar(aes(ymin = CI.L/1000000, ymax = CI.R/1000000),position=position_dodge(width=1), width = 0.3,size=1.2) +
  geom_jitter(data = df_reps, aes(ACC, val/1000000,group=genotype, col = replicate),alpha=0.5,
              size = 1, position = position_dodge(width = 1)) +
  labs(title='Complex I subunits',x = "", y = ("LFQ Intensity" ~ "(?95% CI) [Million Area Count]"), col = "Rep") +
  facet_wrap(~rowname,scales='free') +
  scale_fill_manual(values=c('#339900','#0066cc','#990033'))+
  #scale_y_continuous(label=unit_format(suffix=' M'))+
  theme(axis.text.x = element_blank(),axis.title.y = element_text(size=16,face='bold'),axis.text.y = element_text(size=12),
        strip.text = element_text(size=12,face='bold'),legend.text=element_text(size=14,face='bold'),legend.key.size = unit(2,'cm'))

ggsave(filename = paste0('complex I survey','.png'),path = 'images',width=16,height = 10,device = 'png',dpi=1080,plot = p)





##mitochondrial encoded

data_long %>% filter(str_detect(ID,'ATM')==T) %>% View()

mito_agi <- data_long%>% 
  filter(str_detect(ID,'ATM')==T)
mito_agi <- mito_agi$ID


dep <- add_rejections(data_diff, alpha = 1, lfc = log2(1))



subset <- dep[unique(filter(data_long,ID %in% mito_agi)$name)]
means <- rowMeans(assay(subset), na.rm = TRUE)
df_reps <- data.frame(assay(subset)) %>% rownames_to_column() %>% 
  gather(ID, val, -rowname) %>% left_join(., data.frame(colData(subset)),by = "ID") 

desc <- dplyr::select(data_long,name,ID) %>% 
  filter(ID %in% mito_agi) %>% 
  na.omit() %>% 
  dplyr::rename(ACC=ID) %>% 
  distinct(ACC, .keep_all = T)

df_reps <- left_join(df_reps,desc,by=c('rowname'='name')) %>% distinct() 


df_reps$replicate <- as.factor(df_reps$replicate) 
df_reps <- df_reps %>% 
  mutate(condition=sapply(strsplit(condition,'_'),'[',1),
         fraction=ifelse(condition %in% c("clp1p","clp2p", "wtp"),'membrane','soluble'),
         genotype=substr(condition,1,nchar(condition)-1)) %>% 
  mutate(val=2^val,
         facet=paste(ACC,rowname,sep = ' - '))


df <- df_reps %>% 
  group_by(condition, rowname) %>% 
  summarize(mean = mean(val,na.rm = TRUE), sd = sd(val, na.rm = TRUE), n = n()) %>% 
  mutate(error = qnorm(0.975) * sd/sqrt(n), CI.L = mean -  error, CI.R = mean + error) %>% 
  as.data.frame() %>% 
  mutate(condition=sapply(strsplit(condition,'_'),'[',1),
         fraction=ifelse(condition %in% c("clp1p","clp2p", "wtp"),'membrane','soluble'),
         genotype=substr(condition,1,nchar(condition)-1)) %>% 
  left_join(desc,by=c('rowname'='name')) %>% 
  mutate(facet=paste(ACC,rowname,sep = ' - ')) %>%
  distinct()




rowname <- df$rowname[1]
#levels
df$genotype <- factor(df$genotype, levels = c('wt','clp1','clp2'))
df$fraction <- factor(df$fraction, levels = c('membrane','soluble'))
df_reps$genotype <- factor(df_reps$genotype, levels = c('wt','clp1','clp2'))
df_reps$fraction <- factor(df_reps$fraction, levels = c('membrane','soluble'))

facet=unique(df$facet)

df$facet <- factor(df$facet, levels=unique(df$facet[order(df$rowname)]))
df_reps$facet <- factor(df_reps$facet, levels=unique(df_reps$facet[order(df_reps$rowname)]))

#all single plots

p <- ggplot(filter(df,fraction=='soluble'), aes(fraction, mean/1000000, fill=genotype)) + geom_hline(yintercept = 0) +
  geom_bar(stat='identity',position=position_dodge(width=1)) +
  geom_errorbar(aes(ymin = CI.L/1000000, ymax = CI.R/1000000),position=position_dodge(width=1), width = 0.3,size=1.2) +
  geom_jitter(data = filter(df_reps,fraction=='soluble'), aes(fraction, val/1000000,group=genotype, col = replicate),alpha=0.5,
              size = 1, position = position_dodge(width = 1)) +
  labs(title='Mitochondrial encoded Proteins',x = "", y = ("LFQ Intensity" ~ "(?95% CI) [Million Area Count]"), col = "Rep") +
  facet_wrap(~rowname,scales='free') +
  scale_fill_manual(values=c('#339900','#0066cc','#990033'))+
  #scale_y_continuous(label=unit_format(suffix=' M'))+
  theme(axis.text.x = element_blank(),axis.title.y = element_text(size=16,face='bold'),axis.text.y = element_text(size=12),
        strip.text = element_text(size=12,face='bold'),legend.text=element_text(size=14,face='bold'),legend.key.size = unit(2,'cm'))

ggsave(filename = paste0('mito encoded survey','.png'),path = 'images',width=16,height = 10,device = 'png',dpi=1080,plot = p)




##protease AGI
protease_agi <- fread('data/protease_agi.csv') %>% 
  mutate(Gene=toupper(substr(Gene,1,9)))

dep <- add_rejections(data_diff, alpha = 1, lfc = log2(1))



subset <- dep[unique(filter(data_long,ID %in% protease_agi$Gene)$name)]
means <- rowMeans(assay(subset), na.rm = TRUE)
df_reps <- data.frame(assay(subset)) %>% rownames_to_column() %>% 
  gather(ID, val, -rowname) %>% left_join(., data.frame(colData(subset)),by = "ID") 

desc <- dplyr::select(data_long,name,ID) %>% 
  left_join(protease_agi, by=c('ID'='Gene')) %>% 
  na.omit() %>% 
  dplyr::rename(ACC=ID) %>% 
  distinct(ACC, .keep_all = T)

df_reps <- left_join(df_reps,desc,by=c('rowname'='name')) %>% distinct() 


df_reps$replicate <- as.factor(df_reps$replicate) 
df_reps <- df_reps %>% 
  mutate(condition=sapply(strsplit(condition,'_'),'[',1),
         fraction=ifelse(condition %in% c("clp1p","clp2p", "wtp"),'membrane','soluble'),
         genotype=substr(condition,1,nchar(condition)-1)) %>% 
  mutate(val=2^val,
         facet=paste(ACC,Name,sep = ' - '))


df <- df_reps %>% 
  group_by(condition, rowname) %>% 
  summarize(mean = mean(val,na.rm = TRUE), sd = sd(val, na.rm = TRUE), n = n()) %>% 
  mutate(error = qnorm(0.975) * sd/sqrt(n), CI.L = mean -  error, CI.R = mean + error) %>% 
  as.data.frame() %>% 
  mutate(condition=sapply(strsplit(condition,'_'),'[',1),
         fraction=ifelse(condition %in% c("clp1p","clp2p", "wtp"),'membrane','soluble'),
         genotype=substr(condition,1,nchar(condition)-1)) %>% 
  left_join(desc,by=c('rowname'='name')) %>% 
  mutate(facet=paste(ACC,Name,sep = ' - ')) %>%
  distinct()




rowname <- df$rowname[1]
#levels
df$genotype <- factor(df$genotype, levels = c('wt','clp1','clp2'))
df$fraction <- factor(df$fraction, levels = c('membrane','soluble'))
df_reps$genotype <- factor(df_reps$genotype, levels = c('wt','clp1','clp2'))
df_reps$fraction <- factor(df_reps$fraction, levels = c('membrane','soluble'))

facet=unique(df$facet)

df$facet <- factor(df$facet, levels=unique(df$facet[order(df$Name)]))
df_reps$facet <- factor(df_reps$facet, levels=unique(df_reps$facet[order(df_reps$Name)]))

#all single plots

p <- ggplot(filter(df,fraction=='soluble'), aes(fraction, mean/1000000, fill=genotype)) + geom_hline(yintercept = 0) +
  geom_bar(stat='identity',position=position_dodge(width=1)) +
  geom_errorbar(aes(ymin = CI.L/1000000, ymax = CI.R/1000000),position=position_dodge(width=1), width = 0.3,size=1.2) +
  geom_jitter(data = filter(df_reps,fraction=='soluble'), aes(fraction, val/1000000,group=genotype, col = replicate),alpha=0.5,
             size = 1, position = position_dodge(width = 1)) +
  labs(title='LFQ Raw intensities CLP',x = "", y = ("LFQ Intensity" ~ "(?95% CI) [Million Area Count]"), col = "Rep") +
  facet_wrap(~facet,scales='free') +
  scale_fill_manual(values=c('#339900','#0066cc','#990033'))+
  #scale_y_continuous(label=scales::number_format(suffix=' M'))+
  theme(axis.text.x = element_blank(),axis.title.y = element_text(size=16,face='bold'),axis.text.y = element_text(size=12),
        strip.text = element_text(size=12,face='bold'),legend.text=element_text(size=14,face='bold'),legend.key.size = unit(2,'cm'))

ggsave(filename = paste0('protease survey','.png'),path = 'images',width=16,height = 10,device = 'png',dpi=1080,plot = p)





#stats
library(broom)
all_ttest <- data.frame()
for (prot in unique(df$ACC)){
  for(frac in unique(df$fraction)){
    stat=filter(df_reps,ACC==prot,fraction==frac) 
    stat=tidy(pairwise.t.test(stat$val,stat$genotype,p.adjust.method = 'BH')) %>%
      as.data.frame() %>% 
      mutate(fraction=frac,ACC=prot,
             p.value=round(p.value,digits=3)) %>% 
      filter(group2=='wt') %>% 
      dplyr::select(-group2) %>% 
      dplyr::rename(genotype=group1)
    all_ttest=bind_rows(all_ttest,stat)
  }
}
df <- left_join(df,all_ttest) %>% 
  mutate(p.value=ifelse(is.na(p.value),'1',p.value),sig=ifelse(p.value <= 0.05,'*',''))

df$genotype <- factor(df$genotype, levels = c('wt','clp1','clp2'))
df$fraction <- factor(df$fraction, levels = c('membrane','soluble'))
df_reps$genotype <- factor(df_reps$genotype, levels = c('wt','clp1','clp2'))
df_reps$fraction <- factor(df_reps$fraction, levels = c('membrane','soluble'))


#single all plots
for (chart in unique(df$facet)){
  filename=gsub('/','_',chart)
  p <-ggplot(filter(df, facet==chart), aes(fraction, mean/1000000, fill=genotype)) + geom_hline(yintercept = 0) + 
    geom_bar(stat='identity',position=position_dodge(width=1)) + 
    geom_errorbar(aes(ymin = CI.L/1000000, ymax = CI.R/1000000),position=position_dodge(width=1), width = 0.3,size=1.2) +
    geom_jitter(data = filter(df_reps, facet==chart), aes(fraction, val/1000000,group=genotype, col = replicate),alpha=0.5,
                size = 3, position = position_dodge(width = 1)) + 
    geom_text(aes(y=CI.R/1000000+CI.R/1000000*0.1,label=sig),stat='identity',position=position_dodge(width=1))+
    labs(title=chart,x = "", y = ("LFQ Intensity" ~ "(?95% CI) [Million Area Count]"), col = "Rep") + 
    scale_fill_manual(values=c('#339900','#0066cc','#990033'))+
    scale_y_continuous(label=number_format(suffix=' M'))+
    theme_DEP2()+
    theme(axis.text.x = element_text(angle=0,hjust=0.5))
  ggsave(filename = paste0(filename,'.png'),path = 'images',device = 'png',dpi=1080,plot = p)
  message(paste0(filename,' done!'))
}




##predicted targets
predicted_agi <- fread('data/mito_clp_targets.tsv') %>% 
  mutate(Gene=toupper(substr(AGI,1,9)))

dep <- add_rejections(data_diff, alpha = 1, lfc = log2(1))



subset <- dep[unique(filter(data_long,ID %in% predicted_agi$Gene)$name)]
means <- rowMeans(assay(subset), na.rm = TRUE)
df_reps <- data.frame(assay(subset)) %>% rownames_to_column() %>% 
  gather(ID, val, -rowname) %>% left_join(., data.frame(colData(subset)),by = "ID") 

desc <- dplyr::select(data_long,name,ID) %>% 
  left_join(predicted_agi, by=c('ID'='Gene')) %>% 
  na.omit() %>% 
  dplyr::rename(ACC=ID) %>% 
  distinct(ACC, .keep_all = T)

df_reps <- left_join(df_reps,desc,by=c('rowname'='name')) %>% distinct() 


df_reps$replicate <- as.factor(df_reps$replicate) 
df_reps <- df_reps %>% 
  mutate(condition=sapply(strsplit(condition,'_'),'[',1),
         fraction=ifelse(condition %in% c("clp1p","clp2p", "wtp"),'membrane','soluble'),
         genotype=substr(condition,1,nchar(condition)-1)) %>% 
  mutate(val=2^val,
         facet=paste(ACC,rowname,sep = ' - '))


df <- df_reps %>% 
  group_by(condition, rowname) %>% 
  summarize(mean = mean(val,na.rm = TRUE), sd = sd(val, na.rm = TRUE), n = n()) %>% 
  mutate(error = qnorm(0.975) * sd/sqrt(n), CI.L = mean -  error, CI.R = mean + error) %>% 
  as.data.frame() %>% 
  mutate(condition=sapply(strsplit(condition,'_'),'[',1),
         fraction=ifelse(condition %in% c("clp1p","clp2p", "wtp"),'membrane','soluble'),
         genotype=substr(condition,1,nchar(condition)-1)) %>% 
  left_join(desc,by=c('rowname'='name')) %>% 
  mutate(facet=paste(ACC,rowname,sep = ' - ')) %>%
  distinct()




rowname <- df$rowname[1]
#levels
df$genotype <- factor(df$genotype, levels = c('wt','clp1','clp2'))
df$fraction <- factor(df$fraction, levels = c('membrane','soluble'))
df_reps$genotype <- factor(df_reps$genotype, levels = c('wt','clp1','clp2'))
df_reps$fraction <- factor(df_reps$fraction, levels = c('membrane','soluble'))

facet=unique(df$facet)

df$facet <- factor(df$facet, levels=unique(df$facet[order(df$Name)]))
df_reps$facet <- factor(df_reps$facet, levels=unique(df_reps$facet[order(df_reps$Name)]))

#all single plots

p <- ggplot(filter(df,fraction=='soluble'), aes(fraction, mean/1000000, fill=genotype)) + geom_hline(yintercept = 0) +
  geom_bar(stat='identity',position=position_dodge(width=1)) +
  geom_errorbar(aes(ymin = CI.L/1000000, ymax = CI.R/1000000),position=position_dodge(width=1), width = 0.3,size=1.2) +
  geom_jitter(data = filter(df_reps,fraction=='soluble'), aes(fraction, val/1000000,group=genotype, col = replicate),alpha=0.5,
              size = 1, position = position_dodge(width = 1)) +
  labs(title='LFQ Raw intensities CLP',x = "", y = ("LFQ Intensity" ~ "(?95% CI) [Million Area Count]"), col = "Rep") +
  facet_wrap(~rowname,scales='free') +
  scale_fill_manual(values=c('#339900','#0066cc','#990033'))+
  #scale_y_continuous(label=scales::number_format(suffix=' M'))+
  theme(axis.text.x = element_blank(),axis.title.y = element_text(size=16,face='bold'),axis.text.y = element_text(size=12),
        strip.text = element_text(size=10),legend.text=element_text(size=14,face='bold'),legend.key.size = unit(2,'cm'))

ggsave(filename = paste0('predicted survey','.png'),path = 'images/predicted',width=16,height = 10,device = 'png',dpi=1080,plot = p)



df <- na.omit(df)

#stats
library(broom)
all_ttest <- data.frame()
for (prot in unique(df$ACC)){
  for(frac in unique(df$fraction)){
    stat=filter(df_reps,ACC==prot,fraction==frac) 
    stat=tidy(pairwise.t.test(stat$val,stat$genotype,p.adjust.method = 'BH')) %>%
      as.data.frame() %>% 
      mutate(fraction=frac,ACC=prot,
             p.value=round(p.value,digits=3)) %>% 
      filter(group2=='wt') %>% 
      dplyr::select(-group2) %>% 
      dplyr::rename(genotype=group1)
    all_ttest=bind_rows(all_ttest,stat)
  }
}
df <- left_join(df,all_ttest) %>% 
  mutate(p.value=ifelse(is.na(p.value),'1',p.value),sig=ifelse(p.value <= 0.05,'*',''))

df$genotype <- factor(df$genotype, levels = c('wt','clp1','clp2'))
df$fraction <- factor(df$fraction, levels = c('membrane','soluble'))
df_reps$genotype <- factor(df_reps$genotype, levels = c('wt','clp1','clp2'))
df_reps$fraction <- factor(df_reps$fraction, levels = c('membrane','soluble'))


#single all plots
for (chart in unique(df$facet)){
  filename=gsub('/','_',chart)
  filename=gsub('%','_',filename)
  p <-ggplot(filter(df, facet==chart), aes(fraction, mean/1000000, fill=genotype)) + geom_hline(yintercept = 0) + 
    geom_bar(stat='identity',position=position_dodge(width=1)) + 
    geom_errorbar(aes(ymin = CI.L/1000000, ymax = CI.R/1000000),position=position_dodge(width=1), width = 0.3,size=1.2) +
    geom_jitter(data = filter(df_reps, facet==chart), aes(fraction, val/1000000,group=genotype, col = replicate),alpha=0.5,
                size = 3, position = position_dodge(width = 1)) + 
    geom_text(aes(y=CI.R/1000000+CI.R/1000000*0.1,label=sig),stat='identity',position=position_dodge(width=1))+
    labs(title=chart,x = "", y = ("LFQ Intensity" ~ "(?95% CI) [Million Area Count]"), col = "Rep") + 
    scale_fill_manual(values=c('#339900','#0066cc','#990033'))+
    scale_y_continuous(label=number_format(suffix=' M'))+
    theme_DEP2()+
    theme(axis.text.x = element_text(angle=0,hjust=0.5))
  ggsave(filename = paste0(filename,'.png'),path = 'images/predicted',device = 'png',dpi=1080,plot = p)
  message(paste0(filename,' done!'))
}


## plot CI kd



kd <- fread('data/TPC2016-00768-LSBR1_Supplemental_Data_Set_2b.csv', skip = 1) 
kd <- kd[,c(2,3,4,144:147)]

kd <- kd %>% 
  dplyr::rename(AGI=`1st AGI`) %>% 
  mutate(AGI=toupper(substr(AGI,1,9)))

#filter for dep proteins
kd <- kd %>% filter(AGI %in% df$AGI)



#intersect in Complex I
kd <- kd %>% left_join(cI_anno, by=c('AGI'='ACC')) %>% 
  na.omit()

#levels
levels <- arrange(kd,`Average KD (d-1)`) %>% distinct(Name)

kd$Name <- factor(kd$Name, levels=levels$Name)

ggplot(kd, aes(Name,`Average KD (d-1)`))+
  geom_bar(stat='identity',fill='darkgreen')+
  coord_flip()+
  theme(axis.title.y = element_blank())+
  labs(title = 'Complex I kd - Lei')


#intersect in predicted targets
kd <- fread('data/TPC2016-00768-LSBR1_Supplemental_Data_Set_2b.csv', skip = 1) 
kd <- kd[,c(2,3,4,144:147)]

kd <- kd %>% 
  dplyr::rename(AGI=`1st AGI`) %>% 
  mutate(AGI=toupper(substr(AGI,1,9)))

#intersect with predicted AGi
kd <- kd %>% left_join(predicted_agi, by=c('AGI'='Gene')) %>% 
  na.omit()


#filter for dep proteins
kd <- kd %>% filter(AGI %in% df$AGI)


#better name
kd <- kd %>% mutate(desc=sapply(strsplit(description,': '),'[',2)) %>% 
  mutate(desc=sapply(strsplit(desc,','),'[',1)) %>% 
  mutate(desc=substr(desc,1,20)) %>% 
  mutate(desc=sapply(strsplit(desc,' '),'[',1))

#custom fix some names

kd <- kd %>% 
  mutate(desc=ifelse(AGI=='AT2G20420','ATP citrate lyase',
                     ifelse(AGI=='AT1G54220','Dihydrolipoamide acetyltransferase',
                            ifelse(AGI=='AT3G55410','2-oxoglutarate dehydrogenase',
                                   ifelse(AGI=='AT4G02580','24 kDa subunit',
                                          ifelse(AGI=='AT5G20080','FAD/NAD(P)-binding oxidoreductase',
                                                 ifelse(AGI=='AT1G24360','NAD(P)-binding',
                                                        ifelse(AGI=='AT5G08670','ATP synthase alpha/beta family protein',
                                                               ifelse(AGI=='AT4G02930','GTP binding Elongation factor Tu',desc)))))))))
    
    

#levels
levels <- arrange(kd,`Average KD (d-1)`) %>% distinct(desc)

kd$desc <- factor(kd$desc, levels=levels$desc)


ggplot(kd, aes(desc,`Average KD (d-1)`))+
  geom_bar(stat='identity',fill='darkgreen')+
  coord_flip()+
  theme(axis.title.y = element_blank())+
  labs(title = 'predicted targets kd - Lei')



## plot CI kd from Lon1 paper lei
## mitochondria solution


kd <- fread('data/lei_degradation_rate_lon_mito.csv',skip=2) %>% 
  as.tibble() %>% 
  separate_rows(AGI, convert=T) %>% 
  mutate(AGI=substr(AGI,1,9)) %>% 
  distinct(AGI, .keep_all = T)


#filter for dep proteins
kd <- kd %>% filter(AGI %in% data_long$ID)



#intersect in Complex I
kd <- kd %>% left_join(cI_anno, by=c('AGI'='ACC')) %>% 
  na.omit()

#levels
levels <- arrange(kd,`WT KD (Average)`) %>% distinct(Name)

kd$Name <- factor(kd$Name, levels=levels$Name)

ggplot(kd, aes(Name,`WT KD (Average)`))+
  geom_bar(stat='identity',fill='darkgreen')+
  coord_flip()+
  theme(axis.title.y = element_blank())+
  labs(title = 'Complex I kd - Lei - lon paper')



## plot CI kd from Lon1 paper lei
## mitochondria BN fractions

#load data, convert to long format
kd <- fread('data/lei_lon1_bn_deg_rates.csv',header = F) %>% 
  as.tibble()
colnames(kd) <- paste(kd[1,],gsub(' ','_',kd[2,]),sep='-')

for (n in which(duplicated(colnames(kd))==T)){
  colnames(kd)[n]=paste('wt',colnames(kd)[n],sep='%')
}
kd <- kd[3:nrow(kd),]
colnames(kd)[1:5] <- gsub('-','',colnames(kd)[1:5])


kd <- kd %>% 
  separate_rows(AGI, convert=T) %>% 
  mutate(AGI=substr(AGI,1,9)) %>% 
  distinct(AGI, .keep_all = T)

kd <- kd %>% gather(sample,value,6:ncol(.))

#filter for average and wild type & split sample column into meta data
kd <- kd %>% 
  filter(str_detect(sample,regex('wt', ignore_case = T)) 
         & str_detect(sample,regex('average', ignore_case = T))) %>% 
  mutate(fraction=sapply(strsplit(sample,'_'),'[',2),
         mem_sol=sapply(strsplit(sapply(strsplit(sample,'-'),'[',2),'_'),'[',1))

#replace empty values with 0 (not found)

kd <- kd %>% 
  mutate(value=as.numeric(value),
         value=ifelse(is.na(value),0,value))
  


#filter for dep proteins
kd <- kd %>% filter(AGI %in% data_long$ID)


#intersect in Complex I
kd <- kd %>% left_join(cI_anno, by=c('AGI'='ACC')) %>% 
  na.omit()


##adjust soluble and insoluble fraction positions


kd <- kd %>% 
  mutate(fraction=as.numeric(gsub('F','',fraction)),
         fraction=ifelse(mem_sol=='Soluble',fraction+3,fraction),
         fraction=paste0('F',fraction))


#levels
kd <- kd %>%
  group_by(fraction,mem_sol,Name) %>% 
  mutate(max=max(value)) %>% 
  ungroup()
levels <- arrange(kd,mem_sol,max) %>% distinct(Name)

kd$Name <- factor(kd$Name, levels=levels$Name)
kd$fraction <- factor(kd$fraction, levels=rev(c('F1','F2','F3','F4','F5','F6','F7','F8',
                                            'F9','F10','F11','F12','F13','F14','F15')))

kd <- kd %>% 
  mutate(value=ifelse(value < 0, value==0.01,value))


#ccombined kd
#library(RColorBrewer)
p=ggplot(kd, aes(fraction,value,fill=mem_sol))+
  geom_bar(stat='identity',alpha=0.5)+
  coord_flip()+
  facet_wrap(~Name,scales='free_y')+
  theme(axis.title.y = element_blank(),axis.text.y = element_text(size=8,face = 'bold'))+
  labs(title = 'Complex I kd - Lei - lon paper - bn gel')+
  theme(legend.position = c(0.5,0.1),axis.title.x = element_blank())+
  ylim(0,0.4)+
  scale_fill_brewer(palette = 'Dark2')

ggsave(filename = paste0('Complex I kd - Lei - lon paper - bn gel','.svg'),path = 'images/kd',device = 'svg',dpi=1080,plot = p)



#add protein abundance

quant <- fread('data/lei_abundance_lon_mito.csv') 
colnames(quant) <- c(as.character(quant[2,1:3]),paste(quant[1,4:ncol(quant)],quant[2,4:ncol(quant)],sep='-'))
quant <- quant[3:nrow(quant),]  

quant <- quant %>% 
  as_tibble(.name_repair = 'universal') %>% 
  filter(Identified.Pep.NO..p.0.95.=='WT') %>% 
  dplyr::select(-2) %>% 
  gather(sample,count,3:ncol(.))

#fix AGI
quant <- quant %>% 
  separate_rows(Protein,convert = T) %>% 
  mutate(Protein=substr(Protein,1,9)) %>% 
  distinct(Protein,sample, .keep_all = T)

#match Kd frame for left join

quant <- quant %>% 
  mutate(fraction=sapply(strsplit(sample, '[.]'),'[',2),
         mem_sol=ifelse(str_detect(sample,'Insoluble')==T,'Insoluble','Soluble'),
         fraction=as.numeric(gsub('F','',fraction)),
         fraction=ifelse(mem_sol=='Soluble',fraction+3,fraction),
         fraction=paste0('F',fraction))

quant <- quant %>% 
  dplyr::rename(AGI=Protein) %>% 
  dplyr::select(-2,-3)

#join

kd <- kd %>% left_join(quant)

#ccombined quant
#library(RColorBrewer)


#heaps of eleves (taken from kd)
kd <- kd %>%
  group_by(fraction,mem_sol,Name) %>% 
  mutate(max=max(value)) %>% 
  ungroup()
levels <- arrange(kd,mem_sol,max) %>% distinct(Name)

kd$Name <- factor(kd$Name, levels=levels$Name)
kd$fraction <- factor(kd$fraction, levels=rev(c('F1','F2','F3','F4','F5','F6','F7','F8',
                                                'F9','F10','F11','F12','F13','F14','F15')))

kd <- kd %>% 
  mutate(value=ifelse(value < 0, value==0.01,value))


p=ggplot(kd, aes(fraction,as.numeric(count),fill=mem_sol))+
  geom_bar(stat='identity',alpha=0.5)+
  coord_flip()+
  facet_wrap(~Name,scales='free')+
  theme(axis.title.y = element_blank(),axis.text.y = element_text(size=8,face = 'bold'))+
  labs(title = 'Complex I count - Lei - lon paper - bn gel')+
  theme(legend.position = c(0.5,0.1),axis.title.x = element_blank())+
  #ylim(0,0.4)+
  scale_fill_brewer(palette = 'Dark2')

ggsave(filename = paste0('Complex I kd - Lei - lon paper - bn gel','.svg'),path = 'images/kd',device = 'svg',dpi=1080,plot = p)











################              Figure 1c clp paper               #########################################

##data prep

# # Generate a results table
# # no pvalue cutoff
# dep <- add_rejections(data_diff, alpha = 1, lfc = log2(1))
# 
# data_results <- get_results(dep)
# data_results <- data_results[,which(str_detect(colnames(data_results),'significant')==F)]

data_results <- fread('data/fig3a_data.csv')

#long format(gather)
data_long <- data_results %>% gather(sample,value,3:ncol(.)) %>% 
  mutate(type=sapply(strsplit(sample,'__'),'[',3),
         type=ifelse(is.na(type)==T,sapply(strsplit(sample,'__'),'[',2),type),
         sample=sapply(strsplit(sample,'__'),'[',1),
         fraction=ifelse(sample %in% c("clp1p","clp2p", "wtp"),'membrane','soluble'),
         genotype=substr(sample,1,nchar(sample)-1)) %>% 
  dplyr::select(-sample) %>%
  rowwise() %>% 
  spread(key = type,value = value) %>% 
  mutate(ID=sapply(strsplit(ID,'[.]'),'[',1)) %>% 
  filter(str_detect(genotype,'signif')==F)



#Gene selection
#clpp
clpp <- 'AT5G23140'
#most stabel house keeper overall
#https://link.springer.com/article/10.1007/s13562-017-0403-0
Lon1 <- 'AT5G26860'
#mito housekeeper vdac1
E1alpha <- 'AT1G59900'

#selction

sel <- c(clpp,Lon1,E1alpha)
sel <- dplyr::filter(data_long,ID %in% sel) %>% 
  distinct(name)


#Plot counts
#custom selection, Plotcounts() doesn't work on multiple Genes
subset <- assay(dep[sel$name]) %>% 
  as.data.frame() %>% 
  rownames_to_column(var='gene') %>% 
  as_tibble() %>% 
  gather(sample,count,2:24) %>% 
  mutate(genotype=ifelse(str_detect(sample,'clp1'),'clp1',
                         ifelse(str_detect(sample,'clp2'),'clp2','WT')),
         fraction=ifelse(str_detect(sample,'p__')==T,'membrane','soluble'),
         desc=ifelse(gene=='CLPP2','CLPP2',
                     ifelse(gene=='LON1','LON1',
                            ifelse(gene=='E1 ALPHA','PDC E1 \U221D','failsave'))),
         raw_count=2^count)




#little stars for pvalues :)


res <- data_long %>% 
  dplyr::filter(name %in% sel$name, fraction == 'soluble') %>% 
  dplyr::select(name,fraction,genotype,p.adj) %>% 
  mutate(p.adj=ifelse(genotype=='wt',1,p.adj),
         sig_level=ifelse(p.adj > 0.05,'',
                          ifelse(p.adj <= 0.05 & p.adj > 0.005,'*',
                                 ifelse(p.adj <= 0.005,'*\n*','failsave'))),
         genotype=ifelse(genotype=='wt','WT',genotype),
         name=ifelse(name=='E1 ALPHA','PDC E1 \U221D',name)) %>% 
  dplyr::rename(desc=name)


#add median as ref y axis point 
y_ref <- subset %>% group_by(genotype,desc) %>% 
  summarise(median=median(raw_count/1000000),max=max(raw_count/1000000))
res <- res %>% 
  left_join(y_ref)

#levels

subset$genotype <- factor(subset$genotype, levels = c('WT','clp1','clp2'))
subset$desc <- factor(subset$desc, levels = c('CLPP2','LON1','PDC E1 \U221D'))
res$genotype <- factor(res$genotype, levels = c('WT','clp1','clp2'))
res$desc <- factor(res$desc, levels = c('CLPP2','LON1','PDC E1 \U221D'))





#plot
g <- ggplot(dplyr::filter(subset,fraction=='soluble'), aes(genotype, raw_count/1000000, bg=genotype)) + 
  facet_wrap(~desc,scales = 'free_y')+
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar",col='black', size = 0.3)+
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "bar",col='black', size = 0.15,alpha= 0.6)+
  geom_point(pch = 21,size=2,color='black',alpha=0.5)+
  geom_text(data=res,aes(genotype,c(0.65,0.6,rep(1,7)),label=sig_level),size=6, lineheight = 0.25)+
  expand_limits(y=0)+
  scale_colour_manual(values=c('#339900','#3399cc','#3366cc'))+
  scale_fill_manual(values=c('#339900','#3399cc','#3366cc'))+
  labs(title='Mitochondrial protein abundance',y='LFQ intensity [M]')+
  theme(axis.title.x = element_blank(),legend.position = 'none',
        axis.text.x = element_text(face=c('plain','italic','italic'),size=8, angle = 30),
        axis.title.y = element_text(face='bold',size='8'),
        axis.text.y = element_text(face='bold',size=8),
        strip.text = element_text(face='bold',size=8),
        title=element_text(size=10))



#save
ggsave('Prot_KO_figure1.pdf',device = 'pdf',dpi=1080,plot = g,height = 6.52,width = 6,units = 'cm')



################              Figure 3a clp paper               #########################################

library(data.table)
library(tidyverse)



## custom mutant overlap volcano


# #only run if pre existing data is not good
# #data changed to p=0.05
# dep <- add_rejections(data_diff, alpha = 0.05, lfc = log2(1))
# data_results <- get_results(dep)
# data_results %>% filter(significant) %>% nrow()

#write table s dep seems to perform differently every time
#write_csv(data_results,'data/fig3a_data.csv')


#read data rom premade DEP
data_results <- fread('data/fig3a_data.csv')




#long format(gather)
data_long <- data_results %>% gather(sample,value,3:ncol(.)) %>% 
  mutate(type=sapply(strsplit(sample,'__'),'[',3),
         type=ifelse(is.na(type)==T,sapply(strsplit(sample,'__'),'[',2),type),
         sample=sapply(strsplit(sample,'__'),'[',1),
         fraction=ifelse(sample %in% c("clp1p","clp2p", "wtp"),'membrane','soluble'),
         genotype=substr(sample,1,nchar(sample)-1)) %>% 
  dplyr::select(-sample) %>%
  rowwise() %>% 
  filter(type != 'significant') %>% 
  spread(key = type,value = value) %>% 
  mutate(ID=sapply(strsplit(ID,'[.]'),'[',1))


#filter for both mutants under 0.1pvalue, and use lowest FC and highest pvalue
data_volcano <-  
  data_long %>% 
  filter(genotype !='wt') %>% 
  group_by(ID,fraction) %>% 
  mutate(colour_p=ifelse(max(p.adj) <= 0.05,'red','blue')) %>% 
  mutate(min_ratio=min(abs(ratio)),
         max_p=max(p.adj)) %>% 
  filter(abs(ratio)==min_ratio) %>% 
  dplyr::select(-min_ratio) %>% 
  mutate(colour_r=ifelse(ratio <=-0.4 | ratio >= 0.4,'red','blue')) %>% 
  mutate(sig=ifelse(colour_p=='blue'|colour_r=='blue','non_sig','sig')) %>% 
  distinct(ID,fraction, .keep_all = T)

# #add good annotation
# desc_volcano <- data2 %>% 
#   dplyr::select(Gene.names,`Majority protein IDs`) %>% 
#   dplyr::rename(desc=Gene.names,ACC=`Majority protein IDs`)

# #hand curate some annotation
# 
# desc_volcano <- desc_volcano %>% 
#   mutate(desc=ifelse(ACC=='AT5G08690','ATP synth beta2',
#                      ifelse(ACC=='AT5G62270','mucin related AT5G62270',
#                             ifelse(ACC=='AT1G26460','PPR AT1G26460',
#                                    ifelse(ACC=='AT3G02650','PPR AT3G02650',
#                                           ifelse(ACC=='AT5G64670','PPR AT5G64670',desc))))))
# 
# desc_volcano <- desc_volcano %>% 
#   mutate(desc=ifelse(ACC=='AT3G62530','ARM repeat',
#                      ifelse(ACC=='AT4G21020','Late embryo abundant',
#                             ifelse(ACC=='AT5G55200','CoChaperone GrpE',
#                                    ifelse(ACC=='AT3G19740','p-loop hydrolase',
#                                           ifelse(ACC=='AT4G36680','PPR AT4G36680',
#                                                  ifelse(ACC=='AT1G70190','Ribo L7/L12',desc)))))))
# 
# #write table, hand curate and reimport
# #filter desk for significant proteins
# 
#sig_ids <- filter(data_volcano, sig=='sig') %>% ungroup() %>% distinct(ID)
# 
# 
#write_csv(filter(desc_volcano, ACC %in% sig_ids$ID),'data/desc_volcano3.csv')
desc_volcano <- read_csv('data/desc_volcano3.csv')


#merge annnotation into volcano data

data_volcano <- data_volcano %>% 
  left_join(desc_volcano, by=c('ID'='ACC')) %>% 
  mutate(desc=ifelse(is.na(desc)==T | sig=='non_sig','',desc),
         group=ifelse(is.na(group)==T | sig=='non_sig','',group))


#levels
data_volcano$sig <- factor(data_volcano$sig, levels=c('sig','non_sig'))

## format data for plot
#change-log10 pvalue  high numbers to a maximum of 7.5
#change symbol of  >7.5 to op arrow

data_volcano <- data_volcano %>% 
  mutate(max_p_adj=ifelse(-log10(max_p) > 6,1e-6,max_p),
         pch=ifelse(-log10(max_p) > 6,24,21),
         ratio_adj=ifelse(ratio > 3.5,3.5,ifelse(ratio < -3.5,-3.5,ratio)))


#scale colour based on group
data_volcano <- data_volcano %>% 
  mutate(group=ifelse(sig=='non_sig','non_sig',group),
         alpha=ifelse(group=='non_sig',0.4,0.75))
data_volcano$group <- factor(data_volcano$group , levels=c('non_sig',
                                                           'Energy Metabolism',
                                                           'Mitochondrial Protein Synthesis',
                                                           'Chaperone / Protease',
                                                           'Knock Out',
                                                           'Ungrouped'))

#rename facets
data_volcano <-data_volcano %>% ungroup() %>% 
  mutate(fraction=ifelse(fraction=='membrane','Membrane fraction','Soluble fraction'))


library(ggrepel)

p <- ggplot(data_volcano, aes(x=ratio_adj,y=-log10(max_p_adj),fill=group))+
  geom_point(pch=data_volcano$pch,alpha=data_volcano$alpha,size=3)+
  geom_text_repel(data=filter(data_volcano,sig=='sig'),aes(label=table_nr),col='black',size=2.5, fontface='bold')+
  geom_hline(yintercept = -log10(0.05),size=0.3, alpha=0.5,linetype="dashed",color='#0033ff')+
  geom_vline(xintercept = c(-0.4,0.4),size=0.3, alpha=0.5,linetype="dashed",color='#0033ff')+
  geom_text(aes(x=3.5,y=-log10(0.05)-0.1,label='P = 0.05'), size = 2.5, colour='#003366')+
  geom_text(aes(x=-0.4-0.5,y=5.5,label='Log2FC \u2264 -0.4'), size = 2.5, colour='#003366')+
  geom_text(aes(x=0.4+0.5,y=5.5,label='Log2FC \u2265 -0.4'), size = 2.5, colour='#003366')+
  facet_wrap(~fraction)+
  scale_fill_manual(values=c('#006699','#FF3300','#cc9966','#006600','#660000','#9933ff'))+
  theme(legend.position = 'none', axis.title = element_text(face='bold',size = 18),
        axis.text = element_text(face='bold',size = 16), strip.text = element_text(face='bold',size=18),
        title = element_text(face='bold',size=18))+
  labs(title='', x='log2FC clp1/clp2 VS WT', y=expression(paste("-Lo", g[10]," P",sep="")))
p



#save
ggsave('Prot_volcano_figure3a.pdf',device = 'pdf',dpi=2160,plot = p,height = 8,width = 13,units = 'in')

################              Figure 3b clp paper   (table)            ##################################


#start with desc volcano
fig3b <- desc_volcano %>% 
  left_join(dplyr::select(data_long,2,3,4,6,8), by=c('ACC'='ID')) %>% 
  distinct() %>% 
  dplyr::filter(genotype !='wt') %>% 
  dplyr::rename(AGI=ACC)
  
fig3b <- fig3b[,c(5,4,1,2,3,6,7,8)]

fig3b <- fig3b %>% 
  mutate(p_r=paste(p.adj,ratio,sep = '_')) %>% 
  dplyr::select(-p.adj,-ratio) %>% 
  spread(genotype,p_r) %>%
  mutate(`clp1 ratio`=round(as.numeric(sapply(strsplit(clp1,'_'),'[',2)),digits=2),
         `clp2 ratio`=round(as.numeric(sapply(strsplit(clp2,'_'),'[',2)),digits=2),
         `clp1 p.adj`=round(as.numeric(sapply(strsplit(clp1,'_'),'[',1)),digits=2),
         `clp2 p.adj`=round(as.numeric(sapply(strsplit(clp2,'_'),'[',1)),digits=2)) %>% 
  dplyr::select(-clp1,-clp2)


write.csv(fig3b,'data/fig3b.csv')
  





################              Figure 5a clp paper               #########################################

##data prep

# Generate a results table
# no pvalue cutoff
data_results <- fread('data/fig3a_data.csv')


#long format(gather)
data_long <- data_results %>% gather(sample,value,3:ncol(.)) %>% 
  mutate(type=sapply(strsplit(sample,'__'),'[',3),
         type=ifelse(is.na(type)==T,sapply(strsplit(sample,'__'),'[',2),type),
         sample=sapply(strsplit(sample,'__'),'[',1),
         fraction=ifelse(sample %in% c("clp1p","clp2p", "wtp"),'membrane','soluble'),
         genotype=substr(sample,1,nchar(sample)-1)) %>% 
  dplyr::select(-sample) %>%
  rowwise() %>% 
  spread(key = type,value = value) %>% 
  mutate(ID=sapply(strsplit(ID,'[.]'),'[',1),
         name=gsub('Â','',name)) %>% 
  dplyr::filter(str_detect(genotype,'signif')==F)


#Gene selection
#B8
B8 <- 'AT5G47890'
#B14
B14 <- 'AT3G12260'
#24KDa
kd24 <- 'AT4G02580'
#51KDa
kd51 <- 'AT5G08530'
#75KDa
kd75 <- 'AT5G37510'

#selction

sel <- c(B8,B14,kd24,kd51,kd75)
sel <- dplyr::filter(data_long,ID %in% sel) %>% 
  distinct(name)

sel$name <- c('B8','B14',rownames(dep)[c(47,58,61)])


#Plot counts
#custom selection, Plotcounts() doesn't work on multiple Genes
subset <- assay(dep[sel$name]) %>% 
  as.data.frame() %>% 
  rownames_to_column(var='gene') %>% 
  as_tibble() %>% 
  gather(sample,count,2:24) %>% 
  mutate(genotype=ifelse(str_detect(sample,'clp1'),'clp1',
                         ifelse(str_detect(sample,'clp2'),'clp2','WT')),
         desc=gene,
         fraction=ifelse(str_detect(sample,'p__')==T,'membrane','soluble'),
         raw_count=2^count)


#stats
library(broom)
all_ttest <- data.frame()
for (prot in unique(subset$gene)){
  stat=filter(subset,gene==prot) 
  stat=tidy(pairwise.t.test(stat$raw_count,stat$genotype,p.adjust.method = 'BH')) %>%
    as.data.frame() %>% 
    mutate(p.value=round(p.value,digits=3),
           gene=prot) %>% 
    filter(group1=='WT') %>% 
    dplyr::select(-group1) %>% 
    dplyr::rename(genotype=group2)
  all_ttest=bind_rows(all_ttest,stat)
  
}
subset <- left_join(subset,all_ttest) %>% 
  mutate(p.value=ifelse(is.na(p.value),'1',p.value),
         sig=ifelse(p.value > 0.05,'',
                    ifelse(p.value <= 0.05 & p.value > 0.001,'*',
                           ifelse(p.value <= 0.001,'*\n*','failsave'))))


#stat_data

stat_data <- subset %>% 
  distinct(gene,genotype,p.value,.keep_all = T)


subset$genotype <- factor(subset$genotype, levels = c('WT','clp1','clp2'))
subset$desc <- factor(subset$desc, levels = c('B8','B14',unique(subset$gene)[3:5]))



stat_data$genotype <- factor(stat_data$genotype, levels = c('WT','clp1','clp2'))
stat_data$desc <- factor(stat_data$desc, levels = c('B8','B14',unique(stat_data$desc)[3:5]))
stat_data$gene <- factor(stat_data$gene, levels = c('B8','B14',unique(stat_data$gene)[3:5]))



#plot
g <- ggplot(dplyr::filter(subset), aes(genotype, raw_count/1000000, bg=genotype)) + 
  facet_wrap(~desc,scales = 'free_y',nrow = 1)+
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "bar",col='black', size = 0.15,alpha= 0.6,,width=0.8)+
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar",col='black', size = 0.3,width=0.8)+
  geom_jitter(pch = 21,size=1,color='black',alpha=0.2,width=0.15)+
  geom_text(data=stat_data,aes(genotype,c(1.6,6.5,8.5,21.5,28.5,1.7,8,9.5,23,28,1,1,1,1,1),label=sig),size=6, lineheight = 0.25)+
  expand_limits(y=0)+
  scale_colour_manual(values=c('#339900','#3399cc','#3366cc'))+
  scale_fill_manual(values=c('#339900','#3399cc','#3366cc'))+
  labs(title='Complex I protein abundance',y='LFQ intensity [M]')+
  theme(axis.title.x = element_blank(),legend.position = 'none',
        axis.text.x = element_text(face=c('plain','italic','italic'),size=8, angle = 90),
        axis.title.y = element_text(face='bold',size='8'),
        axis.text.y = element_text(face='bold',size=8),
        strip.text = element_text(face='bold',size=8),
        title=element_text(size=10))
g


#save
ggsave('Prot_CIabundacne_overall_figure4a.pdf',device = 'pdf',dpi=1080,plot = g,height = 6,width = 8,units = 'cm')



#overlap shotguna nd chafradic


data_results <- data_results %>% 
  mutate(gene=substr(ID,1,9))

chafradic <- fread('data/n_terms.csv')

combine <- left_join(chafradic,data_results, by='gene') %>% na.omit()


write.csv(combine,'data/combined_shotgun_chafradic.csv')





################              Figure 5 clp paper               #########################################



#combine prot volcano data with RNA seq


library(data.table)
library(tidyverse)



## custom mutant overlap volcano


# #only run if pre existing data is not good
# #data changed to p=0.05
# dep <- add_rejections(data_diff, alpha = 0.05, lfc = log2(1))
# data_results <- get_results(dep)
# data_results %>% filter(significant) %>% nrow()

#write table s dep seems to perform differently every time
#write_csv(data_results,'data/fig3a_data.csv')



#read data rom premade DEP
data_results <- fread('data/fig3a_data.csv')




#long format(gather)
data_long <- data_results %>% gather(sample,value,3:ncol(.)) %>% 
  mutate(type=sapply(strsplit(sample,'__'),'[',3),
         type=ifelse(is.na(type)==T,sapply(strsplit(sample,'__'),'[',2),type),
         sample=sapply(strsplit(sample,'__'),'[',1),
         fraction=ifelse(sample %in% c("clp1p","clp2p", "wtp"),'membrane','soluble'),
         genotype=substr(sample,1,nchar(sample)-1)) %>% 
  dplyr::select(-sample) %>%
  rowwise() %>% 
  dplyr::filter(type != 'significant') %>% 
  spread(key = type,value = value) %>% 
  mutate(ID=sapply(strsplit(ID,'[.]'),'[',1))


#filter for both mutants under 0.1pvalue, and use lowest FC and highest pvalue
data_volcano <-  
  data_long %>% 
  dplyr::filter(genotype !='wt') %>% 
  group_by(ID,fraction) %>% 
  mutate(colour_p=ifelse(max(p.adj) <= 0.05,'red','blue')) %>% 
  mutate(min_ratio=min(abs(ratio)),
         max_p=max(p.adj)) %>% 
  dplyr::filter(abs(ratio)==min_ratio) %>% 
  dplyr::select(-min_ratio) %>% 
  mutate(colour_r=ifelse(ratio <=-0.4 | ratio >= 0.4,'red','blue')) %>% 
  mutate(sig=ifelse(colour_p=='blue'|colour_r=='blue','non_sig','sig')) %>% 
  distinct(ID,fraction, .keep_all = T)

# #add good annotation
# desc_volcano <- data2 %>% 
#   dplyr::select(Gene.names,`Majority protein IDs`) %>% 
#   dplyr::rename(desc=Gene.names,ACC=`Majority protein IDs`)

# #hand curate some annotation
# 
# desc_volcano <- desc_volcano %>% 
#   mutate(desc=ifelse(ACC=='AT5G08690','ATP synth beta2',
#                      ifelse(ACC=='AT5G62270','mucin related AT5G62270',
#                             ifelse(ACC=='AT1G26460','PPR AT1G26460',
#                                    ifelse(ACC=='AT3G02650','PPR AT3G02650',
#                                           ifelse(ACC=='AT5G64670','PPR AT5G64670',desc))))))
# 
# desc_volcano <- desc_volcano %>% 
#   mutate(desc=ifelse(ACC=='AT3G62530','ARM repeat',
#                      ifelse(ACC=='AT4G21020','Late embryo abundant',
#                             ifelse(ACC=='AT5G55200','CoChaperone GrpE',
#                                    ifelse(ACC=='AT3G19740','p-loop hydrolase',
#                                           ifelse(ACC=='AT4G36680','PPR AT4G36680',
#                                                  ifelse(ACC=='AT1G70190','Ribo L7/L12',desc)))))))
# 
# #write table, hand curate and reimport
# #filter desk for significant proteins
# 
#sig_ids <- filter(data_volcano, sig=='sig') %>% ungroup() %>% distinct(ID)
# 
# 
#write_csv(filter(desc_volcano, ACC %in% sig_ids$ID),'data/desc_volcano3.csv')
desc_volcano <- read_csv('data/desc_volcano3.csv')


#merge annnotation into volcano data

data_volcano <- data_volcano %>% 
  left_join(desc_volcano, by=c('ID'='ACC')) %>% 
  mutate(desc=ifelse(is.na(desc)==T | sig=='non_sig','',desc),
         group=ifelse(is.na(group)==T | sig=='non_sig','',group))


#levels
data_volcano$sig <- factor(data_volcano$sig, levels=c('sig','non_sig'))

## format data for plot
#change-log10 pvalue  high numbers to a maximum of 7.5
#change symbol of  >7.5 to op arrow

data_volcano <- data_volcano %>% 
  mutate(max_p_adj=ifelse(-log10(max_p) > 6,1e-6,max_p),
         pch=ifelse(-log10(max_p) > 6,24,21),
         ratio_adj=ifelse(ratio > 3.5,3.5,ifelse(ratio < -3.5,-3.5,ratio)))


#scale colour based on group
data_volcano <- data_volcano %>% 
  mutate(group=ifelse(sig=='non_sig','non_sig',group),
         alpha=ifelse(group=='non_sig',0.4,0.75))
data_volcano$group <- factor(data_volcano$group , levels=c('non_sig',
                                                           'Energy Metabolism',
                                                           'Mitochondrial Protein Synthesis',
                                                           'Chaperone / Protease',
                                                           'Knock Out',
                                                           'Ungrouped'))

#rename facets
data_volcano <-data_volcano %>% ungroup() %>% 
  mutate(fraction=ifelse(fraction=='membrane','Membrane fraction','Soluble fraction'))


data_volcano_prot <- data_volcano
data_volcano_rna <- fread('data/data_voclano_rnaseq.csv') %>% 
  dplyr::select(rowname,ratio_adj)

combined <- data_volcano_prot %>% 
  left_join(data_volcano_rna,by=c('ID'='rowname'))

#shrink ratio.y

combined <- combined %>% 
  mutate(ratio_adj.y=ifelse(ratio_adj.y >= 1,1,
                            ifelse(ratio_adj.y <= -1,-1,ratio_adj.y)))


g <- ggplot(dplyr::filter(combined,sig=='sig'|abs(ratio_adj.x) >= 1 & max_p_adj <= 0.05|abs(ratio_adj.y) >= 0.5 ),aes(ratio_adj.x,ratio_adj.y))+
  geom_point()+
  geom_text_repel(data=dplyr::filter(combined,sig=='sig'|abs(ratio_adj.x) >= 1 & max_p_adj <= 0.05|abs(ratio_adj.y) >= 0.5 ),aes(label=name),col='black',size=2.5, fontface='bold')



  
ggsave('rna_prot_combine.pdf',device = 'pdf',dpi=1080,plot = g,height = 15,width = 15,units = 'cm')



################              Figure mito encodes expression sublemental clp paper               #########################################
#########################################################################################################

##data prep

# # Generate a results table
# # no pvalue cutoff
# dep <- add_rejections(data_diff, alpha = 1, lfc = log2(1))
# 
# data_results <- get_results(dep)
# data_results <- data_results[,which(str_detect(colnames(data_results),'significant')==F)]

data_results <- fread('data/fig3a_data.csv')

#long format(gather)
data_long <- data_results %>% gather(sample,value,3:ncol(.)) %>% 
  mutate(type=sapply(strsplit(sample,'__'),'[',3),
         type=ifelse(is.na(type)==T,sapply(strsplit(sample,'__'),'[',2),type),
         sample=sapply(strsplit(sample,'__'),'[',1),
         fraction=ifelse(sample %in% c("clp1p","clp2p", "wtp"),'membrane','soluble'),
         genotype=substr(sample,1,nchar(sample)-1)) %>% 
  dplyr::select(-sample) %>%
  rowwise() %>% 
  spread(key = type,value = value) %>% 
  mutate(ID=sapply(strsplit(ID,'[.]'),'[',1)) %>% 
  dplyr::filter(str_detect(genotype,'signif')==F)



#Gene selection
#everything ATM

#selction

sel <- dplyr::filter(data_long,str_detect(ID,'ATM')==T) %>% 
  distinct(name)


#Plot counts
#custom selection, Plotcounts() doesn't work on multiple Genes
subset <- assay(dep[sel$name]) %>% 
  as.data.frame() %>% 
  rownames_to_column(var='gene') %>% 
  as_tibble() %>% 
  gather(sample,count,2:24) %>% 
  mutate(genotype=ifelse(str_detect(sample,'clp1'),'clp1',
                         ifelse(str_detect(sample,'clp2'),'clp2','WT')),
         fraction=ifelse(str_detect(sample,'p__')==T,'membrane','soluble'),
         desc=ifelse(gene=='CLPP2','CLPP2',
                     ifelse(gene=='LON1','LON1',
                            ifelse(gene=='E1 ALPHA','PDC E1 \U221D','failsave'))),
         raw_count=2^count)




#little stars for pvalues :)


res <- data_long %>% 
  dplyr::filter(name %in% sel$name, fraction == 'soluble') %>% 
  dplyr::select(name,fraction,genotype,p.adj) %>% 
  mutate(p.adj=ifelse(genotype=='wt',1,p.adj),
         sig_level=ifelse(p.adj > 0.05,'',
                          ifelse(p.adj <= 0.05 & p.adj > 0.005,'*',
                                 ifelse(p.adj <= 0.005,'*\n*','failsave'))),
         genotype=ifelse(genotype=='wt','WT',genotype),
         name=ifelse(name=='E1 ALPHA','PDC E1 \U221D',name)) %>% 
  dplyr::rename(desc=name)


#add median as ref y axis point 
y_ref <- subset %>% group_by(genotype,desc) %>% 
  summarise(median=median(raw_count/1000000),max=max(raw_count/1000000))
res <- res %>% 
  left_join(y_ref)

#levels

subset$genotype <- factor(subset$genotype, levels = c('WT','clp1','clp2'))
subset$desc <- factor(subset$desc, levels = c('CLPP2','LON1','PDC E1 \U221D'))
res$genotype <- factor(res$genotype, levels = c('WT','clp1','clp2'))
res$desc <- factor(res$desc, levels = c('CLPP2','LON1','PDC E1 \U221D'))





#plot
g <- ggplot(dplyr::filter(subset,fraction=='soluble'), aes(genotype, raw_count/1000000, bg=genotype)) + 
  facet_wrap(~desc,scales = 'free_y')+
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar",col='black', size = 0.3)+
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "bar",col='black', size = 0.15,alpha= 0.6)+
  geom_point(pch = 21,size=2,color='black',alpha=0.5)+
  geom_text(data=res,aes(genotype,c(0.65,0.6,rep(1,7)),label=sig_level),size=6, lineheight = 0.25)+
  expand_limits(y=0)+
  scale_colour_manual(values=c('#339900','#3399cc','#3366cc'))+
  scale_fill_manual(values=c('#339900','#3399cc','#3366cc'))+
  labs(title='Mitochondrial protein abundance',y='LFQ intensity [M]')+
  theme(axis.title.x = element_blank(),legend.position = 'none',
        axis.text.x = element_text(face=c('plain','italic','italic'),size=8, angle = 30),
        axis.title.y = element_text(face='bold',size='8'),
        axis.text.y = element_text(face='bold',size=8),
        strip.text = element_text(face='bold',size=8),
        title=element_text(size=10))



#save
ggsave('Prot_KO_figure1.pdf',device = 'pdf',dpi=1080,plot = g,height = 6.52,width = 6,units = 'cm')


