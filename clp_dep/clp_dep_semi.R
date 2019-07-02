library(DEP)
library(tidyverse)
library(data.table)
library(fitdistrplus)
library(Biostrings)
library(SummarizedExperiment)

data <- fread('data/proteinGroups_semi.txt')
data <- filter(data, Reverse != "+", `Potential contaminant` != "+")
colnames(data)

##############################
##### Data preparation #######
##############################



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

#data3 <- data3[,c(1:3,221,4:220)] 
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




####################################################
##### Generate a SummarizedExperiment object #######
####################################################


# Generate a SummarizedExperiment object by parsing condition information from the column names
LFQ_columns <- grep("LFQ.", colnames(data_unique)) # get LFQ column numbers
data_se <- make_se_parse(data_unique, LFQ_columns)
data_se


#normal distribution ?
#before log2
data_dist <- data_unique[,which(str_detect(colnames(data_unique),'LFQ')==T)] %>% 
  gather(column,value)
descdist(data_dist$value)
hist(data_dist$value,col = 'royalblue',breaks = 500,xlim = c(0,50000000))
#after
data_dist_after <- assay(data_se) %>% as.data.frame() %>% 
  gather(column,value) %>% 
  na.omit()
descdist(data_dist_after$value)
hist(data_dist_after$value, col='royalblue')

######################################
##### Filter on missing values #######
######################################

plot_frequency(data_se)

# Filter for proteins that are identified in all replicates of at least one condition
data_filt <- filter_missval(data_se, thr = 0)

# Less stringent filtering:
# Filter for proteins that are identified in 2 out of 3 replicates of at least one condition
data_filt2 <- filter_missval(data_se, thr = 1)
data_filt <- data_filt2
# Plot a barplot of the number of identified proteins per samples
plot_numbers(data_filt)
# Plot a barplot of the protein identification overlap between samples
plot_coverage(data_filt)
# Normalize the data
data_norm <- normalize_vsn(data_filt)
# Visualize normalization by boxplots for all samples before and after normalization
plot_normalization(data_filt, data_norm)

############################################
##### Impute data for missing values #######
############################################

# Plot a heatmap of proteins with missing values
plot_missval(data_filt)
# Plot intensity distributions and cumulative fraction of proteins with and without missing values
plot_detect(data_filt)

# All possible imputation methods are printed in an error, if an invalid function name is given.
impute(data_norm, fun = "")
## Error in match.arg(fun): 'arg' should be one of "bpca", "knn", "QRILC", "MLE", "MinDet", "MinProb", "man", "min", "zero", "mixed", "nbavg"
# Impute missing data using random draws from a Gaussian distribution centered around a minimal value (for MNAR)
data_imp <- impute(data_norm, fun = "MinProb", q = 0.01)

# Impute missing data using random draws from a manually defined left-shifted Gaussian distribution (for MNAR)
data_imp_man <- impute(data_norm, fun = "man", shift = 1.8, scale = 0.3)
data_imp <- data_imp_man
# Impute missing data using the k-nearest neighbour approach (for MAR)
data_imp_knn <- impute(data_norm, fun = "knn", rowmax = 0.9)
data_imp <- data_imp_knn

# Plot intensity distributions before and after imputation
plot_imputation(data_norm, data_imp)

##############################################
##### Differential enrichment analysis #######
##############################################

# Differential enrichment analysis  based on linear models and empherical Bayes statistics

# Test every sample versus control
data_diff <- test_diff(data_imp, type = "control", control = "wt_")
## Tested contrasts: Ubi4_vs_Ctrl, Ubi6_vs_Ctrl, Ubi1_vs_Ctrl
# Test all possible comparisons of samples
data_diff_all_contrasts <- test_diff(data_imp, type = "all")
## Tested contrasts: Ubi4_vs_Ubi6, Ubi4_vs_Ctrl, Ubi4_vs_Ubi1, Ubi6_vs_Ctrl, Ubi6_vs_Ubi1, Ctrl_vs_Ubi1
# Test manually defined comparisons
data_diff_manual <- test_diff(data_imp, type = "manual", 
                              test = c("clp1s__vs_wts_", "clp2s__vs_wts_",'clp1p__vs_wtp_','clp2p__vs_wtp_'))
data_diff <- data_diff_manual
## Tested contrasts: Ubi4_vs_Ctrl, Ubi6_vs_Ctrl

# Denote significant proteins based on user defined cutoffs
dep <- add_rejections(data_diff, alpha = 0.1, lfc = log2(1))

######################
##### PCA PLOT #######
######################

# Plot the first and second principal components
plot_pca(dep, x = 1, n=nrow(dep),y = 2, point_size = 4)


########################
##### Cor Matrix #######
########################

# Plot the Pearson correlation matrix
plot_cor(dep, significant = TRUE, lower = 0, upper = 1, pal = "Reds")

#####################
##### Heatmap #######
#####################

# Plot a heatmap of all significant proteins with the data centered per protein
plot_heatmap(dep, type = "centered", kmeans = TRUE, 
             k = 6, col_limit = 4, show_row_names = F,
             indicate = c("condition", "replicate"))

# Plot a heatmap of all significant proteins (rows) and the tested contrasts (columns)
plot_heatmap(dep, type = "contrast", kmeans = TRUE, 
             k = 6, col_limit = 10, show_row_names = FALSE)


##########################
##### Volcano Plot #######
##########################

# Plot a volcano plot for the contrast "Ubi6 vs Ctrl""
plot_volcano(dep, contrast = "clp1s__vs_wts_", label_size = 3, add_names = TRUE)
plot_volcano(dep, contrast = "clp2s__vs_wts_", label_size = 3, add_names = TRUE)

###############################################
##### Barplots of a protein of interest #######
###############################################

#change levels
dep_save <- dep
dep@colData$condition <- gsub('_', '',dep@colData$condition)
dep@colData$condition <- factor(dep@colData$condition, 
                                levels = c('wtp','clp1p','clp2p','wts','clp1s','clp2s',4))
# Plot a barplot for USP15 and IKBKG
plot_single(dep, proteins = ' presequence protease 1 ')
# Plot a barplot for the protein USP15 with the data centered
plot_single(dep, proteins = ' lon protease 1 ', type = "centered")



###########################
##### Results table #######
###########################


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

prot_of_interest <- 'AT5G23140'

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
  labs(title=paste(prot_of_interest,rowname, sep=' - '),x = "sample", y = expression(log[2] ~ "Centered intensity" ~ "(±95% CI)"), col = "Rep") + 
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
#labs(title='LFQ Raw intensities CLP',x = "", y = ("LFQ Intensity" ~ "(±95% CI) [Million Area Count]"), col = "Rep") + 
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
    labs(title=chart,x = "", y = ("LFQ Intensity" ~ "(±95% CI) [Million Area Count]"), col = "Rep") + 
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

library(ggrepel)

p <- ggplot(data_volcano, aes(x=ratio,y=-log10(max_p),col=sig))+
  geom_point(pch=18,alpha=0.75)+
  geom_text_repel(data=filter(data_volcano,sig=='sig'),aes(label=desc),col='black',size=2.5, fontface='bold')+
  facet_wrap(~fraction)+
  scale_colour_manual(values=c('#990000','#99ccff'))+
  theme(legend.position = 'none')+
  labs(title='Clp vs WT', x='log2 fold change', y='-log10 p-value')

ggsave(filename = paste0('volcano','.png'),path = 'images',device = 'png',dpi=1080,plot = p)



###custom table significant genes
dep <- add_rejections(data_diff, alpha = 0.1, lfc = log2(1))

data_results <- get_results(dep)
data_results <- data_results[,which(str_detect(colnames(data_results),'significant')==F)]
data_results <- data_results[,which(str_detect(colnames(data_results),'p.val')==F)]
data_results <- data_results[,which(str_detect(colnames(data_results),'center')==F)]

data_results <- data_results[,c(2,1,3,5,7,9,4,6,8,10)]

write_csv(data_results,'data/results_clp.csv')

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
  labs(title='Complex I subunits',x = "", y = ("LFQ Intensity" ~ "(±95% CI) [Million Area Count]"), col = "Rep") +
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
  labs(title='Mitochondrial encoded Proteins',x = "", y = ("LFQ Intensity" ~ "(±95% CI) [Million Area Count]"), col = "Rep") +
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
  labs(title='LFQ Raw intensities CLP',x = "", y = ("LFQ Intensity" ~ "(±95% CI) [Million Area Count]"), col = "Rep") +
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
    labs(title=chart,x = "", y = ("LFQ Intensity" ~ "(±95% CI) [Million Area Count]"), col = "Rep") + 
    scale_fill_manual(values=c('#339900','#0066cc','#990033'))+
    scale_y_continuous(label=number_format(suffix=' M'))+
    theme_DEP2()+
    theme(axis.text.x = element_text(angle=0,hjust=0.5))
  ggsave(filename = paste0(filename,'.png'),path = 'images',device = 'png',dpi=1080,plot = p)
  message(paste0(filename,' done!'))
}



