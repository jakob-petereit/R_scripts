library(data.table)
library(tidyverse)
library(Biostrings)
library(stringi)
library(msa)
library(progress)
library(cowplot)


peptides <- fread('data/peptides_semi.txt')

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


### fix accession

peptides <- peptides %>% 
  filter(str_detect(`Leading razor protein`, 'AT')==T) %>% 
  mutate(ACC=substr(`Leading razor protein`,1,9)) 

#remove unnecessary columns
peptides <- peptides[,which(str_detect(colnames(peptides),'Count')==F)]
peptides <- peptides[,which(str_detect(colnames(peptides),'Experiment')==F)]
peptides <- peptides[,which(str_detect(colnames(peptides),'Intensity')==F)]

# #filter for complexI
# 
# CI_peptides <- filter(peptides,ACC %in% cI_anno$ACC)
# #add description
# CI_peptides <- CI_peptides %>% 
#   left_join(cI_anno)

#gather

peptides <- peptides %>% 
  gather(sample,intensity,`LFQ intensity clp1p_1`:`LFQ intensity wts_4`)

#remove more columns

peptides <- peptides %>% 
  select(-c(Mass,Proteins,`Leading razor protein`,`Mod. peptide IDs`,`Fraction Average`,`Fraction Std. Dev.`,
            Reverse,`Potential contaminant`,`Evidence IDs`,`MS/MS IDs`,`LFQ intensity blank`))
#add annotation

peptides <- peptides %>% 
  left_join(all_anno)

#reorder for visibility

peptides <- peptides[,c(24:27,1:23)]

#read sample meta data from sample name

peptides <- peptides %>% 
  mutate(sample=sapply(strsplit(sample,' '),'[',3)) %>% 
  filter(str_detect(sample,'pool')==F) %>% 
  mutate(replicate=sapply(strsplit(sample,'_'),'[',2),
         genotype=sapply(strsplit(sample,'._'),'[',1),
         fraction=ifelse(str_detect(sample,'p_')==T,'p','s')) 


#remove peptides with all zero

peptides <- peptides %>% 
  group_by(Sequence,fraction) %>% 
  filter(max(intensity)!=0)

peptides <- peptides[,c(1:5,28:30,6:27)]

#order table

peptides <- peptides %>% 
  arrange(ACC,fraction,Sequence)

## msa for matching fuzzy peptides
fuzzy <- peptides %>% 
  ungroup() %>% 
  distinct(Sequence) %>% 
  mutate(rev=stri_reverse(Sequence),
         no_tryp=ifelse(substr(rev,1,2)=='AK',sapply(strsplit(rev,'R'),'[',2),
                        ifelse(substr(rev,1,3)=='KAK',sapply(strsplit(rev,'R'),'[',3),
                               ifelse(substr(rev,1,1)!='K',sapply(strsplit(rev,'K'),'[',1),
                        sapply(strsplit(rev,'K'),'[',2)))),
         no_tryp=ifelse(substr(no_tryp,1,2)=='RR',sapply(strsplit(no_tryp,'R'),'[',3),
                        ifelse(substr(no_tryp,1,1)!='R',sapply(strsplit(no_tryp,'R'),'[',1),
                        sapply(strsplit(no_tryp,'R'),'[',2)))) %>% 
  arrange(no_tryp) %>% 
  mutate(match=substr(no_tryp,1,5),
         variable=ifelse(duplicated(match)==T,'yes','no')) 

variable_peptides <- filter(fuzzy, variable=='yes')

fuzzy <- filter(fuzzy, match  %in% variable_peptides$match)

match <- fuzzy %>% 
  select(Sequence,match)





fuzzy <- peptides %>% 
  left_join(match) %>% 
  ungroup() %>% 
  na.omit () %>% 
  distinct(ACC,genotype,Sequence,replicate,.keep_all = T) %>% 
  filter(nchar(match) > 1)



main_frame <- data.frame()
total=878
x=0
pb <- progress_bar$new(  format = "  completed [:bar] :percent eta: :eta",
                         total = 878, clear = FALSE, width= 60)

for(prot in unique(fuzzy$ACC)){
  for (seq in unique(filter(fuzzy,ACC==prot)$match)){
    tryCatch({
      x=x+1
      prot_frame=filter(fuzzy,ACC==prot)
      msa_frame=msa(prot_frame$Sequence,order='input',type='protein')
      msa_frame=msa_frame@unmasked %>% as.data.frame()
      prot_frame = bind_cols(prot_frame,msa_frame)
      main_frame = bind_rows(main_frame,prot_frame)
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    pb$tick()
    
  }
}
close(pb)

main_frame <- main_frame %>% 
  arrange(ACC,match)
#save main_frame
#write_tsv(main_frame,'data/main_semi_tryp_peptides.tsv')
#load main_frame
main_frame <- fread('data/main_semi_tryp_peptides.tsv')

## single plot for 75KDa


  p = ggplot(filter(main_frame,ACC=='AT5G37510'),
             aes(x,intensity/1000000,fill=genotype))+
    geom_boxplot(position = position_dodge(width=1))+
    #geom_point(data=filter(CI_peptides, Name == prot),
    #           aes(Sequence,intensity/1000000,fill=genotype),size=2,pch=21,position = position_dodge(width=1))+
    #facet_wrap(~fraction)+
    ggtitle(label=prot)+
    facet_wrap(~match,scales='free')+
    labs(y='Intensity [M counts]')
  ggsave(filename = paste0(filename,'.png'),path = 'images/petide_level',device = 'png',dpi=1080,plot = p)
  message(paste0(filename,' done!'))

## single plot for each protein (350)
  
pb <- progress_bar$new(  format = "  completed [:bar] :percent eta: :eta",
                         total = length(unique(main_frame$ACC)), clear = FALSE, width= 60)
  
for(prot in unique(main_frame$ACC)){
  desc = filter(main_frame, ACC==prot) %>% distinct(desc)
  ifelse(nrow(filter(main_frame,ACC==prot,fraction=='p'))==0,
         assign('p1',ggplot(mtcars, aes(wt, mpg))),
         assign('p1',ggplot(filter(main_frame,ACC==prot,fraction =='p'),
                            aes(x,intensity/1000000,fill=genotype))+
                  geom_boxplot(position = position_dodge(width=1))+
                  ggtitle(label=paste(prot,'- membrane','\n',desc[1,1]))+
                  facet_wrap(~match,scales='free')+
                  labs(y='Intensity [M counts]')+
                  theme(axis.text.x = element_text(size=8,angle = 15,vjust = .7))))
  ifelse(nrow(filter(main_frame,ACC==prot,fraction=='s'))==0,
         assign('p2',ggplot(mtcars, aes(wt, mpg))),
         assign('p2',ggplot(filter(main_frame,ACC==prot,fraction =='s'),
                            aes(x,intensity/1000000,fill=genotype))+
                  geom_boxplot(position = position_dodge(width=1))+
                  ggtitle(label=paste(prot,'- soluble','\n',desc[1,1]))+
                  facet_wrap(~match,scales='free')+
                  labs(y='Intensity [M counts]')+
                  theme(axis.text.x = element_text(size=8,angle = 15,vjust = .7))))
  p = plot_grid(p1, p2, labels = 'AUTO')
                                                                                     
  ggsave(filename = paste0(prot,'.png'),path = 'images/semi_tryp',device = 'png',dpi=1080,plot = p,height = 9,width = 16)
  message(paste0('\n',prot,' done!'))
  pb$tick()
  
}







