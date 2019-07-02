library(data.table)
library(tidyverse)

peptides <- fread('data/peptides.txt')

#complex I annotation
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


### fix accession

peptides <- peptides %>% 
  filter(str_detect(`Leading razor protein`, 'AT')==T) %>% 
  mutate(ACC=substr(`Leading razor protein`,1,9)) 

#remove unnecessary columns
peptides <- peptides[,which(str_detect(colnames(peptides),'Count')==F)]
peptides <- peptides[,which(str_detect(colnames(peptides),'Experiment')==F)]
peptides <- peptides[,which(str_detect(colnames(peptides),'Intensity')==F)]

#filter for complexI

CI_peptides <- filter(peptides,ACC %in% cI_anno$ACC)
#add description
CI_peptides <- CI_peptides %>% 
  left_join(cI_anno)

#gather

CI_peptides <- CI_peptides %>% 
  gather(sample,intensity,`LFQ intensity clp1p_1`:`LFQ intensity wts_4`)

#remove more columns

CI_peptides <- CI_peptides %>% 
  dplyr::select(-c(Mass,Proteins,`Leading razor protein`,`Mod. peptide IDs`,`Fraction Average`,`Fraction Std. Dev.`,
            Reverse,`Potential contaminant`,`Evidence IDs`,`MS/MS IDs`,`LFQ intensity blank`))

#reorder for visibility

CI_peptides <- CI_peptides[,c(24:27,1:23)]

#read sample meta data from sample name

CI_peptides <- CI_peptides %>% 
  mutate(sample=sapply(strsplit(sample,' '),'[',3)) %>% 
  filter(str_detect(sample,'pool')==F) %>% 
  mutate(replicate=sapply(strsplit(sample,'_'),'[',2),
         genotype=sapply(strsplit(sample,'._'),'[',1),
         fraction=ifelse(str_detect(sample,'p_')==T,'p','s')) 


#remove peptides with all zero

CI_peptides <- CI_peptides %>% 
  group_by(Sequence,fraction) %>% 
  filter(max(intensity)!=0)

CI_peptides <- CI_peptides[,c(1:5,28:30,6:27)]

# CI_peptides <- CI_peptides %>%
#   group_by(Sequence,genotype,fraction) %>% 
#   mutate(max=max(intensity)) %>% 
#   dplyr::dplyr::select(c(1:8,ncol(.))) %>% 
#   ungroup() %>% 
#   filter(!(intensity ==0 & max !=0)) 



#median up replicates

CI_peptides_median <- CI_peptides %>% 
  group_by(Sequence,genotype,Name,fraction) %>% 
  summarize(median=median(intensity))

#replace backslash with underscore for names
CI_peptides_median <- CI_peptides_median %>% 
  ungroup() %>% 
  mutate(Name=gsub('/','_',Name))


#stacked bars

ggplot(CI_peptides_median,aes(genotype,median,fill=Sequence))+
  geom_bar(stat='identity')+
  facet_wrap(~Name,scales='free')+
  theme(legend.position = 'none')

#single plots for each protein, with bar for each peptide

for (prot in unique(CI_peptides_median$Name)){
  filename=prot
  p = ggplot(filter(CI_peptides_median,Name==prot),
             aes(Sequence,median/1000000,fill=genotype))+
    geom_bar(stat='identity', position = position_dodge(width=1))+
    geom_point(data=filter(CI_peptides, Name == prot),
               aes(Sequence,intensity/1000000,fill=genotype),size=2,pch=21,position = position_dodge(width=1))+
    facet_wrap(~fraction)+
    theme(legend.position = 'none')+
    ggtitle(label=prot)+
    labs(y='Intensity [M counts]')
  ggsave(filename = paste0(filename,'.png'),path = 'images/petide_level',device = 'png',dpi=1080,plot = p)
  message(paste0(filename,' done!'))
}







