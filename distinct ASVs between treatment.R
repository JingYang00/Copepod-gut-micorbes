setwd("/Users/yangjing/Desktop/lab/Dinoflagellate project/data/data copy/")
df1 <- read.delim("cop-16s-fitered.txt",sep="\t",check.names = FALSE, row.names = 1)
library(tidyverse)
library(scales)
library(ggpubr)

asv_table1 = df1 %>% t()
asv_table1 %>% dim()

# rarefy to same depth ----------------------------------------------------

asv_table1 %>% rowSums() %>% min()
asv_table1 = asv_table1[which(rowSums(asv_table1)>1000),]
asv_table1 %>% dim()  


library(vegan)
asv_table1 = rrarefy(asv_table1, min(rowSums(asv_table1)))
asv_table1 %>% rowSums()
asv_table1 %>%  dim
## remove non-exist asv
asv_table1 = asv_table1[,which(colSums(asv_table1)>0)]
asv_table1 %>% dim



# ### distinct ASVs between treatment ---------------------------------------

tax = read.delim("taxonomy.tsv")
tax %>% as_tibble()

## Acartia. ASV level
asv_acar = asv_table1[str_detect(rownames(asv_table1), "Acar"),]
asv_acar = asv_acar[,which(colSums(asv_acar)>0)]
asv_acar = rrarefy(asv_acar, min(rowSums(asv_acar)))
asv_acar_ra = asv_acar/min(rowSums(asv_acar))

simper(asv_acar_ra, asv_acar %>% rownames %>% str_split("_") %>% sapply('[',2)) -> acar_asv_table

#find the bacteria

acar_aa_at_up = acar_asv_table$aa_at %>% as_tibble %>% 
  mutate(taxon = tax$Taxon[match(acar_asv_table$aa_at$species, tax$Feature.ID)]) %>% 
  mutate(Class = taxon %>% str_split(";") %>% sapply('[',3)) %>% 
  filter(p < 0.05) %>% 
  filter(ava < avb) %>%
  arrange(desc(ava))

acar_aa_at_down = acar_asv_table$aa_at %>% as_tibble %>% 
  mutate(taxon = tax$Taxon[match(acar_asv_table$aa_at$species, tax$Feature.ID)]) %>% 
  mutate(Class = taxon %>% str_split(";") %>% sapply('[',3)) %>% 
  filter(p < 0.05) %>% 
  filter(ava > avb) %>% 
  arrange(desc(ava))

write_csv(acar_aa_at_up,"acar_aa_at_up.csv")
write_csv(acar_aa_at_down,"acar_aa_at_down.csv")


#figure

acar_aa_at = acar_asv_table$aa_at %>% as_tibble()
acar_aa_at$avb[acar_aa_at$avb == 0] <- 0.0001
acar_aa_at_asv = acar_aa_at %>% 
  mutate(taxon = tax$Taxon[match(acar_asv_table$aa_at$species, tax$Feature.ID)]) %>% 
  mutate(Class = taxon %>% str_split(";") %>% sapply('[',3)) %>% 
  ggplot(aes(x = log(avb/ava,2) , y = -log(p,10)))+
  geom_point(aes(color = Class), size = 4)+
  geom_hline(yintercept = -log(0.05, 10), color = "black", linetype = "dashed")+
  geom_vline(xintercept = 1, color = "black", linetype = "dashed")+
  geom_vline(xintercept = -1, color = "black", linetype = "dashed")+
  labs(x = "log(%Toxic / %Non-toxic, 2)", 
       y = "- lg p-value", 
       title ="Acartia: Toxic / Non-toxic")+
  theme_bw()+
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))+
  theme(axis.text = element_text(color="black", size=15),
        title = element_text(color="black", size=15))+
  theme(legend.title=element_text(size=13,face = "bold"), 
        legend.text=element_text(size=8,face = "bold"),
        axis.title=element_text(size=15))
acar_aa_at_asv


# para
asv_para = asv_table1[str_detect(rownames(asv_table1), "Para"),]
asv_para = asv_para[,which(colSums(asv_para)>0)]
asv_para = rrarefy(asv_para, min(rowSums(asv_para)))
asv_para_ra = asv_para/min(rowSums(asv_para))

simper(asv_para_ra, asv_para %>% rownames %>% str_split("_") %>% sapply('[',2)) ->para_asv_table

#find the bacteria

para_aa_at_up = para_asv_table$aa_at %>% as_tibble %>% 
  mutate(taxon = tax$Taxon[match(para_asv_table$aa_at$species, tax$Feature.ID)]) %>% 
  mutate(Class = taxon %>% str_split(";") %>% sapply('[',3)) %>% 
  filter(p < 0.05) %>% 
  filter(ava < avb) %>%
  arrange(desc(ava))

para_aa_at_down = para_asv_table$aa_at %>% as_tibble %>% 
  mutate(taxon = tax$Taxon[match(para_asv_table$aa_at$species, tax$Feature.ID)]) %>% 
  mutate(Class = taxon %>% str_split(";") %>% sapply('[',3)) %>% 
  filter(p < 0.05) %>% 
  filter(ava > avb) %>% 
  arrange(desc(ava))

write_csv(para_aa_at_up,"para_aa_at_up.csv")
write_csv(para_aa_at_down,"para_aa_at_down.csv")

#figure

para_aa_at = para_asv_table$aa_at %>% as_tibble()
para_aa_at$avb[para_aa_at$avb == 0] <- 0.0001
para_aa_at_asv = para_aa_at %>% 
  mutate(taxon = tax$Taxon[match(para_asv_table$aa_at$species, tax$Feature.ID)]) %>% 
  mutate(Class = taxon %>% str_split(";") %>% sapply('[',3)) %>% 
  ggplot(aes(x = log(avb/ava,2) , y = -log(p,10)))+
  geom_point(aes(color = Class), size = 4)+
  geom_hline(yintercept = -log(0.05, 10), color = "black", linetype = "dashed")+
  geom_vline(xintercept = 1, color = "black", linetype = "dashed")+
  geom_vline(xintercept = -1, color = "black", linetype = "dashed")+
  labs(x = "log(%Toxic / %Non-toxic, 2)", 
       y = "- lg p-value", 
       title ="Paracalanus: Toxic / Non-toxic")+
  theme_bw()+
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))+
  theme(axis.text = element_text(color="black", size=15),
        title = element_text(color="black", size=15))+
  theme(legend.title=element_text(size=13,face = "bold"), 
        legend.text=element_text(size=8,face = "bold"),
        axis.title=element_text(size=15))
para_aa_at_asv

ggarrange(acar_aa_at_asv, NULL, para_aa_at_asv, ncol = 3, widths = c(1,0.05,1),legend = "bottom", common.legend = T)
