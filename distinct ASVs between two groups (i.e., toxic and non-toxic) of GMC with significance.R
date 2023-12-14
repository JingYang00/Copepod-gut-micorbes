setwd("/Users/yangjing/Desktop/lab/Dinoflagellate project/data/data copy/")
df1 <- read.delim("cop-16s-fitered.txt",sep="\t",check.names = FALSE, row.names = 1)
library(tidyverse)
library(scales)
library(ggpubr)
library(vegan)

# ### distinct ASVs between treatment ---------------------------------------

tax = read.delim("taxonomy.tsv") %>% as_tibble()
#tax[(str_detect(tax$Taxon,"Chloroplast")|str_detect(tax$Taxon,"Mitochondria")),] %>% as_tibble() %>% pull(Taxon) %>% length()
tax1 = tax[!(str_detect(tax$Taxon,"Chloroplast")|str_detect(tax$Taxon,"Mitochondria")),]

asv_table = df1[match(tax$Feature.ID,rownames(df1)),] %>% t()
asv_table1 = df1[match(tax1$Feature.ID,rownames(df1)),] %>% t()


# rarefy to same depth ----------------------------------------------------

asv_table = asv_table[which(rowSums(asv_table)>1000),]
asv_table = rrarefy(asv_table, min(rowSums(asv_table)))
asv_table = asv_table[,which(colSums(asv_table)>0)]
asv_table = asv_table[!str_detect(rownames(asv_table),"d7"),]


asv_table1 %>% rowSums() %>% min()
asv_table1 = asv_table1[which(rowSums(asv_table1)>1000),]
asv_table1 = rrarefy(asv_table1, min(rowSums(asv_table1)))
asv_table1 %>% rowSums()
## remove non-exist asv
asv_table1 = asv_table1[,which(colSums(asv_table1)>0)]
asv_table1 = asv_table1[!str_detect(rownames(asv_table1),"d7"),]



## Acartia. ASV level

asv_acar = asv_table[str_detect(rownames(asv_table), "Acar"),]
asv_acar = rrarefy(asv_acar, min(rowSums(asv_acar)))
asv_acar_ra = asv_acar/min(rowSums(asv_acar))

set.seed(996)
simper(asv_acar_ra, asv_acar %>% rownames %>% str_split("_") %>% sapply('[',2))->acar_asv_table

acar_aa_at = acar_asv_table$aa_at %>% as_tibble()
acar_aa_at$ava[acar_aa_at$ava == 0] <- 1
acar_aa_at$ava[acar_aa_at$ava == 1] <- min(acar_aa_at$ava)
acar_aa_at_asv = acar_aa_at %>% 
  mutate(taxon = tax$Taxon[match(acar_asv_table$aa_at$species, tax1$Feature.ID)]) %>% 
  mutate(Class = taxon %>% str_split(";") %>% sapply('[',3)) %>%
  mutate(Phylum = taxon %>% str_split(";") %>% sapply('[',2)) %>% 
  mutate(all = paste(Phylum, Class, sep = ";")) %>% 
  mutate(a = log(avb/ava,2)) %>% 
  filter(avb > 0.001) %>% 
  mutate(Bacteria = if_else(p > 0.05 | a >= -1 & a <= 1, "No change", if_else(is.na(all), "Unassigned", all))) %>%
  ggplot(aes(x = log(avb/ava,2) , y = -log(p,10)))+
  geom_point(aes(color = Bacteria), size = 4) +
  scale_color_manual(values = c("brown2","deepskyblue","gold","purple","lightgray"))+
  geom_hline(yintercept = -log(0.05, 10), color = "black", linetype = "dashed")+
  geom_vline(xintercept = 1, color = "black", linetype = "dashed")+
  geom_vline(xintercept = -1, color = "black", linetype = "dashed")+
  labs(x = expression(Log[2]~(Fold~Change)), 
       y = "- lg p-value")+
  theme_bw()+
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"),
        axis.text = element_text(color="black", size=15),
        title = element_text(color="black", size=15),
        legend.title=element_text(size=15), 
        legend.text=element_text(size=12),
        axis.title=element_text(size=15))+
  scale_y_continuous(limits = c(0,2), breaks=c(0, 1, 2))
acar_aa_at_asv


acar_aa_at_up = acar_aa_at %>% 
  mutate(taxon = tax$Taxon[match(acar_asv_table$aa_at$species, tax1$Feature.ID)]) %>% 
  mutate(Class = taxon %>% str_split(";") %>% sapply('[',3)) %>% 
  filter(avb > 0.001) %>% 
  filter(p <= 0.05) %>% 
  mutate(a = log(avb/ava,2)) %>% 
  filter(a > 1) %>% 
  arrange(desc(ava))

acar_aa_at_down = acar_aa_at %>% 
  mutate(taxon = tax$Taxon[match(acar_asv_table$aa_at$species, tax1$Feature.ID)]) %>% 
  mutate(Class = taxon %>% str_split(";") %>% sapply('[',3)) %>% 
  filter(avb > 0.001) %>% 
  filter(p <= 0.05) %>% 
  mutate(a = log(avb/ava,2)) %>% 
  filter(a < -1) %>% 
  arrange(desc(ava))

write_csv(acar_aa_at_up,"acar_aa_at_up.csv")
write_csv(acar_aa_at_down,"acar_aa_at_down.csv")




# Paracalanus. ASV level

asv_para = asv_table1[str_detect(rownames(asv_table1), "Para"),]
asv_para = rrarefy(asv_para, min(rowSums(asv_para)))
asv_para_ra = asv_para/min(rowSums(asv_para))

set.seed(456)
simper(asv_para_ra, asv_para %>% rownames %>% str_split("_") %>% sapply('[',2)) ->para_asv_table

para_aa_at = para_asv_table$aa_at %>% as_tibble()
para_aa_at$ava[para_aa_at$ava == 0] <- 1
para_aa_at$ava[para_aa_at$ava == 1] <- min(para_aa_at$ava)
para_aa_at_asv = para_aa_at %>% 
  mutate(taxon = tax$Taxon[match(para_asv_table$aa_at$species, tax$Feature.ID)]) %>% 
  mutate(Class = taxon %>% str_split(";") %>% sapply('[',3)) %>% 
  mutate(Phylum = taxon %>% str_split(";") %>% sapply('[',2)) %>% 
  mutate(all = paste(Phylum, Class, sep = ";")) %>% 
  mutate(a = log(avb/ava,2)) %>% 
  filter(avb > 0.001) %>% 
  mutate(Bacteria = if_else(p > 0.05 | a >= -1 & a <= 1, "No change", if_else(is.na(all), "Unassigned", all))) %>%
  ggplot(aes(x = log(avb/ava,2) , y = -log(p,10)))+
  geom_point(aes(color = Bacteria), size = 4) +
  scale_color_manual(values = c("brown2", "deepskyblue", "lightpink",
                                "yellowgreen","violet","yellow","lightgray"))+
  geom_hline(yintercept = -log(0.05, 10), color = "black", linetype = "dashed")+
  geom_vline(xintercept = 1, color = "black", linetype = "dashed")+
  geom_vline(xintercept = -1, color = "black", linetype = "dashed")+
  labs(x = expression(Log[2]~(Fold~Change)), 
       y = "- lg p-value")+
  theme_bw()+
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"),
        axis.text = element_text(color="black", size=15),
        title = element_text(color="black", size=15),
        legend.title=element_text(size=15), 
        legend.text=element_text(size=12),
        axis.title=element_text(size=15))
para_aa_at_asv


para_aa_at_up = para_aa_at %>% 
  mutate(taxon = tax$Taxon[match(para_asv_table$aa_at$species, tax$Feature.ID)]) %>% 
  mutate(Class = taxon %>% str_split(";") %>% sapply('[',3)) %>% 
  filter(avb > 0.001) %>% 
  filter(p <= 0.05) %>% 
  mutate(a = log(avb/ava,2)) %>% 
  filter(a > 1) %>% 
  arrange(desc(ava))

para_aa_at_down = para_aa_at %>% 
  mutate(taxon = tax1$Taxon[match(para_asv_table$aa_at$species, tax1$Feature.ID)]) %>% 
  mutate(Class = taxon %>% str_split(";") %>% sapply('[',3)) %>%
  filter(avb > 0.001) %>% 
  filter(p <= 0.05) %>% 
  mutate(a = log(avb/ava,2)) %>% 
  filter(a < -1) %>% 
  arrange(desc(ava))

write_csv(para_aa_at_up,"para_aa_at_up.csv")
write_csv(para_aa_at_down,"para_aa_at_down.csv")


ggarrange(acar_aa_at_asv, NULL, para_aa_at_asv, nrow = 3, heights = c(1,0.05,1),legend = "right", common.legend = F)


