setwd("/Users/yangjing/Desktop/lab/Dinoflagellate project/data/data copy/")
df <- read.delim("cop-16s-fitered.txt",sep="\t",check.names = FALSE, row.names = 1)
library(tidyverse)
library(scales)
library(ggpubr)

tax = read.delim("taxonomy.tsv") %>% as_tibble()
tax1 = tax[!(str_detect(tax$Taxon,"Chloroplast")|str_detect(tax$Taxon,"Mitochondria")),]

asv_table = df[match(tax1$Feature.ID,rownames(df)),] %>% t()
asv_table = asv_table[which(rowSums(asv_table)>1000),]
asv_table = asv_table[,which(colSums(asv_table)>0)]
asv_table = asv_table[!str_detect(rownames(asv_table),"d7"),]



# alpha diversity dataframe---------------------------------------------------------

df_alpha = asv_table %>% alpha() %>% 
  rownames_to_column(., var = "Sample") %>% 
  mutate(Source = Sample %>% str_split("_") %>% sapply('[', 2) %>%
           str_replace_all(c("aa"="Non-toxic algae","at"="Toxic algae")) %>%
           str_replace("sw", "Field"))

alpha_acar = df_alpha[!str_detect(df_alpha$Sample,"Para"),]
alpha_para = df_alpha[!str_detect(df_alpha$Sample,"Acar"),]

my_comparisons = list( c("Field", "Non-toxic algae"), c("Field", "Toxic algae"), c("Non-toxic algae", "Toxic algae"))

# shannon -----------------------------------------------------------------

qshannona=
  alpha_acar %>%
  ggplot(aes(x = Source, y = Shannon)) +
  geom_boxplot(aes(color = Source))+
  scale_color_manual(values = c("grey","#009900","#FF0000"))+
  labs(title = '(a)')+
  theme_classic()+
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"),
        axis.text = element_text(color="black", size=17),
        axis.title.x = element_blank(),
        title = element_text(color="black", size=17),
        legend.position = "none")+
  scale_y_continuous(limits = c(0.5,6),labels = label_number(accuracy = 0.01))+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", size=7) 
qshannona

qshannonp=
  alpha_para %>%
  ggplot(aes(x = Source, y = Shannon)) +
  geom_boxplot(aes(color = Source))+
  scale_color_manual(values = c("grey","#009900","#FF0000"))+
  labs(title = '(b)')+
  theme_classic()+
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"),
        axis.text = element_text(color="black", size=17),
        axis.title.x = element_blank(),
        title = element_text(color="black", size=17),
        legend.position = "none")+
  scale_y_continuous(limits = c(0.5,6), labels = label_number(accuracy = 0.01))+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", size=7) 
qshannonp



# richness ---------------------------------------------------------------

qrichnessa=
  alpha_acar %>%
  ggplot(aes(x = Source, y = Richness)) +
  geom_boxplot(aes(color = Source))+
  scale_color_manual(values = c("grey","#009900","#FF0000"))+
  labs(title = '(c)')+
  theme_classic()+
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"),
        axis.text = element_text(color="black", size=17),
        axis.title.x = element_blank(),
        title = element_text(color="black", size=17),
        legend.position = "none")+
  scale_y_continuous(limits = c(0,400), labels = label_number(accuracy = 1))+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", size=7) 
qrichnessa

qrichnessp=
  alpha_para %>%
  ggplot(aes(x = Source, y = Richness)) +
  geom_boxplot(aes(color = Source))+
  scale_color_manual(values = c("grey","#009900","#FF0000"))+
  labs(title = '(d)')+
  theme_classic()+
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"),
        axis.text = element_text(color="black", size=17),
        axis.title.x = element_blank(),
        title = element_text(color="black", size=17),
        legend.position = "none")+
  scale_y_continuous(limits = c(0,350), labels = label_number(accuracy = 1))+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", size=7) 
qrichnessp

ggarrange(qshannona, qshannonp, qrichnessa, qrichnessp, ncol = 2, nrow = 2)

  
# beta diversity dataframe and anosim test----------------------------------------------------------

##Copepod
grp_cop <- data.frame(Sample = rownames(asv_table)) %>% 
  mutate(Copepod = Sample %>% str_split("_") %>% sapply('[',1))
dd = vegdist(asv_table)
anosim(dd,grp_cop$Copepod) # R: 0.8234, Significance: 0.001


##Acartia
df_acar = asv_table[str_detect(rownames(asv_table),"Acar"),]
grp_acar <- data.frame(Sample = rownames(df_acar)) %>% 
  mutate(Alage= Sample %>% str_split("_") %>% sapply('[',2))

d1 = vegdist(df_acar)
anosim(d1, grp_acar$Alage) # R: 0.3626, Significance: 0.001


##Paracalanus
df_para = asv_table[str_detect(rownames(asv_table),"Para"),]
grp_para <- data.frame(Sample = rownames(df_para)) %>% 
  mutate(Alage= Sample %>% str_split("_") %>% sapply('[',2))

d2 = vegdist(df_para) 
anosim(d2, grp_para$Alage) # R: 0.3456, Significance: 0.001



# beta diversity of Acartia----------------------------------------------------------

grp_acar1 = grp_acar %>%  mutate(Source = Sample %>% str_split("_") %>% sapply('[', 2) %>%
                                   str_replace_all(c("aa"="Non-toxic algae","at"="Toxic algae")) %>%
                                   str_replace("sw", "Field"))
nmds = monoMDS(d1)
dfa = data.frame(MDS1 = nmds$points[,1],
                 MDS2 = nmds$points[,2]) %>% as.data.frame()
dfa = cbind(dfa, grp_acar1)

pa <- dfa %>% 
  ggplot(aes(x = MDS1, y = MDS2))+
  geom_point(aes(color = Source), size = 8)+
  scale_color_manual(values = c("grey","#009900","#FF0000"))+
  labs(title = '(b)',
       x = 'NMDS1',
       y = 'NMDS2')+
  theme_classic()+
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"),
        title = element_text(color="black", size=15),
        legend.title=element_blank(), 
        legend.text=element_text(size=17),
        axis.title=element_text(size=17),
        axis.text = element_text(color="black", size=17))+
  stat_ellipse(aes(color = Source),level=0.6)
pa

# beta diversity of Paracalanus -------------------------------------------

grp_para1 = grp_para %>%  mutate(Source = Sample %>% str_split("_") %>% sapply('[', 2) %>%
                                   str_replace_all(c("aa"="Non-toxic algae","at"="Toxic algae")) %>%
                                   str_replace("sw", "Field"))

nmds = monoMDS(d2)
dfp = data.frame(MDS1 = nmds$points[,1],
                 MDS2 = nmds$points[,2]) %>% as.data.frame() 
dfp = cbind(dfp, grp_para1)

pp <- dfp %>% 
  ggplot(aes(x = MDS1, y = MDS2))+
  geom_point(aes(color = Source), size = 8)+
  labs(title = '(c)',
       x = 'NMDS1',
       y = 'NMDS2')+
  scale_color_manual(values = c("grey","#009900","#FF0000"))+
  theme_classic()+
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"),
        title = element_text(color="black", size=15),
        legend.title=element_blank(), 
        legend.text=element_text(size=17),
        axis.title=element_text(size=17),
        axis.text = element_text(color="black", size=17))+
  stat_ellipse(aes(color = Source),level=0.6)
pp

ggarrange(pa,pp, common.legend = TRUE, legend = "right")


# beta-comparison ---------------------------------------------------------

Para_aa = asv_table[str_detect(rownames(asv_table), "Para_aa"),]
Para_at = asv_table[str_detect(rownames(asv_table), "Para_at"),]
Para_sw = asv_table[str_detect(rownames(asv_table), "Para_sw"),]

Acar_aa = asv_table[str_detect(rownames(asv_table), "Acar_aa"),]
Acar_at = asv_table[str_detect(rownames(asv_table), "Acar_at"),]
Acar_sw = asv_table[str_detect(rownames(asv_table), "Acar_sw"),]

Para_aa = Para_aa[,which(colSums(Para_aa)>0)]
Para_at = Para_at[,which(colSums(Para_at)>0)]
Para_sw = Para_sw[,which(colSums(Para_sw)>0)]

Acar_aa = Acar_aa[,which(colSums(Acar_aa)>0)]
Acar_at = Acar_at[,which(colSums(Acar_at)>0)]
Acar_sw = Acar_sw[,which(colSums(Acar_sw)>0)]

# d1 = as.matrix(d1)
# d1 = d1[lower.tri(d1)]

d1 = vegdist(Acar_sw)
d2 = vegdist(Acar_aa)
d3 = vegdist(Acar_at)
d4 = vegdist(Para_sw)
d5 = vegdist(Para_aa)
d6 = vegdist(Para_at)

df = data.frame(Dissimilarity = c(d1,d2,d3,d4,d5,d6), 
                Treatment = c(rep("Acar_sw", length(d1)),
                              rep("Acar_aa", length(d2)),
                              rep("Acar_at", length(d3)),
                              rep("Para_sw", length(d4)),
                              rep("Para_aa", length(d5)),
                              rep("Para_at", length(d6))))
  
df %>% 
  mutate(Algae = Treatment %>% str_split("_") %>% sapply('[',2) %>%
           str_replace_all(c('aa'="Non-toxic algae",'at'="Toxic algae")) %>%
           str_replace('sw',"Field")) %>% 
  mutate(Copepod = Treatment %>% str_split("_") %>% sapply('[',1) %>% 
           str_replace_all(c('Acar'="Acartia",'Para'="Paracalanus"))) %>% 
  ggplot(aes(x = Copepod, y = Dissimilarity)) +
  geom_boxplot(aes(colour = Algae))+
  scale_color_manual(values = c("grey","#009900","#FF0000"))+
  theme_classic()+
  theme(panel.border = element_rect(fill = NA,color = "black", size = 2),
        axis.text = element_text(color = "black", size = 17),
        axis.text.x = element_text(face = "italic"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = "black",size = 17),
        legend.title = element_blank(),
        legend.text=element_text(size=14))

