library(ggplot2)
library(tidyverse)
library(vegan)
setwd("/Users/yangjing/Desktop/lab/Dinoflagellate project/data/")

com = read.csv("level-3.csv",check.names = F)
com = com[,1:61] %>% column_to_rownames(., var = "index")
com %>% rownames()
rowSums(com) %>% min -> depth
com = rrarefy(com, depth)
com = com/depth

# sort most abundant groups -----------------------------------------------
com1 = com %>% colSums() %>% sort(decreasing = TRUE) %>% .[1:10] %>% data.frame()
com1 %>% rownames()
m <- match(rownames(com1), colnames(com))

# made dataframe ----------------------------------------------------------

df = 
  cbind(com[,m], data.frame(apply(com[,-m],1,sum)))%>% 
  rename(Others = `apply.com....m...1..sum.`) %>% t() %>% 
  data.frame() %>% 
  rownames_to_column(., var = "Dominant bacteria") %>% as_tibble() %>% 
  pivot_longer(!`Dominant bacteria`, names_to = "Sample", values_to = "Seqs") %>% 
  mutate(Species = Sample %>% str_split("_") %>% sapply('[', 1) %>% 
           str_replace("Acar", "Acartia") %>% str_replace("Para", "Paracalanus")) %>% 
  mutate(Day = Sample %>% str_split("_") %>% sapply('[', 3) %>% 
           str_replace("d0", "Day 0") %>% str_replace("d2", "Day 2")) %>% 
  mutate(Food = Sample %>% str_split("_") %>% sapply('[', 2) %>% 
           str_replace("at","At") %>% str_replace("aa", "Aa") %>% str_replace("sw", "Seawater")) 


# relative abundance plot -------------------------------------------------

df$Sample <- factor(df$Sample, levels = unique(df$Sample),ordered = T)
df$`Dominant bacteria` <- factor(df$`Dominant bacteria`, levels = unique(df$`Dominant bacteria`))
g = ggplot(df, aes(x = Sample, y = Seqs, fill = `Dominant bacteria`))+
  geom_bar(stat="identity", width = 1)+
  theme_classic()+
  theme(axis.text.x = element_text(color="black", size=12,angle = 90, hjust = 1, vjust = .5),
        axis.text.y = element_text(color="black", size=15),
        title = element_text(color="black", size=20),
        legend.title=element_text(size=15), 
        legend.text=element_text(size=12),
        legend.position="right",
        axis.title=element_text(size=20))+
  ylab("Relative abundance")+
  ggtitle("Gut microbiota community (16S rDNA)")+
  scale_fill_brewer(palette="Set3")

g
