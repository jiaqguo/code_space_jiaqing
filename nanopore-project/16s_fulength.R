setwd("D:/OneDrive - Oklahoma A and M System/tranfer documents/16s_fulength/cluster97")
setwd("C:/Users/jiaguo/OneDrive - Oklahoma A and M System/tranfer documents/16s_fulength/cluster99")

# down stream analysis

# modify tax table
tax <- read.csv("taxonomy_kmer.csv",header = T,row.names = 1)

for (i in 1:nrow(tax)){
  if (tax$Domain[i] == "Unassigned"){
    tax$Phylum[i] = tax$Class[i] = tax$Order[i] =  tax$Family[i] =  tax$Genus[i] =  tax$Species[i] = "Unassigned"
  } else if (tax$Phylum[i] == ""){
    tax$Phylum[i] = tax$Class[i] = tax$Order[i] =  tax$Family[i] =  tax$Genus[i] =  tax$Species[i] = paste(tax$Domain[i],"Unassigned",sep = "_")
  } else if (tax$Class[i] == ""){
    tax$Class[i] = tax$Order[i] =  tax$Family[i] =  tax$Genus[i] =  tax$Species[i] = paste(tax$Phylum[i],"Unassigned",sep = "_")
  } else if (tax$Order[i] == ""){
    tax$Order[i] =  tax$Family[i] =  tax$Genus[i] =  tax$Species[i] = paste(tax$Class[i],"Unassigned",sep = "_")
  } else if (tax$Family[i] == ""){
    tax$Family[i] =  tax$Genus[i] =  tax$Species[i] = paste(tax$Order[i],"Unassigned",sep = "_")
  } else if (tax$Genus[i] == ""){
    tax$Genus[i] =  tax$Species[i] = paste(tax$Family[i],"Unassigned",sep = "_")
  } else if (tax$Species[i] == ""){
    tax$Species[i] = paste(tax$Genus[i],"Unassigned",sep = "_")
  }
}

write.csv(tax,"taxonomy_kmer.csv")

# a diversity------------------
library(tidyverse)
library(vegan)
library(rstatix)
library(ggpubr)
library(car) # equal variance

feature_table <- read.csv("otutab.csv", header = T,  row.names = 1) %>% 
  rownames_to_column()

tax <- read.csv("taxonomy_kmer.csv",header = T,row.names = 1) %>% 
  rownames_to_column()

sample_data <- read.csv("metadata.csv",header = T,row.names = 1) %>% 
  rownames_to_column()

otu_tab <- left_join(feature_table,tax,by="rowname") %>% 
  select(2:89,ncol(.)) %>% 
  group_by(Species) %>% 
  summarise(across(1:ncol(.)-1, sum)) %>% 
  column_to_rownames(.,var = "Species") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  left_join(sample_data,.,by="rowname") %>% 
  filter(!is.na(Sample))

# cal a-diversity
group_info <- otu_tab$seg
id_info <- otu_tab$Sample

df2 <- otu_tab %>% 
  select(-c(1:3))

shan_index <- data.frame(id=id_info,group=group_info,shannon=diversity(df2, index = "shannon")) 

obs_spec <- data.frame(id=id_info,group=group_info,obs=specnumber(df2))

pielou <- data.frame(id=id_info,group=group_info,pielou=shan_index$shannon / log(obs_spec$obs)) 

a_diversity <- left_join(shan_index,obs_spec[,c(1,3)],by = "id") %>% 
  left_join(.,pielou[,c(1,3)],by="id") %>% 
  filter(!group %in% c("Neg","Pos"))

a_diversity$group <- factor(a_diversity$group,levels = c("Crop","Gizzard","Duodenum","Jejunum","Ileum","Cecum",    "Colon"))


# Normality test and Homogeneity of variance test
check_function <- function(formula, data){
  fit <- aov(formula, data = data)
  res <- residuals(fit)
  par(mfrow=c(1,2))
  qqnorm(res); qqline(res)
  cat("\nNormality test:\n")
  print(shapiro.test(res))
  cat("\nHomogeneity of variance test:\n")
  plot(a_diversity$group,res)
  abline(h=0)
  # using car package
  print(leveneTest(a_diversity$shannon ~ a_diversity$group))
  par(mfrow=c(1,1))
}

# statistical analysis
kw_function <- function(formula, data){
  kw_test <- kruskal.test(formula,data = data)
  print(kw_test)
  posthoc <- dunn_test(data = data,formula,p.adjust.method = "fdr")
  posthoc$compar <- paste(posthoc$group1,"-",posthoc$group2,sep="")
  print(posthoc)
  return(posthoc)
}

aov_function <- function(formula, data){
  aov <- aov(formula,data = data)
  print(summary(aov))
  posthoc <- as.data.frame(TukeyHSD(aov)[[1]])
  posthoc$compar <- rownames(posthoc)
  colnames(posthoc)[colnames(posthoc) == "p adj"] <- "p.adj"
  posthoc<-separate(data = posthoc,col = compar,into = c("group1", "group2"),sep = "-")
  posthoc$compar <- paste(posthoc$group2,"-",posthoc$group1,sep = "")
  print(posthoc)
  return(posthoc)
}

compar <- pielou ~ group

check_function(compar, a_diversity)

posthoc <- kw_function(compar, data=a_diversity)

posthoc <- aov_function(compar, data=a_diversity)


# get the compact letter display for lists of comparisons
library(rcompanion)
max <- a_diversity %>% 
  group_by(group) %>% 
  summarise(max=max(pielou)+max(pielou)*0.01)
names(max)[1] <- "Group"
sig_mark <- cldList(p.adj ~ compar, data = posthoc) %>% 
  left_join(.,max, by = "Group")
sig_mark

# plot
library(paletteer)
library(ggnewscale)
col1 <- paletteer_d("ggsci::category20_d3")
light_col <- col1[c(14,12,19,13,20,18,15)]
dark_col <- col1[c(4,2,9,3,1,8,5)]

# Shannon Index
# Observed Species
# Pielou's Eveness
p1 <- ggplot(a_diversity,aes(x = group,y = pielou)) +
  geom_boxplot(aes(colour = group,fill = group),outliers = F) +
  scale_color_manual(values = dark_col) +
  scale_fill_manual(values = light_col) +
  #geom_text(aes(label = id),position = position_jitter(width = 0.1), size = 3, vjust = -0.5) +
  labs(x="",y="Pielou's Eveness") +
  scale_x_discrete(expand = c(0.1, 0)) +  # default 0.2
  scale_y_continuous(expand = c(0.1, 0)) +  # default 0.05
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 8,color = "black",angle = 45,hjust = 1),
    axis.text.y = element_text(size = 8,color = "black"),
    axis.title = element_text(size = 9,color = "black"),
    legend.text = element_text(size =9, color = "black"), 
    legend.title = element_blank(),
    legend.position = "none"
  )+
  geom_text(aes(label = Letter, x=Group, y =(max)), data=sig_mark,vjust = -0.5,size=3)+
  new_scale_fill() +
  geom_jitter(aes(colour = group,fill = group),shape = 21, width = 0.1,size=0.5) +
  scale_fill_manual(values = dark_col)

p1
library(ggh4x)
p1 <- p1 + force_panelsizes(total_width=unit(1.7, "in"),total_height=unit(1.7, "in"))
p1 
ggsave("pielou.pdf",p1)


# b diversity------------------
library(tidyverse)
library(vegan)
library(ape)

feature_table <- read.csv("otutab.csv", header = T,  row.names = 1) %>% 
  rownames_to_column()

tax <- read.csv("taxonomy_kmer.csv",header = T,row.names = 1) %>% 
  rownames_to_column()

sample_data <- read.csv("metadata.csv",header = T,row.names = 1) %>% 
  rownames_to_column()

otu_tab <- left_join(feature_table,tax,by="rowname") %>% 
  select(2:89,ncol(.)) %>% 
  group_by(Species) %>% 
  summarise(across(1:ncol(.)-1, sum)) %>% 
  column_to_rownames(.,var = "Species") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  left_join(sample_data,.,by="rowname") %>% 
  filter(!is.na(Sample))%>% 
  filter(!seg %in% c("Neg","Pos"))

group_info <- otu_tab$seg
id_info <- otu_tab$Sample

df2 <- otu_tab %>% 
  select(2,4:ncol(.)) %>% 
  column_to_rownames(var = "Sample")

# pcoa analysis
dis_bray <- vegdist(df2, method = "bray")
dis_jaccard <- vegdist(df2, method = "jaccard")

pcoa_list <- pcoa(dis_bray)

# get data
df3 <- pcoa_list$vectors[,1:2] %>% 
  as.data.frame() %>% 
  mutate(group = group_info) %>% 
  mutate(rel_eig = pcoa_list$values[,2]) %>% 
  mutate(id = id_info)

# plot
#library(ggh4x)   #stat_centroid
#library(ggalt)  #geom_encircle 
#library(ggpubr)  #link the points
# stat_centroid(aes(xend = axis1, yend = axis2, colour = group),geom = "segment", crop_other = F,alpha=0.2,size = 1,show.legend = F)+
# stat_chull(geom = "polygon", aes(group = group, fill = group), alpha = 0.1,show.legend = F)+
# geom_encircle(aes(fill=group), alpha = 0.05, show.legend = F,expand=0,spread=0.5,s_shape=0.8)+ 
library(ggpubr)

df3$group <- factor(df3$group,levels = c("Crop","Gizzard","Duodenum","Jejunum","Ileum","Cecum",    "Colon"))

p1 <- ggplot(df3,aes(Axis.1,Axis.2,group = group)) +
  geom_point(size=2,aes(color=group)) + 
  scale_color_manual(values = dark_col) +
  labs(x='PCoA1(41.06%)',y='PCoA2(16.57%)',color="",title = "Bray") +
  theme_bw() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(size = 8,face = "bold"),
    axis.title = element_text(size = 8),
    legend.text = element_text(size =8, color = "black"), 
    legend.title = element_blank()
  ) 
#stat_chull(geom = "polygon", aes(group = group, fill = group), alpha = 0.2,show.legend = F)+
#geom_text(aes(label = id), size = 3, vjust = 1)

p1

library(ggh4x)
p1 <- p1 + force_panelsizes(total_width=unit(1.7, "in"),total_height=unit(1.7, "in"))
p1 
ggsave("bray.pdf",p1)

## statistical analysis of b-diversity --------
set.seed(111)
aov_bray <- adonis2(dis_bray~group_info, permutations = 1000)
aov_jaccard <- adonis2(dis_jaccard~group_info, permutations = 1000)
# Export statistical analysis result (R square and P value)
aov_bray$source <- 'bray'
aov_jaccard$source <- 'jaccard'
# empty data.frame used to separate different result
empty <- data.frame(matrix(ncol = 6, nrow = 1))
colnames(empty) <- colnames(aov_bray)
rownames(empty) <- '####'
empty[] <- '####'

result_aov <- rbind(aov_bray,empty,aov_jaccard)
write.csv(result_aov,file = "premanova.csv")

# pairwise comparison
co_group <- combn(unique(group_info),m=2)
pair_bray <- data.frame(matrix(ncol = 6,nrow = ncol(co_group)))
colnames(pair_bray) <- c('type','pairs','R2','p.value','p.adjust','sig')
pair_jaccard <- data.frame(matrix(ncol = 6,nrow = ncol(co_group)))
colnames(pair_jaccard) <- c('type','pairs','R2','p.value','p.adjust','sig')
for (co in 1:ncol(co_group)){
  sub_data <- subset(df, seg %in% co_group[,co])
  sub_sample <- sub_data$seg
  sub_data <- sub_data[,-1] %>% 
    mutate_all(as.numeric)
  dis_bray <- vegdist(sub_data, method = "bray")
  dis_jaccard <- vegdist(sub_data, method = "jaccard")
  
  aov_jaccard <- adonis2(dis_jaccard~sub_sample, permutations = 1000)
  pair_jaccard$type = "jaccard"
  pair_jaccard$pairs[co] = paste(co_group[,co][1],co_group[,co][2],sep = " vs ")
  pair_jaccard$R2[co] = aov_jaccard$R2[1]
  pair_jaccard$p.value[co] = aov_jaccard$`Pr(>F)`[1]
  
  set.seed(111)
  aov_bray <- adonis2(dis_bray~sub_sample, permutations = 1000)
  pair_bray$type = "bray"
  pair_bray$pairs[co] = paste(co_group[,co][1],co_group[,co][2],sep = " vs ")
  pair_bray$R2[co] = aov_bray$R2[1]
  pair_bray$p.value[co] = aov_bray$`Pr(>F)`[1]
}

#adjust p value and get sig sign 
padj_function <- function(data){
  data$p.adjust <- p.adjust(data$p.value,method="fdr")
  for (i in 1:nrow(data)){
    if (data$p.adjust[i] < 0.05){
      data$sig[i] = "*"
      if (data$p.adjust[i] < 0.01){
        data$sig[i] = "**"
        if (data$p.adjust[i] < 0.001){
          data$sig[i] = "***"
        }
      }
    }
    else {
      data$sig[i] = " "
    }
  }
  return(data)
}

pair_bray <- padj_function(pair_bray) 
pair_jaccard <- padj_function(pair_jaccard) 

pairwise_info <- as.data.frame(rbind(pair_bray,pair_jaccard))

write.csv(pairwise_info,file = "B_diversity_pairwise_info.csv")

# taxonomy composition --------
library(tidyverse)
library(paletteer)
library(ggtext)
library(glue)

col1 <- paletteer_d("ggsci::category20_d3")
col2 <- paletteer_d("Polychrome::kelly")[7:15]
my_color <- c(col1[1:7],col2,col1[8:20])

# check color
data.frame(
  color_id = seq_along(my_color),
  color_code = my_color
) %>% ggplot(aes(x = factor(color_id), y = 1, fill = color_code)) +
  geom_bar(stat = "identity") +
  scale_fill_identity()

# individual composition
feature_table <- read.csv("otutab.csv", header = T,  row.names = 1) %>% 
  apply(.,2,function(x){
    x/sum(x)*100
  }) %>% 
  as.data.frame() %>% 
  rownames_to_column()

tax <- read.csv("taxonomy_kmer.csv",header = T,row.names = 1) %>% 
  rownames_to_column()

sample_data <- read.csv("metadata.csv",header = T,row.names = 1) %>% 
  rownames_to_column()

otu_tab <- left_join(tax,feature_table,by="rowname")

## phylum -------------
library(glue)
library(ggtext)
df2 <- otu_tab %>% 
  select(3,9:ncol(.)) %>% 
  group_by(Phylum) %>% 
  mutate(across(contains("barcode"), ~ replace(as.numeric(.), is.na(as.numeric(.)), 0))) %>% 
  summarise(across(where(is.numeric), sum)) %>% 
  column_to_rownames("Phylum") %>% 
  arrange(.,rowSums(.)) %>% 
  t() %>% 
  as.data.frame() %>% 
  select(Unassigned,everything(.)) %>% 
  rownames_to_column() %>% 
  left_join(sample_data,.,by="rowname") %>% 
  filter(.,!seg %in% c("Neg","Pos")) %>% 
  select(2,4:ncol(.)) %>% 
  pivot_longer(.,cols = -1, names_to = "bacteria", values_to = "abun") %>% 
  separate(bacteria,into=c("V1","V2") ,sep = "_") %>% 
  mutate(V2 = replace(V2, is.na(V2), "")) %>% 
  mutate(font_name=ifelse(V2=="",glue("<i>{V1}"),glue("<i>{V1}_</i>{V2}")))
df2$font_name <- gsub("<i>Unassigned","Unassigned",df2$font_name) 

df2$font_name <- factor(df2$font_name, levels = unique(df2$font_name))
df2$Sample <- factor(df2$Sample, levels = unique(df2$Sample))
col1 <- c(my_color[1:length(unique(df2$font_name))-1],"gray")

p1 <-  ggplot(df2, aes(Sample,abun)) +
  geom_col(aes(fill = font_name),color="black")+
  #scale_color_manual(values = my_color) +
  scale_fill_manual(values = rev(col1)) +
  labs(x="",y="Relative Abundance (%)",fill="Phylum") +
  scale_x_discrete(expand = c(0.01, 0)) +  # default 0.2
  scale_y_continuous(expand = c(0.005, 0)) +  # default 0.05
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 7,color = "black",angle = 90,hjust = 1),
    axis.text.y = element_text(size = 7,color = "black"),
    axis.title = element_text(size = 8,color = "black"),
    legend.text = element_markdown(size =7, color = "black"), 
    legend.title = element_text(size =7, color = "black"), 
    panel.grid = element_blank(),
    legend.key.size = unit(0.3, "cm"),
    panel.border = element_blank()
  )

p1
library(ggh4x)
p1 <- p1 + force_panelsizes(total_width=unit(7, "in"),total_height=unit(3, "in"))
p1 
ggsave("phylum_individual.pdf",p1)


## family -------------
df2 <- otu_tab %>% 
  select(6,9:ncol(.)) %>% 
  group_by(Family) %>% 
  mutate(across(contains("barcode"), ~ replace(as.numeric(.), is.na(as.numeric(.)), 0))) %>% 
  summarise(across(where(is.numeric), sum)) %>% 
  column_to_rownames("Family") %>% 
  arrange(.,rowSums(.)) %>%
  t() %>% 
  as.data.frame() %>% 
  select(Unassigned,(ncol(.)-9):ncol(.)) %>% 
  mutate(Others=100-rowSums(.)) %>% 
  select(Others,everything()) %>% 
  rownames_to_column() %>% 
  left_join(sample_data,.,by="rowname") %>% 
  filter(.,!seg %in% c("Neg","Pos")) %>% 
  select(2,4:ncol(.)) %>% 
  pivot_longer(.,cols = -1, names_to = "bacteria", values_to = "abun") %>% 
  separate(bacteria,into=c("V1","V2") ,sep = "_U") %>% 
  mutate(V2 = replace(V2, is.na(V2), "")) %>% 
  mutate(font_name=ifelse(V2=="",glue("<i>{V1}"),glue("<i>{V1}</i>_U{V2}")))

df2$font_name <- gsub("<i>Unassigned","Unassigned",df2$font_name) 
df2$font_name <- gsub("<i>Others","Others",df2$font_name) 

df2$font_name <- factor(df2$font_name, levels = unique(df2$font_name))
df2$Sample <- factor(df2$Sample, levels = unique(df2$Sample))
col1 <- c(my_color[1:(length(unique(df2$font_name))-2)],"gray","snow")

p1 <-  ggplot(df2, aes(Sample,abun)) +
  geom_col(aes(fill = font_name),color="black")+
  #scale_color_manual(values = my_color) +
  scale_fill_manual(values = rev(col1)) +
  labs(x="",y="Relative Abundance (%)",fill="Family") +
  scale_x_discrete(expand = c(0.01, 0)) +  # default 0.2
  scale_y_continuous(expand = c(0.005, 0)) +  # default 0.05
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 7,color = "black",angle = 90,hjust = 1),
    axis.text.y = element_text(size = 7,color = "black"),
    axis.title = element_text(size = 8,color = "black"),
    legend.text = element_markdown(size =7, color = "black"), 
    legend.title = element_text(size =7, color = "black"), 
    panel.grid = element_blank(),
    legend.key.size = unit(0.3, "cm"),
    panel.border = element_blank()
  )

p1
library(ggh4x)
p1 <- p1 + force_panelsizes(total_width=unit(7, "in"),total_height=unit(3, "in"))
p1 
ggsave("family_individual.pdf",p1)


## genus -------------
df2 <- otu_tab %>% 
  select(7,9:ncol(.)) %>% 
  group_by(Genus) %>% 
  mutate(across(contains("barcode"), ~ replace(as.numeric(.), is.na(as.numeric(.)), 0))) %>% 
  summarise(across(where(is.numeric), sum)) %>% 
  column_to_rownames("Genus") %>% 
  arrange(.,rowSums(.)) %>%
  t() %>% 
  as.data.frame() %>% 
  select(Unassigned,(ncol(.)-19):ncol(.)) %>% 
  mutate(Others=100-rowSums(.)) %>% 
  select(Others,everything()) %>% 
  rownames_to_column() %>% 
  left_join(sample_data,.,by="rowname") %>% 
  filter(.,!seg %in% c("Neg","Pos")) %>% 
  select(2,4:ncol(.)) %>% 
  pivot_longer(.,cols = -1, names_to = "bacteria", values_to = "abun") %>% 
  separate(bacteria,into=c("V1","V2") ,sep = "_U") %>% 
  mutate(V2 = replace(V2, is.na(V2), "")) %>% 
  mutate(font_name=ifelse(V2=="",glue("<i>{V1}"),glue("<i>{V1}</i>_U{V2}")))

df2$font_name <- gsub("<i>Unassigned","Unassigned",df2$font_name) 
df2$font_name <- gsub("<i>Others","Others",df2$font_name) 

df2$font_name <- factor(df2$font_name, levels = unique(df2$font_name))
df2$Sample <- factor(df2$Sample, levels = unique(df2$Sample))
col1 <- c(my_color[1:(length(unique(df2$font_name))-2)],"gray","snow")

p1 <-  ggplot(df2, aes(Sample,abun)) +
  geom_col(aes(fill = font_name),color="black")+
  #scale_color_manual(values = my_color) +
  scale_fill_manual(values = rev(col1)) +
  labs(x="",y="Relative Abundance (%)",fill="Family") +
  scale_x_discrete(expand = c(0.01, 0)) +  # default 0.2
  scale_y_continuous(expand = c(0.005, 0)) +  # default 0.05
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 7,color = "black",angle = 90,hjust = 1),
    axis.text.y = element_text(size = 7,color = "black"),
    axis.title = element_text(size = 8,color = "black"),
    legend.text = element_markdown(size =7, color = "black"), 
    legend.title = element_text(size =7, color = "black"), 
    panel.grid = element_blank(),
    legend.key.size = unit(0.3, "cm"),
    panel.border = element_blank()
  )+
  guides(color = guide_legend(ncol = 1), fill = guide_legend(ncol = 1))

p1
library(ggh4x)
p1 <- p1 + force_panelsizes(total_width=unit(7, "in"),total_height=unit(3, "in"))
p1 
ggsave("genus_individual.pdf",p1)


## species -------------
df2 <- otu_tab %>% 
  select(8,9:ncol(.)) %>% 
  group_by(Species) %>% 
  mutate(across(contains("barcode"), ~ replace(as.numeric(.), is.na(as.numeric(.)), 0))) %>% 
  summarise(across(where(is.numeric), sum)) %>% 
  column_to_rownames("Species") %>% 
  arrange(.,rowSums(.)) %>%
  t() %>% 
  as.data.frame() %>% 
  select(Unassigned,(ncol(.)-29):ncol(.)) %>% 
  mutate(Others=100-rowSums(.)) %>% 
  select(Others,everything()) %>% 
  rownames_to_column() %>% 
  left_join(sample_data,.,by="rowname") %>% 
  filter(.,!seg %in% c("Neg","Pos")) %>% 
  select(2,4:ncol(.)) %>% 
  pivot_longer(.,cols = -1, names_to = "bacteria", values_to = "abun") %>% 
  separate(bacteria,into=c("V1","V2") ,sep = "_U") %>% 
  mutate(V2 = replace(V2, is.na(V2), "")) %>% 
  mutate(font_name=ifelse(V2=="",glue("<i>{V1}"),glue("<i>{V1}</i>_U{V2}")))

df2$font_name <- gsub("<i>Unassigned","Unassigned",df2$font_name) 
df2$font_name <- gsub("<i>Others","Others",df2$font_name) 

df2$font_name <- factor(df2$font_name, levels = unique(df2$font_name))
df2$Sample <- factor(df2$Sample, levels = unique(df2$Sample))
col1 <- c(my_color[1:(length(unique(df2$font_name))-2)],"gray","snow")

p1 <-  ggplot(df2, aes(Sample,abun)) +
  geom_col(aes(fill = font_name),color="black")+
  #scale_color_manual(values = my_color) +
  scale_fill_manual(values = rev(col1)) +
  labs(x="",y="Relative Abundance (%)",fill="Family") +
  scale_x_discrete(expand = c(0.01, 0)) +  # default 0.2
  scale_y_continuous(expand = c(0.005, 0)) +  # default 0.05
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 7,color = "black",angle = 90,hjust = 1),
    axis.text.y = element_text(size = 7,color = "black"),
    axis.title = element_text(size = 8,color = "black"),
    legend.text = element_markdown(size =7, color = "black"), 
    legend.title = element_text(size =7, color = "black"), 
    panel.grid = element_blank(),
    legend.key.size = unit(0.3, "cm"),
    panel.border = element_blank()
  )+
  guides(color = guide_legend(ncol = 1), fill = guide_legend(ncol = 1))

p1
library(ggh4x)
p1 <- p1 + force_panelsizes(total_width=unit(7, "in"),total_height=unit(3, "in"))
p1 
ggsave("species_individual.pdf",p1)





# color
library(paletteer)
col1 <- paletteer_d("khroma::soil")
col1 <- paletteer_d("Polychrome::kelly")
col1 <- paletteer_d("ggthemes::Classic_20")
col1 <- paletteer_d("ggsci::category20_d3")

data.frame(
  color_id = seq_along(col1),
  color_code = col1
) %>% ggplot(aes(x = factor(color_id), y = 1, fill = color_code)) +
  geom_bar(stat = "identity") +
  scale_fill_identity()+labs(title = "ggsci::category20_d3")



#
col1 <- paletteer_d("ggsci::category20_d3")
col2 <- paletteer_d("Polychrome::kelly")[7:15]
my_color <- c(col1[1:7],col2,col1[8:20])

p <-  ggplot(fam_data,aes(x=Sample,y=value,fill=Species)) +
  geom_col() + 
  facet_wrap(~seg,nrow=1,scales = "free_x") +
  labs(x = NULL, y = 'Relative Abundance (%)',fill='Species') +
  scale_x_discrete(expand = c(0, 0))+
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = rev(c(my_color[1:15],"gray"))) +  
  theme_bw()+
  guides(fill = guide_legend(ncol = 1))+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(linewidth = 0.5,color = 'white'),
        #axis.text.x = element_text(size = 7,angle = 0,vjust = 1,hjust = 0.5,color = "black"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 7,color = "black"),
        axis.title.y = element_text(margin = margin(r = 4), size = 9, color = "black"), 
        axis.ticks.x = element_blank(),  
        axis.ticks.y = element_blank(),  
        legend.key.size = unit(0.3, "cm"),  
        legend.text = element_markdown(size =10, color = "black"), 
        legend.margin = margin(l=-5),
        legend.spacing.x = unit(0.02, "cm"), 
        legend.spacing.y = unit(0.02, "cm"), 
        strip.background = element_rect(linetype = "solid",color = "white", linewidth = 0.5),
        strip.text = element_text(size = 9),
        panel.spacing.x= unit(0.15, "lines")
  ) 

ggsave("p.png", p ,width = 4,height = 5,dpi=600)





# epi2me-----------------------
setwd("D:/OneDrive - Oklahoma A and M System/tranfer documents/16s_fulength/epi2me")
setwd("C:/Users/jiaguo/OneDrive - Oklahoma A and M System/tranfer documents/16s_fulength/epi2me")

library(tidyverse)

sample_data <- read.csv("metadata.csv",header = T,row.names = 1) %>% 
  rownames_to_column()

feature_table <- read.csv("feature_table.csv", header = T,  row.names = 1) %>% 
  rownames_to_column() %>% 
  separate(rowname,into = c("Domain","sub_domain","Phylum","Class","Order","Family","Genus","Species"),sep = ";")



# a diversity------------------
library(tidyverse)
library(vegan)
library(rstatix)
library(ggpubr)
library(car) # equal variance

sample_data <- read.csv("metadata.csv",header = T,row.names = 1) %>% 
  rownames_to_column()

feature_table <- read.csv("feature_table.csv", header = T,  row.names = 1) %>% 
  rownames_to_column() %>% 
  separate(rowname,into = c("Domain","sub_domain","Phylum","Class","Order","Family","Genus","Species"),sep = ";") %>% 
  select(-c(1:7)) %>% 
  column_to_rownames(var = "Species") %>% 
  apply(.,2,function(x){
    x/sum(x)
  }) %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  left_join(.,sample_data,by="rowname") %>% 
  select((ncol(.)-1),ncol(.),everything())

otu_tab <- feature_table

# cal a-diversity
group_info <- otu_tab$seg
id_info <- otu_tab$Sample

df2 <- otu_tab %>% 
  select(-c(1:3))

shan_index <- data.frame(id=id_info,group=group_info,shannon=diversity(df2, index = "shannon")) 

obs_spec <- data.frame(id=id_info,group=group_info,obs=specnumber(df2))

pielou <- data.frame(id=id_info,group=group_info,pielou=shan_index$shannon / log(obs_spec$obs)) 

a_diversity <- left_join(shan_index,obs_spec[,c(1,3)],by = "id") %>% 
  left_join(.,pielou[,c(1,3)],by="id") %>% 
  filter(!group %in% c("Neg","Pos"))

a_diversity$group <- factor(a_diversity$group,levels = c("Crop","Gizzard","Duodenum","Jejunum","Ileum","Cecum",    "Colon"))


# Normality test and Homogeneity of variance test
check_function <- function(formula, data){
  fit <- aov(formula, data = data)
  res <- residuals(fit)
  par(mfrow=c(1,2))
  qqnorm(res); qqline(res)
  cat("\nNormality test:\n")
  print(shapiro.test(res))
  cat("\nHomogeneity of variance test:\n")
  plot(a_diversity$group,res)
  abline(h=0)
  # using car package
  print(leveneTest(formula,data=data))
  par(mfrow=c(1,1))
}

# statistical analysis
kw_function <- function(formula, data){
  kw_test <- kruskal.test(formula,data = data)
  print(kw_test)
  posthoc <- dunn_test(data = data,formula,p.adjust.method = "fdr")
  posthoc$compar <- paste(posthoc$group1,"-",posthoc$group2,sep="")
  print(posthoc)
  return(posthoc)
}

aov_function <- function(formula, data){
  aov <- aov(formula,data = data)
  print(summary(aov))
  posthoc <- as.data.frame(TukeyHSD(aov)[[1]])
  posthoc$compar <- rownames(posthoc)
  colnames(posthoc)[colnames(posthoc) == "p adj"] <- "p.adj"
  posthoc<-separate(data = posthoc,col = compar,into = c("group1", "group2"),sep = "-")
  posthoc$compar <- paste(posthoc$group2,"-",posthoc$group1,sep = "")
  print(posthoc)
  return(posthoc)
}



compar <- shannon ~ group

check_function(compar, a_diversity)

posthoc <- kw_function(compar, data=a_diversity)

posthoc <- aov_function(compar, data=a_diversity)


# get the compact letter display for lists of comparisons
library(rcompanion)
max <- a_diversity %>% 
  group_by(group) %>% 
  summarise(max=max(shannon)+max(shannon)*0.01)
names(max)[1] <- "Group"
sig_mark <- cldList(p.adj ~ compar, data = posthoc) %>% 
  left_join(.,max, by = "Group")
sig_mark

# plot
library(paletteer)
library(ggnewscale)
col1 <- paletteer_d("ggsci::category20_d3")
light_col <- col1[c(14,12,19,13,20,18,15)]
dark_col <- col1[c(4,2,9,3,1,8,5)]

# Shannon Index
# Observed Species
# Pielou's Eveness
p1 <- ggplot(a_diversity,aes(x = group,y = shannon)) +
  geom_boxplot(aes(colour = group,fill = group),outliers = F) +
  scale_color_manual(values = dark_col) +
  scale_fill_manual(values = light_col) +
  #geom_text(aes(label = id),position = position_jitter(width = 0.1), size = 3, vjust = -0.5) +
  labs(x="",y="Shannon Index") +
  scale_x_discrete(expand = c(0.1, 0)) +  # default 0.2
  scale_y_continuous(expand = c(0.1, 0)) +  # default 0.05
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 8,color = "black",angle = 45,hjust = 1),
    axis.text.y = element_text(size = 8,color = "black"),
    axis.title = element_text(size = 9,color = "black"),
    legend.text = element_text(size =9, color = "black"), 
    legend.title = element_blank(),
    legend.position = "none"
  )+
  geom_text(aes(label = Letter, x=Group, y =(max)), data=sig_mark,vjust = -0.5,size=3)+
  new_scale_fill() +
  geom_jitter(aes(colour = group,fill = group),shape = 21, width = 0.1,size=0.5) +
  scale_fill_manual(values = dark_col)

p1
library(ggh4x)
p1 <- p1 + force_panelsizes(total_width=unit(1.7, "in"),total_height=unit(1.7, "in"))
p1 
ggsave("shannon.pdf",p1)


# b diversity------------------
library(tidyverse)
library(vegan)
library(ape)

sample_data <- read.csv("metadata.csv",header = T,row.names = 1) %>% 
  rownames_to_column()

feature_table <- read.csv("feature_table.csv", header = T,  row.names = 1) %>% 
  rownames_to_column() %>% 
  separate(rowname,into = c("Domain","sub_domain","Phylum","Class","Order","Family","Genus","Species"),sep = ";") %>% 
  select(-c(1:7)) %>% 
  column_to_rownames(var = "Species") %>% 
  apply(.,2,function(x){
    x/sum(x)
  }) %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  left_join(.,sample_data,by="rowname") %>% 
  select((ncol(.)-1),ncol(.),everything()) %>% 
  filter(!seg %in% c("Neg","Pos"))

otu_tab <- feature_table

group_info <- otu_tab$seg
id_info <- otu_tab$Sample

df2 <- otu_tab %>% 
  select(1,4:ncol(.)) %>% 
  column_to_rownames(var = "Sample")

# pcoa analysis
dis_bray <- vegdist(df2, method = "bray")
dis_jaccard <- vegdist(df2, method = "jaccard")

pcoa_list <- pcoa(dis_jaccard)

# get data
df3 <- pcoa_list$vectors[,1:2] %>% 
  as.data.frame() %>% 
  mutate(group = group_info) %>% 
  mutate(rel_eig = pcoa_list$values[,2]) %>% 
  mutate(id = id_info)

# plot
#library(ggh4x)   #stat_centroid
#library(ggalt)  #geom_encircle 
#library(ggpubr)  #link the points
# stat_centroid(aes(xend = axis1, yend = axis2, colour = group),geom = "segment", crop_other = F,alpha=0.2,size = 1,show.legend = F)+
# stat_chull(geom = "polygon", aes(group = group, fill = group), alpha = 0.1,show.legend = F)+
# geom_encircle(aes(fill=group), alpha = 0.05, show.legend = F,expand=0,spread=0.5,s_shape=0.8)+ 
library(ggpubr)

df3$group <- factor(df3$group,levels = c("Crop","Gizzard","Duodenum","Jejunum","Ileum","Cecum",    "Colon"))

# Bray–Curtis
# Jaccard
p1 <- ggplot(df3,aes(Axis.1,Axis.2,group = group)) +
  geom_point(size=2,aes(color=group)) + 
  scale_color_manual(values = dark_col) +
  labs(x='PCoA1(67.4%)',y='PCoA2(10.3%)',color="",title = "Jaccard") +
  theme_bw() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(size = 8,face = "bold"),
    axis.title = element_text(size = 8),
    legend.text = element_text(size =8, color = "black"), 
    legend.title = element_blank()
  ) 
#stat_chull(geom = "polygon", aes(group = group, fill = group), alpha = 0.2,show.legend = F)+
#geom_text(aes(label = id), size = 3, vjust = 1)

p1

library(ggh4x)
p1 <- p1 + force_panelsizes(total_width=unit(1.7, "in"),total_height=unit(1.7, "in"))
p1 
ggsave("Jaccard.pdf",p1)

## statistical analysis of b-diversity --------
set.seed(111)
aov_bray <- adonis2(dis_bray~group_info, permutations = 1000)
set.seed(111)
aov_jaccard <- adonis2(dis_jaccard~group_info, permutations = 1000)
# Export statistical analysis result (R square and P value)
aov_bray$source <- 'bray'
aov_jaccard$source <- 'jaccard'
# empty data.frame used to separate different result
empty <- data.frame(matrix(ncol = 6, nrow = 1))
colnames(empty) <- colnames(aov_bray)
rownames(empty) <- '####'
empty[] <- '####'

result_aov <- rbind(aov_bray,empty,aov_jaccard)
write.csv(result_aov,file = "premanova.csv")

# pairwise comparison
co_group <- combn(unique(group_info),m=2)
pair_bray <- data.frame(matrix(ncol = 6,nrow = ncol(co_group)))
colnames(pair_bray) <- c('type','pairs','R2','p.value','p.adjust','sig')
pair_jaccard <- data.frame(matrix(ncol = 6,nrow = ncol(co_group)))
colnames(pair_jaccard) <- c('type','pairs','R2','p.value','p.adjust','sig')
for (co in 1:ncol(co_group)){
  sub_data <- subset(otu_tab, seg %in% co_group[,co])
  sub_sample <- sub_data$seg
  sub_data <- sub_data %>% 
    select(1,4:ncol(.)) %>% 
    remove_rownames() %>%
    column_to_rownames(.,var = "Sample")
  dis_bray <- vegdist(sub_data, method = "bray")
  dis_jaccard <- vegdist(sub_data, method = "jaccard")
  
  set.seed(111)
  aov_jaccard <- adonis2(dis_jaccard~sub_sample, permutations = 1000)
  pair_jaccard$type = "jaccard"
  pair_jaccard$pairs[co] = paste(co_group[,co][1],co_group[,co][2],sep = " vs ")
  pair_jaccard$R2[co] = aov_jaccard$R2[1]
  pair_jaccard$p.value[co] = aov_jaccard$`Pr(>F)`[1]
  
  set.seed(111)
  aov_bray <- adonis2(dis_bray~sub_sample, permutations = 1000)
  pair_bray$type = "bray"
  pair_bray$pairs[co] = paste(co_group[,co][1],co_group[,co][2],sep = " vs ")
  pair_bray$R2[co] = aov_bray$R2[1]
  pair_bray$p.value[co] = aov_bray$`Pr(>F)`[1]
}

#adjust p value and get sig sign 
padj_function <- function(data){
  data$p.adjust <- p.adjust(data$p.value,method="fdr")
  for (i in 1:nrow(data)){
    if (data$p.adjust[i] < 0.05){
      data$sig[i] = "*"
      if (data$p.adjust[i] < 0.01){
        data$sig[i] = "**"
        if (data$p.adjust[i] < 0.001){
          data$sig[i] = "***"
        }
      }
    }
    else {
      data$sig[i] = " "
    }
  }
  return(data)
}

pair_bray <- padj_function(pair_bray) 
pair_jaccard <- padj_function(pair_jaccard) 

pairwise_info <- as.data.frame(rbind(pair_bray,pair_jaccard))

write.csv(pairwise_info,file = "B_diversity_pairwise_info.csv")

# EMU_NEW_database---------------------------
## raw data combine ----------------
setwd("D:/OneDrive - Oklahoma A and M System/tranfer documents/16s_fulength/emu_new/raw_tax")

library(tidyverse)

# modify file name
ori_name <- list.files(full.names = FALSE)
file_names <- strsplit(ori_name, "_")
file_names <- sapply(file_names, function(x) x[1])
file_names <- paste0(file_names,".tsv")
file.copy(ori_name, file_names)

# combine
file_paths <- list.files(full.names = T)

data_list <- lapply(file_paths, function(file) {
  read.table(file, sep = "\t",header = T)
})

file_paths <- strsplit(basename(file_paths),"\\.")
file_paths <- sapply(file_paths,function(x){x[1]})
names(data_list) <- file_paths

data_list1 <- lapply(data_list, function(file) {
  file <- file[1:(nrow(file)-2),]
  file <- file[,-c(9,10,11,12,13)]
})

data_list2 <- lapply(data_list1, function(file) {
  file %>%
    unite("tax", c("species", "genus", "family", "order", "class", "phylum"), sep = ";")
})


data_list3 <- lapply(seq_along(data_list2), function(x) {
  df <- data_list2[[x]]
  names(df)[2] <- file_paths[x]
  df
})

names(data_list3) <- names(data_list2)

data_list4 <- lapply(data_list3,function(file){
  file <- file[,-1]
})

merged_df <- reduce(data_list4, full_join, by = "tax")

merged_df[is.na(merged_df)] <- 0

df <- merged_df %>% 
  select(tax,everything(.)) %>% 
  separate(.,tax,c("species","genus","family", "order", "class", "phylum"),sep=";")


df2 <- as.data.frame(t(df)) %>% 
  mutate(id=rownames(.)) %>% 
  select(ncol(.),everything(.))


metadata <- read.csv("metadata.csv",header = T)
names(metadata)[1] <- "id"
metadata <- metadata[,-2]

df3 <- left_join(df2, metadata,by = "id") %>% 
  select(.,c(1,ncol(.),everything(.)))

df4 <- as.data.frame(t(df3))  

write.csv(df4, "feature_table.csv")

# rarefy curve ----------

library(tidyverse)
library(vegan)


df <- read.csv("feature_table.csv", header = TRUE, skip = 1, row.names = 1) %>%
  select(-(1:6)) %>%
  t() %>%
  as.data.frame() %>% 
  slice(-(68:74))

group_info <- df[,1]

df <- df[,-1]

df <- apply(df,2,function(x){
  x = as.numeric(x)
  round(x*10000)
}) %>% 
  as.data.frame()

df <- df %>% 
  mutate(group = group_info) %>% 
  select(group,everything())

## rarefy sequencing depth---------------
# data: sample(row)*otu   first column: group
# define function
my_fun <- function(data,group){
  # subset data and rarefy
  data1 <- df[df$group == group,] %>% 
    select(-1) 
  
  # rarefy
  depths <- c(seq(1,10000,500))
  data1_list <- list()
  for (i in seq_along(depths)){
    data1_list[[i]] <- rrarefy(data1, depths[i])
  }
  names(data1_list) <- depths
  
  # cal shannon
  #shan1 <- sapply(data1_list, function(df){
  #  diversity(df, index = "shannon")
  #}) %>% 
  #  as.data.frame() 
  
  # cal specnumber
  shan1 <- sapply(data1_list, function(df){
    specnumber(df)
  }) %>% 
    as.data.frame() 
  
  # cal mean and se
  shan2 <- apply(shan1,2,function(x){
    c(mean = mean(x), se = sd(x) / sqrt(length(x)))
  })%>%
    t() %>% 
    as.data.frame() %>% 
    rownames_to_column() %>% 
    mutate(group = group)
  
  # return data
  return(shan2)
}

# get data
shanlist <- list()
shandata <- list()
iteration_seq <- paste0("iteration",seq(1,20))
for (i in 1:length(iteration_seq)){
  for (j in unique(df$group)){
    shanlist[[j]] <- my_fun(df,j)
  }
  shandata[[i]] <- do.call("rbind", shanlist)
}
names(shandata) <- iteration_seq

for (i in names(shandata)){
  shandata[[i]]$iteration <- i
}


shandata <- do.call("rbind", shandata)


pic_df <- shandata %>% 
  group_by(rowname,group) %>% 
  summarise(mean_f <- mean(mean))
names(pic_df)[3] <- "mean"


# insert value (spline)
pic_df2 <- pic_df 
re_list <- list()
for (i in unique(pic_df2$group)){
  re_list[[i]] <- pic_df2[pic_df2$group == i,]
  spline_temp <- spline(re_list[[i]]$rowname,re_list[[i]]$mean,n=10*nrow(re_list[[i]]))
  re_list[[i]] <- data.frame(
    rowname = spline_temp$x,
    mean = spline_temp$y,
    group = i
  )
}

pic_df2 <- do.call("rbind", re_list)
pic_df2$group <- factor(pic_df2$group,levels = c("Crop","Gizzard","Duodenum","Jejunum","Ileum","Cecum",    "Colon"))


# plot
library(paletteer)
col1 <- paletteer_d("ggsci::category20_d3")
light_col <- col1[c(14,12,19,13,20,18,15)]
dark_col <- col1[c(4,2,9,3,1,8,5)]

p <- ggplot(pic_df2,aes(x =as.numeric(rowname),y=mean,group = group)) +
  geom_line(linewidth=0.5,aes(color=group)) +
  scale_color_manual(values = dark_col) +
  labs(x="Sequencing Depth",y="Observed Species")+
  scale_x_continuous(expand = c(0.05, 0)) +
  scale_y_continuous(expand = c(0.05, 0)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 8,color = "black"),
    axis.text.y = element_text(size = 8,color = "black"),
    axis.title = element_text(size = 9,color = "black"),
    legend.text = element_text(size =9, color = "black"), 
    legend.title = element_blank(),
  )
p
ggsave("obs_depth.pdf", p ,width = 4,height = 2.5,dpi=600)



## rarefy sample---------------

# define function for rarefy
my_fun2 <- function(data){
  temp_list <- list()
  result_list <- list()
  for (i in unique(df$group)){
    temp_list[[i]] <- df[df$group == i,]
    temp_data <- temp_list[[i]]
    result_list[[i]] <- matrix(NA,nrow = nrow(temp_data)+1,ncol = 3)
    for (j in 1:nrow(temp_data)){
      result_list[[i]][1,1] = 0
      result_list[[i]][1,2] = 0
      result_list[[i]][,3] = i
      result_list[[i]][j+1,1] = j 
      temp_data2 <-colSums(temp_data[sample(1:nrow(temp_data),j),-1])
      result_list[[i]][j+1,2] = specnumber(temp_data2)
    }
  }
  return(result_list)
}

# iteration 
iteration_seq <- paste0("iteration",seq(1,400))
result_list <- list()
for (i in 1:length(iteration_seq)){
  temp1 <- as.data.frame(do.call("rbind", my_fun2(df)))
  names(temp1) <- c("x","y","group")
  temp1$iteration = iteration_seq[i]
  result_list[[i]] = temp1
}
names(result_list) <- iteration_seq

pic_df <- as.data.frame(do.call("rbind",result_list)) 
pic_df$x <- as.numeric(pic_df$x)
pic_df$y <- as.numeric(pic_df$y)
pic_df <- pic_df %>% 
  group_by(x,group) %>% 
  summarise(mean=mean(y))

pic_df$group <- factor(pic_df$group,levels = c("Crop","Gizzard","Duodenum","Jejunum","Ileum","Cecum",    "Colon"))

# insert value (spline)
pic_df2 <- pic_df 
re_list <- list()
for (i in unique(pic_df2$group)){
  re_list[[i]] <- pic_df2[pic_df2$group == i,]
  spline_temp <- spline(re_list[[i]]$x,re_list[[i]]$mean,n=10*nrow(re_list[[i]]))
  re_list[[i]] <- data.frame(
    x = spline_temp$x,
    y = spline_temp$y,
    group = i
  )
}

pic_df2 <- do.call("rbind", re_list)
pic_df2$group <- factor(pic_df2$group,levels = c("Crop","Gizzard","Duodenum","Jejunum","Ileum","Cecum",    "Colon"))






# plot
library(paletteer)
col1 <- paletteer_d("ggsci::category20_d3")
light_col <- col1[c(14,12,19,13,20,18,15)]
dark_col <- col1[c(4,2,9,3,1,8,5)]

p <- ggplot(pic_df,aes(x = as.numeric(x),y = as.numeric(mean),group = group)) +
  geom_line(linewidth=0.5,aes(color=group)) +
  scale_color_manual(values = dark_col) +
  labs(x="Number of Individuals",y="Observed Species")+
  scale_x_continuous(expand = c(0.05, 0),limits = c(0,12),breaks = seq(0,12,2)) +
  scale_y_continuous(expand = c(0.05, 0)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 8,color = "black"),
    axis.text.y = element_text(size = 8,color = "black"),
    axis.title = element_text(size = 9,color = "black"),
    legend.text = element_text(size =9, color = "black"), 
    legend.title = element_blank(),
  )
p
ggsave("obs_sample.pdf", p ,width = 4,height = 2.5,dpi=600)

# upset venn ------------------
library(tidyverse)
library(UpSetR)


df <- read.csv("feature_table.csv", header = TRUE, skip = 1, row.names = 1) %>%
  select(-(1:6)) %>%
  t() %>%
  as.data.frame() %>% 
  slice(-(68:74))

list1 <- list()
for (i in unique(df$seg)) {
  list1[[i]] <- df[df$seg == i,] %>% 
    select(-1) %>%         
    mutate_all(as.numeric) %>%  # Convert all columns to numeric
    colSums() 
}

pic_df <- do.call("rbind",list1) %>% 
  t() %>% 
  apply(.,c(1,2),function(x){
    ifelse(x > 0, 1, 0)
  }) %>% 
  as.data.frame() %>% 
  rownames_to_column()

# 5*7.5 inch
upset(pic_df,nsets = 7,keep.order = TRUE,set_size.show = TRUE, 
      sets = c("Colon","Cecum","Ileum","Jejunum","Duodenum","Gizzard","Crop"),
      mainbar.y.label = "Intersect Species", sets.x.label = "Observed Species")


# a diversity------------------
library(tidyverse)
library(vegan)
library(rstatix)
library(ggpubr)
library(car) # equal variance

setwd("D:/OneDrive - Oklahoma A and M System/tranfer documents/16s_fulength/emu_new")

df <- read.csv("feature_table.csv", header = TRUE, skip = 1, row.names = 1) %>%
  select(-(1:6)) %>%
  t() %>%
  as.data.frame() %>% 
  slice(-(68:74))

group_info <- df$seg

df2 <- df[,-1] %>% 
  mutate_all(as.numeric)

shan_index <- data.frame(group=group_info,shannon=diversity(df2, index = "shannon")) %>% 
  mutate(id = paste0(group,seq(1,nrow(.))))

obs_spec <- data.frame(group=group_info,obs=specnumber(df2)) %>% 
  mutate(id = paste0(group,seq(1,nrow(.))))

pielou <- data.frame(group=group_info,pielou=shan_index$shannon / log(obs_spec$obs)) %>% 
  mutate(id = paste0(group,seq(1,nrow(.))))

a_diversity <- left_join(shan_index,obs_spec[,2:3],by = "id") %>% 
  left_join(.,pielou[,2:3],by="id")

a_diversity$group <- factor(a_diversity$group,levels = c("Crop","Gizzard","Duodenum","Jejunum","Ileum","Cecum",    "Colon"))


# Normality test and Homogeneity of variance test
check_function <- function(formula, data){
  fit <- aov(formula, data = data)
  res <- residuals(fit)
  par(mfrow=c(1,2))
  qqnorm(res); qqline(res)
  cat("\nNormality test:\n")
  print(shapiro.test(res))
  cat("\nHomogeneity of variance test:\n")
  plot(a_diversity$group,res)
  abline(h=0)
  # using car package
  print(leveneTest(formula,data=data))
  par(mfrow=c(1,1))
}

# statistical analysis
kw_function <- function(formula, data){
  kw_test <- kruskal.test(formula,data = data)
  print(kw_test)
  posthoc <- dunn_test(data = data,formula,p.adjust.method = "fdr")
  posthoc$compar <- paste(posthoc$group1,"-",posthoc$group2,sep="")
  print(posthoc)
  return(posthoc)
}

aov_function <- function(formula, data){
  aov <- aov(formula,data = data)
  print(summary(aov))
  posthoc <- as.data.frame(TukeyHSD(aov)[[1]])
  posthoc$compar <- rownames(posthoc)
  colnames(posthoc)[colnames(posthoc) == "p adj"] <- "p.adj"
  posthoc<-separate(data = posthoc,col = compar,into = c("group1", "group2"),sep = "-")
  posthoc$compar <- paste(posthoc$group2,"-",posthoc$group1,sep = "")
  print(posthoc)
  return(posthoc)
}


compar <- obs ~ group

check_function(compar, a_diversity)

posthoc <- kw_function(compar, data=a_diversity)

posthoc <- aov_function(compar, data=a_diversity)


# get the compact letter display for lists of comparisons
library(rcompanion)
max <- a_diversity %>% 
  group_by(group) %>% 
  summarise(max=max(obs)+max(obs)*0.01)
names(max)[1] <- "Group"
sig_mark <- cldList(p.adj ~ compar, data = posthoc) %>% 
  left_join(.,max, by = "Group")
sig_mark

# plot
library(paletteer)
col1 <- paletteer_d("ggsci::category20_d3")
light_col <- col1[c(14,12,19,13,20,18,15)]
dark_col <- col1[c(4,2,9,3,1,8,5)]

# Shannon Index
# Observed Species
# Pielou's Eveness
p1 <- ggplot(a_diversity,aes(x = group,y = obs)) +
  geom_boxplot(aes(colour = group,fill = group),outliers = F) +
  scale_color_manual(values = dark_col) +
  scale_fill_manual(values = light_col) +
  geom_jitter(aes(colour = group,fill = group),shape=21,width = 0.1) +
  #geom_text(aes(label = id),position = position_jitter(width = 0.1), size = 3, vjust = -0.5) +
  labs(x="",y="Observed Species") +
  scale_x_discrete(expand = c(0.1, 0)) +  # default 0.2
  scale_y_continuous(expand = c(0.1, 0)) +  # default 0.05
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 8,color = "black",angle = 45,hjust = 1),
    axis.text.y = element_text(size = 8,color = "black"),
    axis.title = element_text(size = 9,color = "black"),
    legend.text = element_text(size =9, color = "black"), 
    legend.title = element_blank(),
    legend.position = "none"
  )+
  geom_text(aes(label = Letter, x=Group, y =(max)), data=sig_mark,vjust = -0.5,size=3)

p1
library(ggh4x)
p1 <- p1 + force_panelsizes(total_width=unit(1.7, "in"),total_height=unit(1.7, "in"))
p1 
ggsave("observed.pdf",p1)


# b diversity------------------
library(tidyverse)
library(vegan)
library(ape)

df <- read.csv("feature_table.csv", header = TRUE, skip = 1, row.names = 1) %>%
  select(-(1:6)) %>%
  t() %>%
  as.data.frame() %>% 
  slice(-(68:74))

group_info <- df$seg

df2 <- df[,-1] %>% 
  mutate_all(as.numeric)

# pcoa analysis
dis_bray <- vegdist(df2, method = "bray")
dis_jaccard <- vegdist(df2, method = "jaccard")

pcoa_list <- pcoa(dis_bray)

# get data
df3 <- pcoa_list$vectors[,1:2] %>% 
  as.data.frame() %>% 
  mutate(group = group_info) %>% 
  mutate(rel_eig = pcoa_list$values[,2]) %>% 
  mutate(id = paste0(group,seq(1,nrow(.))))

# plot
#library(ggh4x)   #stat_centroid
#library(ggalt)  #geom_encircle 
#library(ggpubr)  #link the points
# stat_centroid(aes(xend = axis1, yend = axis2, colour = group),geom = "segment", crop_other = F,alpha=0.2,size = 1,show.legend = F)+
# stat_chull(geom = "polygon", aes(group = group, fill = group), alpha = 0.1,show.legend = F)+
# geom_encircle(aes(fill=group), alpha = 0.05, show.legend = F,expand=0,spread=0.5,s_shape=0.8)+ 
library(ggpubr)

df3$group <- factor(df3$group,levels = c("Crop","Gizzard","Duodenum","Jejunum","Ileum","Cecum",    "Colon"))

p1 <- ggplot(df3,aes(Axis.1,Axis.2,group = group)) +
  geom_point(size=2,aes(color=group)) + 
  scale_color_manual(values = dark_col) +
  labs(x='PCoA1 (77.3%)',y='PCoA2 (8.3%)',color="",title = "Bray Curtis") +
  theme_bw() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(size = 8, face = "bold"),  
    axis.title = element_text(size = 8, face = "plain"),
    legend.text = element_text(size =8, color = "black"), 
    legend.title = element_blank()
  ) 
#stat_chull(geom = "polygon", aes(group = group, fill = group), alpha = 0.2,show.legend = F)+
#geom_text(aes(label = id), size = 3, vjust = 1)

p1

library(ggh4x)
p1 <- p1 + force_panelsizes(total_width=unit(1.7, "in"),total_height=unit(1.7, "in"))
p1 
ggsave("Bray.pdf",p1)

## statistical analysis of b-diversity --------
set.seed(111)
aov_bray <- adonis2(dis_bray~group_info, permutations = 1000)
aov_jaccard <- adonis2(dis_jaccard~group_info, permutations = 1000)
# Export statistical analysis result (R square and P value)
aov_bray$source <- 'bray'
aov_jaccard$source <- 'jaccard'
# empty data.frame used to separate different result
empty <- data.frame(matrix(ncol = 6, nrow = 1))
colnames(empty) <- colnames(aov_bray)
rownames(empty) <- '####'
empty[] <- '####'

result_aov <- rbind(aov_bray,empty,aov_jaccard)
write.csv(result_aov,file = "premanova.csv")

# pairwise comparison
co_group <- combn(unique(group_info),m=2)
pair_bray <- data.frame(matrix(ncol = 6,nrow = ncol(co_group)))
colnames(pair_bray) <- c('type','pairs','R2','p.value','p.adjust','sig')
pair_jaccard <- data.frame(matrix(ncol = 6,nrow = ncol(co_group)))
colnames(pair_jaccard) <- c('type','pairs','R2','p.value','p.adjust','sig')
for (co in 1:ncol(co_group)){
  sub_data <- subset(df, seg %in% co_group[,co])
  sub_sample <- sub_data$seg
  sub_data <- sub_data[,-1] %>% 
    mutate_all(as.numeric)
  dis_bray <- vegdist(sub_data, method = "bray")
  dis_jaccard <- vegdist(sub_data, method = "jaccard")
  
  aov_jaccard <- adonis2(dis_jaccard~sub_sample, permutations = 1000)
  pair_jaccard$type = "jaccard"
  pair_jaccard$pairs[co] = paste(co_group[,co][1],co_group[,co][2],sep = " vs ")
  pair_jaccard$R2[co] = aov_jaccard$R2[1]
  pair_jaccard$p.value[co] = aov_jaccard$`Pr(>F)`[1]
  
  set.seed(111)
  aov_bray <- adonis2(dis_bray~sub_sample, permutations = 1000)
  pair_bray$type = "bray"
  pair_bray$pairs[co] = paste(co_group[,co][1],co_group[,co][2],sep = " vs ")
  pair_bray$R2[co] = aov_bray$R2[1]
  pair_bray$p.value[co] = aov_bray$`Pr(>F)`[1]
}

#adjust p value and get sig sign 
padj_function <- function(data){
  data$p.adjust <- p.adjust(data$p.value,method="fdr")
  for (i in 1:nrow(data)){
    if (data$p.adjust[i] < 0.05){
      data$sig[i] = "*"
      if (data$p.adjust[i] < 0.01){
        data$sig[i] = "**"
        if (data$p.adjust[i] < 0.001){
          data$sig[i] = "***"
        }
      }
    }
    else {
      data$sig[i] = " "
    }
  }
  return(data)
}

pair_bray <- padj_function(pair_bray) 
pair_jaccard <- padj_function(pair_jaccard) 

pairwise_info <- as.data.frame(rbind(pair_bray,pair_jaccard))

write.csv(pairwise_info,file = "B_diversity_pairwise_info.csv")

# taxonomy composition --------
library(tidyverse)
library(paletteer)
library(ggtext)
library(glue)

col1 <- paletteer_d("ggsci::category20_d3")
col2 <- paletteer_d("Polychrome::kelly")[7:15]
my_color <- c(col1[1:7],col2,col1[8:20])

# check color
data.frame(
  color_id = seq_along(my_color),
  color_code = my_color
) %>% ggplot(aes(x = factor(color_id), y = 1, fill = color_code)) +
  geom_bar(stat = "identity") +
  scale_fill_identity()

barplot(rep(1,length(my_color))~seq(1,length(my_color)),col=my_color)

# total composition
df <- read.csv("feature_table.csv", header = TRUE, skip = 1, row.names = 1) %>% 
  select(-(74:80))

## phylum -------------
df2 <- df %>% 
  select(6:ncol(.))
group_info <- df2[1,]

df3 <- df2 %>% 
  slice(-1) %>% 
  group_by(phylum) %>% 
  mutate_all(as.numeric) %>% 
  summarise(across(where(is.numeric), sum)) %>% 
  column_to_rownames(.,"phylum") %>% 
  mutate_all(function(x) x * 100)


df4 <- rbind(df3,group_info[,-1]) %>% 
  t() %>% 
  as.data.frame() %>% 
  group_by(seg) %>% 
  mutate_all(as.numeric) %>% 
  summarise(across(where(is.numeric), mean)) %>% 
  column_to_rownames(var = "seg") %>% 
  t() %>% 
  as.data.frame() %>% 
  arrange(rowSums(.)) %>% 
  rownames_to_column() %>% 
  pivot_longer(,cols = -1)

df4$rowname <- factor(df4$rowname, levels = unique(df4$rowname))
df4$name <- factor(df4$name, levels = c("Crop","Gizzard","Duodenum","Jejunum","Ileum","Cecum",    "Colon"))
col1 <- my_color[1:length(unique(df4$rowname))]

p1 <-  ggplot(df4, aes(name,value)) +
  geom_col(aes(fill = rowname))+
  #scale_color_manual(values = my_color) +
  scale_fill_manual(values = rev(col1)) +
  labs(x="",y="Relative Abundance (%)",fill="Phylum") +
  scale_x_discrete(expand = c(0.01, 0)) +  # default 0.2
  scale_y_continuous(expand = c(0.005, 0)) +  # default 0.05
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 8,color = "black",angle = 90,hjust = 1,vjust = 0.5),
    axis.text.y = element_text(size = 8,color = "black"),
    #axis.title = element_text(size = 8,color = "black"),
    legend.text = element_text(size =7, color = "black"), 
    legend.title = element_text(size =7, color = "black"), 
    panel.grid = element_blank(),
    legend.key.size = unit(0.3, "cm"),
    panel.border = element_blank()
  )

p1
library(ggh4x)
p1 <- p1 + force_panelsizes(total_width=unit(1.5, "in"),total_height=unit(3, "in"))
p1 
ggsave("phylum.pdf",p1)


## family -------------
df2 <- df %>% 
  select(c(3,7:ncol(.)))
group_info <- df2[1,]


df3 <- df2 %>% 
  slice(-1) %>% 
  group_by(family) %>% 
  mutate_all(as.numeric) %>% 
  summarise(across(where(is.numeric), sum)) %>% 
  column_to_rownames(.,"family") %>% 
  mutate_all(function(x) x * 100) 

df4 <- rbind(df3,group_info[,-1]) %>% 
  t() %>% 
  as.data.frame() %>% 
  group_by(seg) %>% 
  mutate_all(as.numeric) %>% 
  summarise(across(where(is.numeric), mean)) %>% 
  column_to_rownames(var = "seg") %>% 
  t() %>% 
  as.data.frame() %>% 
  arrange(desc(rowSums(.))) %>% 
  slice(c(1:4,6:11))

df4 <- rbind(df4, 100 - colSums(df4))
rownames(df4)[11] <- 'Others' 

df4 <- df4 %>% 
  rownames_to_column() %>% 
  mutate(., font_name=glue("<i>{rowname}")) %>% 
  select(.,c(1,ncol(.)),everything()) %>% 
  pivot_longer(.,cols = c(-1,-2))

df4$font_name <- gsub("<i>Others","Others",df4$font_name)

df4$rowname <- factor(df4$rowname, levels = rev(unique(df4$rowname)))

df4$name <- factor(df4$name, levels = c("Crop","Gizzard","Duodenum","Jejunum","Ileum","Cecum",    "Colon"))
col1 <- c(my_color[1:10],"gray")


p1 <-  ggplot(df4, aes(name,value)) +
  geom_col(aes(fill = rowname))+
  #scale_color_manual(values = my_color) +
  scale_fill_manual(values = rev(col1),labels = rev(unique(df4$font_name))) +
  labs(x="",y="Relative Abundance (%)",fill="Family") +
  scale_x_discrete(expand = c(0.01, 0)) +  # default 0.2
  scale_y_continuous(expand = c(0.005, 0)) +  # default 0.05
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 8,color = "black",angle = 90,hjust = 1,vjust = 0.5),
    axis.text.y = element_text(size = 8,color = "black"),
    axis.title = element_text(size = 8,color = "black"),
    legend.text = element_markdown(size =7, color = "black"), 
    legend.title = element_text(size = 7,color = "black"),
    panel.grid = element_blank(),
    legend.key.size = unit(0.3, "cm"),
    panel.border = element_blank()
  )

p1
library(ggh4x)
p1 <- p1 + force_panelsizes(total_width=unit(1.5, "in"),total_height=unit(3, "in"))
p1 
ggsave("family.pdf",p1)


## genus -------------
df2 <- df %>% 
  select(c(2,7:ncol(.)))
group_info <- df2[1,]


df3 <- df2 %>% 
  slice(-1) %>% 
  group_by(genus) %>% 
  mutate_all(as.numeric) %>% 
  summarise(across(where(is.numeric), sum)) %>% 
  column_to_rownames(.,"genus") %>% 
  mutate_all(function(x) x * 100) 


df4 <- rbind(df3,group_info[,-1]) %>% 
  t() %>% 
  as.data.frame() %>% 
  group_by(seg) %>% 
  mutate_all(as.numeric) %>% 
  summarise(across(where(is.numeric), mean)) %>% 
  column_to_rownames(var = "seg") %>% 
  t() %>% 
  as.data.frame() %>% 
  arrange(desc(rowSums(.))) %>% 
  slice(1:20)

df4 <- rbind(df4, 100 - colSums(df4))
rownames(df4)[21] <- 'Others' 

df4 <- df4 %>% 
  rownames_to_column() %>% 
  mutate(., font_name=glue("<i>{rowname}")) %>% 
  select(.,c(1,ncol(.)),everything()) %>% 
  pivot_longer(.,cols = c(-1,-2))


df4$font_name <- gsub("<i>Others","Others",df4$font_name)
df4$rowname <- factor(df4$rowname, levels = rev(unique(df4$rowname)))
df4$name <- factor(df4$name, levels = c("Crop","Gizzard","Duodenum","Jejunum","Ileum","Cecum",    "Colon"))
col1 <- c(my_color[1:20],"gray")


p1 <-  ggplot(df4, aes(name,value)) +
  geom_col(aes(fill = rowname))+
  #scale_color_manual(values = my_color) +
  scale_fill_manual(values = rev(col1),labels = rev(unique(df4$font_name))) +
  labs(x="",y="Relative Abundance (%)",fill="Genus") +
  scale_x_discrete(expand = c(0.01, 0)) +  # default 0.2
  scale_y_continuous(expand = c(0.005, 0)) +  # default 0.05
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 8,color = "black",angle = 90,hjust = 1,vjust = 0.5),
    axis.text.y = element_text(size = 8,color = "black"),
    axis.title = element_text(size = 8,color = "black"),
    legend.text = element_markdown(size =7, color = "black"), 
    legend.title = element_text(size = 7,color = "black"),
    panel.grid = element_blank(),
    legend.key.size = unit(0.3, "cm"),
    panel.border = element_blank()
  )+
  guides(fill = guide_legend(ncol = 1))

p1
library(ggh4x)
p1 <- p1 + force_panelsizes(total_width=unit(1.5, "in"),total_height=unit(3, "in"))
p1 
ggsave("genus.pdf",p1)


## species -------------
df2 <- df %>% 
  select(c(1,7:ncol(.)))
group_info <- df2[1,]


df3 <- df2 %>% 
  slice(-1) %>% 
  group_by(species) %>% 
  mutate_all(as.numeric) %>% 
  summarise(across(where(is.numeric), sum)) %>% 
  column_to_rownames(.,"species") %>% 
  mutate_all(function(x) x * 100) 


df4 <- rbind(df3,group_info[,-1]) %>% 
  t() %>% 
  as.data.frame() %>% 
  group_by(seg) %>% 
  mutate_all(as.numeric) %>% 
  summarise(across(where(is.numeric), mean)) %>% 
  column_to_rownames(var = "seg") %>% 
  t() %>% 
  as.data.frame() %>% 
  arrange(desc(rowSums(.))) %>% 
  slice(1:25)

df4 <- rbind(df4, 100 - colSums(df4))
rownames(df4)[26] <- 'Others' 

df4 <- df4 %>% 
  rownames_to_column() %>% 
  mutate(., font_name=glue("<i>{rowname}")) %>% 
  select(.,c(1,ncol(.)),everything()) %>% 
  pivot_longer(.,cols = c(-1,-2))


df4$font_name <- gsub("<i>Others","Others",df4$font_name)
df4$rowname <- factor(df4$rowname, levels = rev(unique(df4$rowname)))
df4$name <- factor(df4$name, levels = c("Crop","Gizzard","Duodenum","Jejunum","Ileum","Cecum",    "Colon"))
col1 <- c(my_color[1:25],"gray")


p1 <-  ggplot(df4, aes(name,value)) +
  geom_col(aes(fill = rowname))+
  #scale_color_manual(values = my_color) +
  scale_fill_manual(values = rev(col1),labels = rev(unique(df4$font_name))) +
  labs(x="",y="Relative Abundance (%)",fill="Species") +
  scale_x_discrete(expand = c(0.01, 0)) +  # default 0.2
  scale_y_continuous(expand = c(0.005, 0)) +  # default 0.05
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 8,color = "black",angle = 90,hjust = 1,vjust = 0.5),
    axis.text.y = element_text(size = 8,color = "black"),
    axis.title = element_text(size = 8,color = "black"),
    legend.text = element_markdown(size =7, color = "black"), 
    legend.title = element_text(size = 7,color = "black"),
    panel.grid = element_blank(),
    legend.key.size = unit(0.3, "cm"),
    panel.border = element_blank()
  )+
  guides(fill = guide_legend(ncol = 1))

p1
library(ggh4x)
p1 <- p1 + force_panelsizes(total_width=unit(1.5, "in"),total_height=unit(3, "in"))
p1 
ggsave("species.pdf",p1)

# positive control -----------
setwd("D:/OneDrive - Oklahoma A and M System/tranfer documents/16s_fulength/emu_new")
df <- read.csv("feature_table.csv",header=T,row.names=1) %>% 
  t() %>% 
  as.data.frame() 

df_taxon <- df[1:2,]

df_data <- df %>% 
  filter(seg == "Pos")

df_final <- rbind(df_taxon,df_data) %>% 
  remove_rownames() %>% 
  column_to_rownames(.,var = 'id') %>% 
  select(-1) %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(barcode83 = as.numeric(barcode83),
        barcode84 = as.numeric(barcode84))

df_species <- df_final %>% 
  select(-2)
