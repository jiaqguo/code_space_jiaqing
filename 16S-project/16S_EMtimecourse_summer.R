#######################################################
####       select until LINE 1754 and run         #####          
####    other codes are for you for exploration   #####
#######################################################

#setwd("D:/OneDrive - Oklahoma A and M System/tranfer documents/workshop/results")

setwd("C:/Users/jiaguo/OneDrive - Oklahoma A and M System/Bioinformatics/work shop/workshop-AMBL/results")


# modify tsv file as csv file -----------
library(tidyverse)
# for feature table
feature_table <- read.table("feature_table.tsv",row.names = 1,sep = "\t",comment.char = "$",skip = 1)
colnames(feature_table) <- feature_table[1,]
feature_table <- feature_table[-1,]
write.csv(feature_table,"feature_table.csv")

# for tax table
tax <- read.table("taxonomy.tsv",row.names = 1,header = T,sep = "\t") %>% 
  select(-2) %>% 
  separate(Taxon,into = c("Domain","Phylum","Class","Order","Family","Genus","Species"),sep = ";")
tax[is.na(tax)] <- ""

# remove "incertae_sedis"
tax <- tax %>%
  mutate(across(where(is.character), 
                ~ gsub("incertae_sedis", "Unclassified", .)))


# add unclassified information
for (row_num in 1:nrow(tax)){
    if (tax$Domain[row_num] == ""){
      tax$Domain[row_num] = tax$Phylum[row_num] = tax$Class[row_num] = tax$Order[row_num] = tax$Family[row_num] = tax$Genus[row_num] = tax$Species[row_num] = "Unclassified"
    } else if (tax$Phylum[row_num] == ""){
      tax$Phylum[row_num] = tax$Class[row_num] = tax$Order[row_num] = tax$Family[row_num] = tax$Genus[row_num] = tax$Species[row_num] = paste0(tax$Domain[row_num],"_Unclassified")
    } else if (tax$Class[row_num] == ""){
      tax$Class[row_num] = tax$Order[row_num] = tax$Family[row_num] = tax$Genus[row_num] = tax$Species[row_num] = paste0(tax$Phylum[row_num],"_Unclassified")
    } else if (tax$Order[row_num] == ""){
      tax$Order[row_num] = tax$Family[row_num] = tax$Genus[row_num] = tax$Species[row_num] = paste0(tax$Class[row_num],"_Unclassified")
    } else if (tax$Family[row_num] == ""){
      tax$Family[row_num] = tax$Genus[row_num] = tax$Species[row_num] = paste0(tax$Order[row_num],"_Unclassified")
    } else if (tax$Genus[row_num] == ""){
      tax$Genus[row_num] = tax$Species[row_num] = paste0(tax$Family[row_num],"_Unclassified")
    }
      else if (tax$Species[row_num] == ""){
      tax$Species[row_num] = paste0(tax$Genus[row_num],"_Unclassified")
    }
}

# remove dobule "Unclassified"
tax <- tax %>%
  mutate(across(where(is.character), 
                ~ gsub("Unclassified_.*", "Unclassified", .)))


# add F1,F2,F3.... information
feature_table <- feature_table %>% 
  mutate(across(everything(),as.numeric)) %>% 
  mutate(across(everything(),~./sum(.))) %>% 
  arrange(desc(rowSums(.)))

tax <- tax[rownames(feature_table),] %>% 
  mutate(feature_num = paste0("F",seq(1:nrow(.)))) %>% 
  unite("Species",c(Species,feature_num),sep = "_")

write.csv(tax,"taxonomy.csv")

# filter feature table 
feature_table_filter <- read.csv("feature_table.csv",header = T,row.names = 1)
feature_table_filter <- feature_table_filter[rownames(feature_table),]
feature_table_filter <- feature_table_filter[rowSums(feature_table_filter>0) >= ncol(feature_table_filter)*0.05,]
write.csv(feature_table_filter,"feature_table_filter.csv")

# Diversity------------
f_table <- "feature_table_filter.csv"
s_table <- "metadata.csv"
t_table <- "taxonomy.csv"

library(phyloseq)
library(vegan)
library(ape)
feature_table <- read.csv(f_table,header = T,row.names = 1)
taxonomy <- read.csv(t_table,header = T,row.names = 1)
sample_data <- read.csv(s_table,header = T,row.names = 1)

#Transform feature table and taxonomy table to matrix
OTU<- otu_table(as.matrix(feature_table),taxa_are_rows = T)
TAX<- tax_table(as.matrix(taxonomy))
SAM<-sample_data(sample_data)
TREE<-read_tree('tree.nwk')
TREE <- multi2di(TREE)
physeq<- phyloseq(OTU,TAX,SAM,TREE)
#rarefy data
rarecurve(t(feature_table), step=200,lty = 2)
physeq<-rarefy_even_depth(physeq,replace = F,rngseed = T)
otu_rarefy <- as.data.frame(otu_table(physeq))
write.csv(otu_rarefy,file = "feature_table_rarefy.csv")
## α-diversity-----
# plot_richness(phy_3,x="treatment",measures = c('Observed','Shannon'))
dir.create("a_diversity")
phy_3 <- subset_samples(physeq,day=="3 dpi")
phy_5 <- subset_samples(physeq,day=="5 dpi")
phy_7 <- subset_samples(physeq,day=="7 dpi")
phy_14 <- subset_samples(physeq,day=="14 dpi")
phy_21 <- subset_samples(physeq,day=="21 dpi")

a_diversity <- estimate_richness(phy_3,measures = c('Observed','Shannon'))
pielou <- a_diversity$Shannon / log(a_diversity$Observed)
a_diversity <- cbind(a_diversity,pielou)
write.csv(a_diversity,file = "a_diversity/a_diversity_3.csv")

a_diversity <- estimate_richness(phy_5,measures = c('Observed','Shannon'))
pielou <- a_diversity$Shannon / log(a_diversity$Observed)
a_diversity <- cbind(a_diversity,pielou)
write.csv(a_diversity,file = "a_diversity/a_diversity_5.csv")

a_diversity <- estimate_richness(phy_7,measures = c('Observed','Shannon'))
pielou <- a_diversity$Shannon / log(a_diversity$Observed)
a_diversity <- cbind(a_diversity,pielou)
write.csv(a_diversity,file = "a_diversity/a_diversity_7.csv")

a_diversity <- estimate_richness(phy_14,measures = c('Observed','Shannon'))
pielou <- a_diversity$Shannon / log(a_diversity$Observed)
a_diversity <- cbind(a_diversity,pielou)
write.csv(a_diversity,file = "a_diversity/a_diversity_14.csv")

a_diversity <- estimate_richness(phy_21,measures = c('Observed','Shannon'))
pielou <- a_diversity$Shannon / log(a_diversity$Observed)
a_diversity <- cbind(a_diversity,pielou)
write.csv(a_diversity,file = "a_diversity/a_diversity_21.csv")

## β-diversity----
library(tidyverse)
unlist(distanceMethodList)  #check the methods to calculate distance matrix
dir.create("b-diversity")
#PCoA analysis  
for (i in c(3,5,7,14,21)){
  if (i == 3){
    name = phy_3
  } else if (i == 5){
    name = phy_5
  } else if (i == 7){
    name = phy_7
  } else if (i == 14){
    name = phy_14
  } else if (i == 21){
    name = phy_21
  }
  e_name=paste("phy",i,sep = "_")
  unifrac_ord <- ordinate(name, method = 'PCoA',distance = 'unifrac')
  wunifrac_ord <- ordinate(name, method = 'PCoA',distance = 'wunifrac')

  #Export distance matrix, coordinates and explained variance of principal components
  unifrac_info <- as.data.frame(cbind(rownames(unifrac_ord$vectors),
                                      unifrac_ord$vectors[,1],
                                      unifrac_ord$vectors[,2],
                                      unifrac_ord$values[,2]))
  colnames(unifrac_info) <- c("SampleID","axis1","axis2","explanation")
  sample_data$SampleID <- rownames(sample_data)
  unifrac_info <- left_join(unifrac_info,sample_data,by="SampleID")
  
  wunifrac_info <- cbind(as.data.frame(rownames(wunifrac_ord$vectors)),
                         as.data.frame(wunifrac_ord$vectors)$Axis.1,
                         as.data.frame(wunifrac_ord$vectors)$Axis.2,
                         as.data.frame(wunifrac_ord$values)$Relative_eig)
  colnames(wunifrac_info) <- c("SampleID","axis1","axis2","explanation")
  wunifrac_info <- left_join(wunifrac_info,sample_data,by="SampleID")
  
  write.csv(unifrac_info,file = paste0("b-diversity/unifrac_info_",e_name,".csv"))
  write.csv(wunifrac_info,file = paste0("b-diversity/wunifrac_info_",e_name,".csv"))
}

## statistical analysis for β-diversity-----------
library(vegan)

for (i in c(3,5,7,14,21)){
  if (i == 3){
    name = phy_3
  } else if (i == 5){
    name = phy_5
  } else if (i == 7){
    name = phy_7
  } else if (i == 14){
    name = phy_14
  } else if (i == 21){
    name = phy_21
  }
  e_name=paste("phy",i,sep = "_")
  dpi = paste(i,"dpi",sep = " ")
  dis_unifrac <- distance(name, method ='unifrac')
  dis_wunifrac <- distance(name, method ='wunifrac')
  
  order <- rownames(sample_data[sample_data$day==dpi,])
  
  dis_unifrac <- as.matrix(dis_unifrac)
  dis_wunifrac <- as.matrix(dis_wunifrac)
  
  dis_unifrac <- dis_unifrac[order,order]
  dis_wunifrac <- dis_wunifrac[order,order]
  
  identical(rownames(as.matrix(dis_unifrac)), rownames(sample_data[sample_data$day==dpi,]))
  identical(rownames(as.matrix(dis_wunifrac)), rownames(sample_data[sample_data$day==dpi,]))
  identical(colnames(as.matrix(dis_unifrac)), rownames(sample_data[sample_data$day==dpi,]))
  identical(colnames(as.matrix(dis_wunifrac)), rownames(sample_data[sample_data$day==dpi,]))
  
  set.seed(240923)
  unifrac_premanova <- adonis2(dis_unifrac ~ treatment,data=sample_data[sample_data$day==dpi,],permutations = 999)
  set.seed(240923)
  wunifrac_premanova <- adonis2(dis_wunifrac ~ treatment,data=sample_data[sample_data$day==dpi,],permutations = 999)
  unifrac_premanova
  wunifrac_premanova
  # add source
  unifrac_premanova$source <- paste0('unifrac_premanova_',e_name)
  wunifrac_premanova$source <- paste0('wunifrac_premanova',e_name)
  
  result <- rbind(unifrac_premanova,NA,wunifrac_premanova)
  
  write.csv(result,file = paste0("b-diversity/premanova_",e_name,".csv"))
}



# Calculate relative abundance-----
feature_table <-read.csv("feature_table_rarefy.csv",header = T,row.names = 1)

rela_abun <- apply(feature_table,2,FUN = function(x){
  x/sum(x)
})

sum(round(colSums(rela_abun))!=1)  #oupt 0

colnames(rela_abun) <- colnames(feature_table)
rownames(rela_abun) <- rownames(feature_table)
write.csv(rela_abun,file = "feature_table_relative.csv")
#Over 

# Taxonomy composition----------
dir.create("composition")
f_table <- "feature_table_relative.csv"
s_table <- "metadata.csv"
t_table <- "taxonomy.csv"
otu <-read.csv(f_table,header = T,row.names = 1)
sample_data <-read.csv(s_table,header = T,row.names = 1)
tax <- read.csv(t_table,header = T,row.names = 1)
otu$rownames <- rownames(otu)
tax$rownames <- rownames(tax)

#TOP Phylum
phylum_otu <- merge(tax, otu, by='rownames')
rownames(phylum_otu) <- phylum_otu$rownames
phylum_otu <- phylum_otu[,c(-1,-2,-4,-5,-6,-7,-8,-9)]
phylum_otu <- aggregate(phylum_otu[,-1], by=list(phylum_otu$Phylum), FUN='sum')
rownames(phylum_otu) <- phylum_otu$Group.1
phylum_otu <- phylum_otu[,-1]
phylum_otu <- phylum_otu[order(rowSums(phylum_otu),decreasing = T),]
write.csv(phylum_otu,file = 'composition/phylum.csv')


#Top 10 Family
number=10
family_otu <- merge(tax, otu, by='rownames')
rownames(family_otu) <- family_otu$rownames
family_otu <- family_otu[,c(-1,-2,-3,-4,-5,-7,-8,-9)]
family_otu <- aggregate(family_otu[,-1], by=list(family_otu$Family), FUN='sum')
rownames(family_otu) <- family_otu$Group.1
family_otu <- family_otu[,-1]
family_otu <- family_otu[order(rowSums(family_otu),decreasing = T),]
family_otu <- family_otu[c(1:number),]   # top 10
family_otu <- rbind(family_otu, 1 - colSums(family_otu))
rownames(family_otu)[number+1] <- 'Others'  # 10 + 1 is others
write.csv(family_otu,file = 'composition/family_top10.csv')


#Top 20 Genus
number=20
genus_otu <- merge(tax, otu, by='rownames')
rownames(genus_otu) <- genus_otu$rownames
genus_otu <- genus_otu[,c(-1,-2,-3,-4,-5,-6,-8,-9)]
genus_otu <- aggregate(genus_otu[,-1], by=list(genus_otu$Genus), FUN='sum')
rownames(genus_otu) <- genus_otu$Group.1
genus_otu <- genus_otu[,-1]
genus_otu <- genus_otu[order(rowSums(genus_otu),decreasing = T),]
genus_otu <- genus_otu[c(1:number),]
genus_otu <- rbind(genus_otu, 1 - colSums(genus_otu))
rownames(genus_otu)[number+1] <- 'Others'
write.csv(genus_otu,file = 'composition/genus_top20.csv')


#TOP number (25/30) ASV
number=25
feature_otu <- merge(tax, otu, by='rownames')
rownames(feature_otu) <-feature_otu$Species
feature_otu <- feature_otu[,-c(1:9)]
feature_otu <- feature_otu[order(rowSums(feature_otu),decreasing = T),]
feature_otu <- feature_otu[1:number,]
feature_otu <- rbind(feature_otu, 1 - colSums(feature_otu))
rownames(feature_otu)[number+1] <- "Others"
write.csv(feature_otu,file = 'composition/feature_top25.csv')


# Visualization--------
library(tidyverse)
library(car)
# a-diversity --------
##statiscal analysis----
# Normality test and Homogeneity of variance test
check_function <- function(data_points,group_var,data){
  # construct the formula from the column name strings
  formula_str <- paste(data_points, "~", group_var)
  formula_obj <- as.formula(formula_str)
  # Fit the ANOVA model using the constructed formula
  fit <- aov(formula_obj, data = data)
  res <- residuals(fit)
  # Get the vectors needed for plotting (safe extraction using [[ operator]])
  group_data <- data[[group_var]]
  par(mfrow=c(1,2))
  # Normality Check
  qqnorm(res, main = "Normal Q-Q Plot of Residuals"); 
  qqline(res)
  cat("\nNormality test (Shapiro-Wilk):\n")
  print(shapiro.test(res))
  # Homogeneity of Variance Check (Residuals vs. Fitted/Group)
  cat("\nResiduals vs. Group Plot (Homogeneity):\n")
  plot(group_data, res, 
       xlab = group_var, 
       ylab = "Residuals",
       main = "Residuals vs. Group")
  abline(h=0, col="red")
  cat("\nHomogeneity of variance test (Levene's Test):\n")
  # Use the constructed formula with leveneTest
  print(leveneTest(formula_obj, data = data))
  par(mfrow=c(1,1))
}

read_file_function <- function(index,data,sample_data){
  if (index == "observed"){
    num = 2
  } else if (index == "shannon"){
    num = 3
  } else if (index == "pielou"){
    num = 4
  }
  a_diversity <- read.csv(data,header = T,row.names = 1)
  sample_data <- read.csv(sample_data,header = T,row.names = 1)
  a_diversity$SampleID <- rownames(a_diversity)
  sample_data$SampleID <- rownames(sample_data)
  a_pic <- merge(a_diversity,sample_data,by='SampleID')
  table(a_pic$treatment)    #check the different treatments of group
  a_pic$treatment <- as.factor(a_pic$treatment)
  df <- a_pic[,c(num,5,6)]
  return(df)
}

###observed ASVs--------------------

# 3 dpi
d_file <- "a_diversity/a_diversity_3.csv"
s_file <- "metadata.csv"
index <- "observed"
# read data
df <- read_file_function(index,d_file,s_file)
# analysis
check_function("Observed","treatment",df)
# statistical analysis
t_results <- t.test(Observed~treatment, data = df)
w_results <- wilcox.test(Observed ~treatment, data = df)
# gather results
stat_results <- as.data.frame(matrix(ncol = 3,nrow = 15))
colnames(stat_results) <- c("index","dpi","p.value")
stat_results$index[1] <- index
stat_results$dpi[1] <- unique(df$day)
stat_results$p.value[1] <- w_results$p.value

# 5 dpi
d_file <- "a_diversity/a_diversity_5.csv"
s_file <- "metadata.csv"
index <- "observed"
# read data
df <- read_file_function(index,d_file,s_file)
# analysis
check_function("Observed","treatment",df)
# statistical analysis
t_results <- t.test(Observed~treatment, data = df)
w_results <- wilcox.test(Observed ~treatment, data = df)
# gather results
stat_results$index[2] <- index
stat_results$dpi[2] <- unique(df$day)
stat_results$p.value[2] <- t_results$p.value

# 7 dpi
d_file <- "a_diversity/a_diversity_7.csv"
s_file <- "metadata.csv"
index <- "observed"
# read data
df <- read_file_function(index,d_file,s_file)
# analysis
check_function("Observed","treatment",df)
# statistical analysis
t_results <- t.test(Observed~treatment, data = df)
w_results <- wilcox.test(Observed ~treatment, data = df)
# gather results
stat_results$index[3] <- index
stat_results$dpi[3] <- unique(df$day)
stat_results$p.value[3] <- t_results$p.value

# 14 dpi
d_file <- "a_diversity/a_diversity_14.csv"
s_file <- "metadata.csv"
index <- "observed"
# read data
df <- read_file_function(index,d_file,s_file)
# analysis
check_function("Observed","treatment",df)
# statistical analysis
t_results <- t.test(Observed~treatment, data = df)
w_results <- wilcox.test(Observed ~treatment, data = df)
# gather results
stat_results$index[4] <- index
stat_results$dpi[4] <- unique(df$day)
stat_results$p.value[4] <- w_results$p.value

# 21 dpi
d_file <- "a_diversity/a_diversity_21.csv"
s_file <- "metadata.csv"
index <- "observed"
# read data
df <- read_file_function(index,d_file,s_file)
# analysis
check_function("Observed","treatment",df)
# statistical analysis
t_results <- t.test(Observed~treatment, data = df)
w_results <- wilcox.test(Observed ~treatment, data = df)
# gather results
stat_results$index[5] <- index
stat_results$dpi[5] <- unique(df$day)
stat_results$p.value[5] <- w_results$p.value


### Shannon--------------------
# 3 dpi
d_file <- "a_diversity/a_diversity_3.csv"
s_file <- "metadata.csv"
index <- "shannon"
# read data
df <- read_file_function(index,d_file,s_file)
# analysis
check_function("Shannon","treatment",df)
# statistical analysis
t_results <- t.test(Shannon~treatment, data = df)
w_results <- wilcox.test(Shannon ~treatment, data = df)
# gather results
stat_results$index[6] <- index
stat_results$dpi[6] <- unique(df$day)
stat_results$p.value[6] <- t_results$p.value

# 5 dpi
d_file <- "a_diversity/a_diversity_5.csv"
s_file <- "metadata.csv"
index <- "shannon"
# read data
df <- read_file_function(index,d_file,s_file)
# analysis
check_function("Shannon","treatment",df)
# statistical analysis
t_results <- t.test(Shannon~treatment, data = df)
w_results <- wilcox.test(Shannon ~treatment, data = df)
# gather results
stat_results$index[7] <- index
stat_results$dpi[7] <- unique(df$day)
stat_results$p.value[7] <- w_results$p.value

# 7 dpi
d_file <- "a_diversity/a_diversity_7.csv"
s_file <- "metadata.csv"
index <- "shannon"
# read data
df <- read_file_function(index,d_file,s_file)
# analysis
check_function("Shannon","treatment",df)
# statistical analysis
t_results <- t.test(Shannon~treatment, data = df)
w_results <- wilcox.test(Shannon ~treatment, data = df)
# gather results
stat_results$index[8] <- index
stat_results$dpi[8] <- unique(df$day)
stat_results$p.value[8] <- w_results$p.value

# 14 dpi
d_file <- "a_diversity/a_diversity_14.csv"
s_file <- "metadata.csv"
index <- "shannon"
# read data
df <- read_file_function(index,d_file,s_file)
# analysis
check_function("Shannon","treatment",df)
# statistical analysis
t_results <- t.test(Shannon~treatment, data = df)
w_results <- wilcox.test(Shannon ~treatment, data = df)
# gather results
stat_results$index[9] <- index
stat_results$dpi[9] <- unique(df$day)
stat_results$p.value[9] <- w_results$p.value

# 21 dpi
d_file <- "a_diversity/a_diversity_21.csv"
s_file <- "metadata.csv"
index <- "shannon"
# read data
df <- read_file_function(index,d_file,s_file)
# analysis
check_function("Shannon","treatment",df)
# statistical analysis
t_results <- t.test(Shannon~treatment, data = df)
w_results <- wilcox.test(Shannon ~treatment, data = df)
# gather results
stat_results$index[10] <- index
stat_results$dpi[10] <- unique(df$day)
stat_results$p.value[10] <- w_results$p.value

### pielou--------------------
# 3 dpi
d_file <- "a_diversity/a_diversity_3.csv"
s_file <- "metadata.csv"
index <- "pielou"
# read data
df <- read_file_function(index,d_file,s_file)
# analysis
check_function("pielou","treatment",df)
# statistical analysis
t_results <- t.test(pielou~treatment, data = df)
w_results <- wilcox.test(pielou ~treatment, data = df)
# gather results
stat_results$index[11] <- index
stat_results$dpi[11] <- unique(df$day)
stat_results$p.value[11] <- t_results$p.value

# 5 dpi
d_file <- "a_diversity/a_diversity_5.csv"
s_file <- "metadata.csv"
index <- "pielou"
# read data
df <- read_file_function(index,d_file,s_file)
# analysis
check_function("pielou","treatment",df)
# statistical analysis
t_results <- t.test(pielou~treatment, data = df)
w_results <- wilcox.test(pielou ~treatment, data = df)
# gather results
stat_results$index[12] <- index
stat_results$dpi[12] <- unique(df$day)
stat_results$p.value[12] <- w_results$p.value

# 7 dpi
d_file <- "a_diversity/a_diversity_7.csv"
s_file <- "metadata.csv"
index <- "pielou"
# read data
df <- read_file_function(index,d_file,s_file)
# analysis
check_function("pielou","treatment",df)
# statistical analysis
t_results <- t.test(pielou~treatment, data = df)
w_results <- wilcox.test(pielou ~treatment, data = df)
# gather results
stat_results$index[13] <- index
stat_results$dpi[13] <- unique(df$day)
stat_results$p.value[13] <- w_results$p.value

# 14 dpi
d_file <- "a_diversity/a_diversity_14.csv"
s_file <- "metadata.csv"
index <- "pielou"
# read data
df <- read_file_function(index,d_file,s_file)
# analysis
check_function("pielou","treatment",df)
# statistical analysis
t_results <- t.test(pielou~treatment, data = df)
w_results <- wilcox.test(pielou ~treatment, data = df)
# gather results
stat_results$index[14] <- index
stat_results$dpi[14] <- unique(df$day)
stat_results$p.value[14] <- t_results$p.value

# 21 dpi
d_file <- "a_diversity/a_diversity_21.csv"
s_file <- "metadata.csv"
index <- "pielou"
# read data
df <- read_file_function(index,d_file,s_file)
# analysis
check_function("pielou","treatment",df)
# statistical analysis
t_results <- t.test(pielou~treatment, data = df)
w_results <- wilcox.test(pielou ~treatment, data = df)
# gather results
stat_results$index[15] <- index
stat_results$dpi[15] <- unique(df$day)
stat_results$p.value[15] <- w_results$p.value

# export 
write.csv(stat_results,file = "a_diversity/a_diversity_stat.csv")



###fig------
level_order <- c("Mock","EM")
a_diversity1 <- read.csv("a_diversity/a_diversity_3.csv",header = T,row.names = 1)
a_diversity2 <- read.csv("a_diversity/a_diversity_5.csv",header = T,row.names = 1)
a_diversity3 <- read.csv("a_diversity/a_diversity_7.csv",header = T,row.names = 1)
a_diversity4 <- read.csv("a_diversity/a_diversity_14.csv",header = T,row.names = 1)
a_diversity5 <- read.csv("a_diversity/a_diversity_21.csv",header = T,row.names = 1)
sample_data <- read.csv("metadata.csv",header = T,row.names = 1)
a_diversity <- rbind(a_diversity1,a_diversity2,a_diversity3,a_diversity4,a_diversity5)
a_diversity$SampleID <- rownames(a_diversity)
sample_data$SampleID <- rownames(sample_data)
a_pic <- merge(a_diversity,sample_data,by='SampleID')
table(a_pic$treatment)    #check the different treatments of group
a_pic$treatment <- factor(a_pic$treatment,levels = level_order)
a_pic$day <- factor(a_pic$day,levels = c("3 dpi","5 dpi","7 dpi","14 dpi","21 dpi"))

my_color_color <- c(Mock="#377483", EM="#B83945")
my_color_fill <- c(Mock="#C7DFF0", EM="#FBDFE2")

# observed
annotation_data <- data.frame(
  day = c("3 dpi","5 dpi","7 dpi","14 dpi","21 dpi"), 
  x = c(1.5,2,1.5,1.5,2),      
  y = c(190,170,190,100,200),       
  label = c(paste0("italic(P) == ", 0.598),
            "***",
            paste0("italic(P) == ", 0.182),
            paste0("italic(P) == ", 0.921),
            "*")
)
annotation_data$day <- factor(annotation_data$day,levels = c("3 dpi","5 dpi","7 dpi","14 dpi","21 dpi"))
# 1. Data that needs plotmath parsing (P-values)
p_value_data <- annotation_data %>% filter(grepl("italic", label)) 
# 2. Data that should be treated as plain text (Stars)
star_data <- annotation_data %>% filter(!grepl("italic", label)) 


#1pt = 0.3527 mm. geom_text's unit is mm
p1 <- ggplot(a_pic,aes(x=treatment,y=Observed))+
  facet_grid(.~day)+
  geom_boxplot(aes(fill=treatment,color=treatment),width=0.6)+
  #scale_y_continuous(limits = c(0,265),breaks = seq(0,250,50))+
  labs(x="",y="",title = "Observed ASVs")+
  scale_fill_manual(values = my_color_fill)+
  scale_color_manual(values = my_color_color)+
  #geom_line(data = med_rep,aes(x=x,y=medvalue),color='white',linewidth=0.6)+
  theme_bw()+
  theme(panel.grid = element_blank(),legend.position = "none",
        axis.text.x = element_text(size = 8,vjust = 1,hjust = 1,
                                   color = "black",angle = 45,
                                   margin = margin(t = 0)),
        axis.text.y = element_text(size = 8,color = "black"),
        panel.border = element_rect(linewidth = 0.5,color = 'black'),
        axis.ticks.x = element_blank(), 
        #axis.ticks.x = element_line(linewidth = 0.3), 
        axis.ticks.y = element_line(linewidth  = 0.3),
        plot.title = element_text(hjust = 0.5,size=9,margin = margin(b = 3)),
        strip.text = element_text(size = 8, color = "black"), #adjust facet_grid text
        panel.spacing.x= unit(0.1, "lines") #adjust facet_grid horizontal distance
  ) +
  geom_text(
    data = p_value_data,
    aes(x = x, y = y, label = label), 
    size = 3, 
    inherit.aes = FALSE, 
    parse = TRUE 
  ) +
  # --- Layer 2: SIGNIFICANCE STARS (No Parsing) ---
  geom_text(
    data = star_data,
    aes(x = x, y = y, label = label), 
    size = 3, 
    inherit.aes = FALSE,
    parse = FALSE 
  )


library(ggh4x)
p1 <- p1 + force_panelsizes(total_width=unit(3.3, "in"),total_height=unit(1.0, "in"))
p1 
ggsave("a_diversity/observed.pdf",p1,width = 3.3*1.2,height = 1*2)
ggsave("a_diversity/observed.png",p1,dpi = 600,width = 3.3*1.2,height = 1*2)


# shannon
annotation_data <- data.frame(
  day = c("3 dpi","5 dpi","7 dpi","14 dpi","21 dpi"), 
  x = c(1.5,2,2,1.5,1.5),      
  y = c(2.5,3.8,3.8,2.5,2),       
  label = c(paste0("italic(P) == ", 0.636),
            "***",
            "*",
            paste0("italic(P) == ", 0.283),
            paste0("italic(P) == ", 0.105))
)
annotation_data$day <- factor(annotation_data$day,levels = c("3 dpi","5 dpi","7 dpi","14 dpi","21 dpi"))
# 1. Data that needs plotmath parsing (P-values)
p_value_data <- annotation_data %>% filter(grepl("italic", label)) 
# 2. Data that should be treated as plain text (Stars)
star_data <- annotation_data %>% filter(!grepl("italic", label)) 


#1pt = 0.3527 mm. geom_text's unit is mm
p1 <- ggplot(a_pic,aes(x=treatment,y=Shannon))+
  facet_grid(.~day)+
  geom_boxplot(aes(fill=treatment,color=treatment),width=0.6)+
  #scale_y_continuous(limits = c(0,265),breaks = seq(0,250,50))+
  labs(x="",y="",title = "Shannon Index")+
  scale_fill_manual(values = my_color_fill)+
  scale_color_manual(values = my_color_color)+
  #geom_line(data = med_rep,aes(x=x,y=medvalue),color='white',linewidth=0.6)+
  theme_bw()+
  theme(panel.grid = element_blank(),legend.position = "none",
        axis.text.x = element_text(size = 8,vjust = 1,hjust = 1,
                                   color = "black",angle = 45,
                                   margin = margin(t = 0)),
        axis.text.y = element_text(size = 8,color = "black"),
        panel.border = element_rect(linewidth = 0.5,color = 'black'),
        axis.ticks.x = element_blank(), 
        #axis.ticks.x = element_line(linewidth = 0.3), 
        axis.ticks.y = element_line(linewidth  = 0.3),
        plot.title = element_text(hjust = 0.5,size=9,margin = margin(b = 3)),
        strip.text = element_text(size = 8, color = "black"), #adjust facet_grid text
        panel.spacing.x= unit(0.1, "lines") #adjust facet_grid horizontal distance
  ) +
  geom_text(
    data = p_value_data,
    aes(x = x, y = y, label = label), 
    size = 3, 
    inherit.aes = FALSE, 
    parse = TRUE 
  ) +
  # --- Layer 2: SIGNIFICANCE STARS (No Parsing) ---
  geom_text(
    data = star_data,
    aes(x = x, y = y, label = label), 
    size = 3, 
    inherit.aes = FALSE,
    parse = FALSE 
  )


library(ggh4x)
p1 <- p1 + force_panelsizes(total_width=unit(3.3, "in"),total_height=unit(1.0, "in"))
p1 
ggsave("a_diversity/shannon.pdf",p1,width = 3.3*1.2,height = 1*2)
ggsave("a_diversity/shannon.png",p1,dpi = 600,width = 3.3*1.2,height = 1*2)


# pielou
annotation_data <- data.frame(
  day = c("3 dpi","5 dpi","7 dpi","14 dpi","21 dpi"), 
  x = c(1.5,2,1.5,1.5,1.5),      
  y = c(0.5,0.75,0.5,0.5,0.4),       
  label = c(paste0("italic(P) == ", 0.687),
            "***",
            paste0("italic(P) == ", 0.080),
            paste0("italic(P) == ", 0.150),
            paste0("italic(P) == ", 0.190))
)
annotation_data$day <- factor(annotation_data$day,levels = c("3 dpi","5 dpi","7 dpi","14 dpi","21 dpi"))
# 1. Data that needs plotmath parsing (P-values)
p_value_data <- annotation_data %>% filter(grepl("italic", label)) 
# 2. Data that should be treated as plain text (Stars)
star_data <- annotation_data %>% filter(!grepl("italic", label)) 


#1pt = 0.3527 mm. geom_text's unit is mm
p1 <- ggplot(a_pic,aes(x=treatment,y=pielou))+
  facet_grid(.~day)+
  geom_boxplot(aes(fill=treatment,color=treatment),width=0.6)+
  #scale_y_continuous(limits = c(0,265),breaks = seq(0,250,50))+
  labs(x="",y="",title = "Pielou's Evenness")+
  scale_fill_manual(values = my_color_fill)+
  scale_color_manual(values = my_color_color)+
  #geom_line(data = med_rep,aes(x=x,y=medvalue),color='white',linewidth=0.6)+
  theme_bw()+
  theme(panel.grid = element_blank(),legend.position = "none",
        axis.text.x = element_text(size = 8,vjust = 1,hjust = 1,
                                   color = "black",angle = 45,
                                   margin = margin(t = 0)),
        axis.text.y = element_text(size = 8,color = "black"),
        panel.border = element_rect(linewidth = 0.5,color = 'black'),
        axis.ticks.x = element_blank(), 
        #axis.ticks.x = element_line(linewidth = 0.3), 
        axis.ticks.y = element_line(linewidth  = 0.3),
        plot.title = element_text(hjust = 0.5,size=9,margin = margin(b = 3)),
        strip.text = element_text(size = 8, color = "black"), #adjust facet_grid text
        panel.spacing.x= unit(0.1, "lines") #adjust facet_grid horizontal distance
  ) +
  geom_text(
    data = p_value_data,
    aes(x = x, y = y, label = label), 
    size = 3, 
    inherit.aes = FALSE, 
    parse = TRUE 
  ) +
  # --- Layer 2: SIGNIFICANCE STARS (No Parsing) ---
  geom_text(
    data = star_data,
    aes(x = x, y = y, label = label), 
    size = 3, 
    inherit.aes = FALSE,
    parse = FALSE 
  )


library(ggh4x)
p1 <- p1 + force_panelsizes(total_width=unit(3.3, "in"),total_height=unit(1.0, "in"))
p1 
ggsave("a_diversity/pielou.pdf",p1,width = 3.3*1.2,height = 1*2)
ggsave("a_diversity/pielou.png",p1,dpi = 600,width = 3.3*1.2,height = 1*2)


# β diversity------------
#library(ggh4x)   #stat_centroid
#library(ggalt)  #geom_encircle
#library(ggpubr)  #link the points
library(ggtext)
##Weighted----------------------------
wunifrac1 <- read.csv('b-diversity/wunifrac_info_phy_3.csv',row.names = 1,header = T)
wunifrac2 <- read.csv('b-diversity/wunifrac_info_phy_5.csv',row.names = 1,header = T)
wunifrac3 <- read.csv('b-diversity/wunifrac_info_phy_7.csv',row.names = 1,header = T)
wunifrac4 <- read.csv('b-diversity/wunifrac_info_phy_14.csv',row.names = 1,header = T)
wunifrac5 <- read.csv('b-diversity/wunifrac_info_phy_21.csv',row.names = 1,header = T)
wunifrac <- rbind(wunifrac1,wunifrac2,wunifrac3,wunifrac4,wunifrac5)

wunifrac$treatment <- factor(wunifrac$treatment,levels = c("Mock","EM"))
wunifrac$day <- factor(wunifrac$day,levels = c("3 dpi","5 dpi","7 dpi","14 dpi","21 dpi"))
my_color_color <- c(Mock="#377483", EM="#B83945")
my_color_fill <- c(Mock="#C7DFF0", EM="#FBDFE2")


# annotation information 
annotation_data <- data.frame(
  day = c("3 dpi", "5 dpi", "7 dpi", "14 dpi", "21 dpi"),
  x = c(-0.05, -0.05, -0.05, -0.05, -0.05),
  y = c(0.5, 0.5, -0.55, -0.55, 0.5),
  label = c(
    "<i>P</i> = 0.421<br>R<sup>2</sup> = 0.046",
    "<i>P</i> = 0.024<br>R<sup>2</sup> = 0.124", 
    "<i>P</i> < 0.001<br>R<sup>2</sup> = 0.348", 
    "<i>P</i> = 0.100<br>R<sup>2</sup> = 0.090", 
    "<i>P</i> = 0.440<br>R<sup>2</sup> = 0.053"
  )
)

annotation_data$day <- factor(annotation_data$day,levels = c("3 dpi","5 dpi","7 dpi","14 dpi","21 dpi"))


p1 <- ggplot(wunifrac,aes(x=axis1,axis2,color=treatment)) +
  geom_point(size=2) + 
  facet_grid(.~day)+
  scale_color_manual(values = my_color_color)+
  labs(x='PC1 (49.29%)',y='PC2 (27.79%)',color="",title = "Weighted Unifrac")+
  scale_x_continuous(limits = c(-0.6,0.6),breaks = c(-0.3,0.3))+
  scale_y_continuous(limits = c(-0.7,0.7),breaks = c(-0.3,0.3))+
  stat_ellipse(level=0.8, show.legend = FALSE)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = 8,vjust = 0.9,color = "black",margin = margin(t=0)),
        axis.text.y = element_text(size = 8,color = "black"),
        panel.border = element_rect(linewidth = 0.5,color = 'black'),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        #axis.ticks.x = element_line(linewidth = 0.3), 
        #axis.ticks.y = element_line(linewidth  = 0.3),
        legend.text = element_text(size = 9,colour = "black",margin = margin(l=-1)),
        legend.spacing.x = unit(0.05, "cm"),
        legend.margin = margin(l=-10),
        #legend.position = c(1.1, 0.6),         
        #legend.justification = c("right"), 
        #plot.margin = unit(c(0, 20, 0, 0), "mm"), # top right bottom left
        plot.title = element_text(hjust = 0.5,size = 9,margin = margin(b=3)),
        axis.title = element_text(size=8),
        strip.text = element_text(size = 9, color = "black"), #adjust facet_grid text
        panel.spacing.x= unit(0.1, "lines") #adjust facet_grid horizontal distance
        )+
  geom_richtext(data = annotation_data,
            aes(x = x, y = y, label = label),
            size = 3, inherit.aes = FALSE,fill = NA,label.colour = NA)
p1

library(ggh4x)
p1 <- p1 + force_panelsizes(total_width=unit(3.3, "in"),total_height=unit(1.0, "in"))
p1 
ggsave("b-diversity/wunifrac.pdf",p1,width = 3.3*1.3,height = 1*2)
ggsave("b-diversity/wunifrac.png",p1,dpi=600,width = 3.3*1.3,height = 1*2)


##Unweighted--------------------------
unifrac1 <- read.csv('b-diversity/unifrac_info_phy_3.csv',row.names = 1,header = T)
unifrac2 <- read.csv('b-diversity/unifrac_info_phy_5.csv',row.names = 1,header = T)
unifrac3 <- read.csv('b-diversity/unifrac_info_phy_7.csv',row.names = 1,header = T)
unifrac4 <- read.csv('b-diversity/unifrac_info_phy_14.csv',row.names = 1,header = T)
unifrac5 <- read.csv('b-diversity/unifrac_info_phy_21.csv',row.names = 1,header = T)
unifrac <- rbind(unifrac1,unifrac2,unifrac3,unifrac4,unifrac5)

unifrac$treatment <- factor(unifrac$treatment,levels = c("Mock","EM"))
unifrac$day <- factor(unifrac$day,levels = c("3 dpi","5 dpi","7 dpi","14 dpi","21 dpi"))

# annotation information
annotation_data <- data.frame(
  day = c("3 dpi","5 dpi","7 dpi","14 dpi","21 dpi"), 
  x = c(-0.05,-0.05,-0.05,-0.05,-0.05),      
  y = c(0.5,0.5,0.5,0.5,0.5),       
  label = c(
    "<i>P</i> = 0.169<br>R<sup>2</sup> = 0.064",
    "<i>P</i> = 0.003<br>R<sup>2</sup> = 0.155", 
    "<i>P</i> < 0.001<br>R<sup>2</sup> = 0.194", 
    "<i>P</i> = 0.013<br>R<sup>2</sup> = 0.094", 
    "<i>P</i> < 0.001<br>R<sup>2</sup> = 0.138"
  )
)
  
annotation_data$day <- factor(annotation_data$day,levels = c("3 dpi","5 dpi","7 dpi","14 dpi","21 dpi"))
  

p1 <- ggplot(unifrac,aes(x=axis1,axis2,color=treatment)) +
  geom_point(size=2) + 
  facet_grid(.~day)+
  scale_color_manual(values = my_color_color)+
  labs(x='PC1 (21.25%)',y='PC2 (16.87%)',color="",title = "Unweighted Unifrac")+
  scale_x_continuous(limits = c(-0.6,0.6),breaks = c(-0.3,0.3))+
  scale_y_continuous(limits = c(-0.7,0.7),breaks = c(-0.3,0.3))+
  stat_ellipse(level=0.8, show.legend = FALSE)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = 8,vjust = 0.9,color = "black",margin = margin(t=0)),
        axis.text.y = element_text(size = 8,color = "black"),
        panel.border = element_rect(linewidth = 0.5,color = 'black'),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        #axis.ticks.x = element_line(linewidth = 0.3), 
        #axis.ticks.y = element_line(linewidth  = 0.3),
        legend.text = element_text(size = 9,colour = "black",margin = margin(l=-1)),
        legend.spacing.x = unit(0.05, "cm"),
        legend.margin = margin(l=-10),
        #legend.position = c(1.1, 0.6),         
        #legend.justification = c("right"), 
        #plot.margin = unit(c(0, 20, 0, 0), "mm"), # top right bottom left
        plot.title = element_text(hjust = 0.5,size = 9,margin = margin(b=3)),
        axis.title = element_text(size=8),
        strip.text = element_text(size = 9, color = "black"), #adjust facet_grid text
        panel.spacing.x= unit(0.1, "lines") #adjust facet_grid horizontal distance
  )+
  geom_richtext(data = annotation_data,
                aes(x = x, y = y, label = label),
                size = 3, inherit.aes = FALSE,fill = NA,label.colour = NA)
p1

library(ggh4x)
p1 <- p1 + force_panelsizes(total_width=unit(3.3, "in"),total_height=unit(1.0, "in"))
p1 
ggsave("b-diversity/unifrac.pdf",p1,width = 3.3*1.3,height = 1*2)
ggsave("b-diversity/unifrac.png",p1,dpi=600,width = 3.3*1.3,height = 1*2)


# Microbiota composition-----------------
library(tidyverse)
library(paletteer)
library(ggtext)
library(glue)
###  Pick the colors  ###
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


col1 <- paletteer_d("ggsci::category20_d3")
col2 <- paletteer_d("Polychrome::kelly")[7:15]
my_color <- c(col1[1:7],col2,col1[8:20])


#Import data
sample_data <- read.csv('metadata.csv',header = T,row.names = 1)
sample_data$rownames <- rownames(sample_data)

###  Phylum  ------------------
phylum <- read.csv('composition/phylum.csv',header = T,row.names = 1)
phy_name <- rownames(phylum)
phylum <- as.data.frame(t(phylum))
phylum$rownames <- rownames(phylum)

phy_pic <- merge(phylum, sample_data, by='rownames')
phy_pic <- pivot_longer(phy_pic,cols = colnames(phylum)[-ncol(phylum)],
                        names_to = 'name',values_to = 'abundance')

phy_pic <- phy_pic %>%
  group_by(treatment,name,day) %>%
  summarise(mean=mean(abundance))

phy_pic$treatment <- factor(phy_pic$treatment,levels = c("Mock","EM"))
phy_pic$day <- factor(phy_pic$day,levels = c("3 dpi","5 dpi","7 dpi","14 dpi","21 dpi"))
phy_pic$name <- factor(phy_pic$name,levels = rev(phy_name))
phy_col<- c(my_color[1:length(unique(phy_pic$name))])
phy_pic$mean <- phy_pic$mean*100

# x_axis info
x_axis_Mock <- scale_x_discrete(
  labels = c("3", "5", "7", "14", "21")
)
x_axis_EM <- scale_x_discrete(
  labels = c("3", "5", "7", "14", "21dpi")
)

# plot
p1 <- ggplot(phy_pic,aes(x=day,y=mean,fill=name))+
  facet_grid(.~treatment,scales = "free_x")+
  geom_col()+
  facetted_pos_scales(
    x = list(
      treatment == "Mock" ~ x_axis_Mock,
      treatment == "EM" ~ x_axis_EM
    )) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = 'Relative Abundance (%)',fill='Phylum') +
  scale_fill_manual(values = rev(phy_col)) +  
  scale_color_manual(values = rev(phy_col)) +
  theme_bw()+
  guides(color = guide_legend(ncol = 1), fill = guide_legend(ncol = 1))+
  theme(
        panel.background = element_rect(fill = "gray90"),
        panel.grid = element_blank(),
        panel.border = element_rect(linewidth = 0.5,color = 'white'),
        axis.text.x = element_text(size = 8,angle = 0,vjust = 1,
                                   hjust = 0.5,color = "black",
                                   margin = margin(t=1)),
        axis.text.y = element_text(size = 8,color = "black"),
        axis.title.y = element_text(margin = margin(r = 4), size = 10, color = "black"), 
        axis.ticks.x = element_blank(),  
        axis.ticks.y = element_blank(),  
        legend.key.size = unit(0.3, "cm"),  
        legend.text = element_markdown(size =7, color = "black",margin = margin(l=1)), 
        legend.margin = margin(l=-8),
        legend.spacing.x = unit(0.02, "cm"), 
        legend.spacing.y = unit(0.02, "cm"), 
        legend.title = element_text(size = 8, margin = margin(b=1)),
        strip.background = element_rect(linetype = "solid",
                                        color = "white",fill = "gray90", 
                                        linewidth = 0.5),
        strip.text = element_text(size = 9),
        panel.spacing.x= unit(0.15, "lines")
  )
library(ggh4x)
p1 <- p1 + force_panelsizes(total_width=unit(2.5, "in"),total_height=unit(3, "in"))
p1 
ggsave("composition/phylum.pdf",p1,width = 2.5*1.9,height = 3*1.2)
ggsave("composition/phylum.png",p1,width = 2.5*1.9,height = 3*1.2,dpi = 600)

###  Family  ------------------- 
family <- read.csv('composition/family_top10.csv',header = T,row.names = 1)
fam_name <- rownames(family)
family <- as.data.frame(t(family))
family$rownames <- rownames(family)
fam_pic <- merge(family, sample_data, by='rownames')

fam_pic <- pivot_longer(fam_pic,cols = colnames(family)[-ncol(family)],
                        names_to = 'name',values_to = 'abundance')


fam_pic <- fam_pic %>%
  group_by(treatment,name,day) %>%
  summarise(mean=mean(abundance))


### Make legend name ###
# all family name
# fam_pic <- mutate(fam_pic, font_name=glue("<i>{name}"))
# fam_pic$font_name <- gsub("<i>Others","Others",fam_pic$font_name)

# have unclassified family name
name <- as.data.frame(t(as.data.frame(strsplit(fam_pic$name,"_"))))
for (i in 1:nrow(name)){
  if (name$V1[i] == name$V2[i]) {
    name$V2[i] = ""
  }
}
name$V2 <- gsub("Unclassified","_Unclassified",name$V2)

fam_pic$pre <- name$V1
fam_pic$post <- name$V2
fam_pic <- mutate(fam_pic, font_name=glue("<i>{pre}</i>{post}"))

fam_pic$font_name <- gsub("<i>Bacillota</i>_Unclassified","Bacillota_Unclassified",fam_pic$font_name)
fam_pic$font_name <- gsub("<i>Eubacteriales</i>_Unclassified","Eubacteriales_Unclassified",fam_pic$font_name)

fam_pic$font_name <- gsub("<i>Others</i>","Others",fam_pic$font_name)
### Over ###

fam_pic$treatment <- factor(fam_pic$treatment,levels = c("Mock","EM"))
fam_pic$day <- factor(fam_pic$day,levels = c("3 dpi","5 dpi","7 dpi","14 dpi","21 dpi"))
fam_pic$name <- factor(fam_pic$name,levels = rev(fam_name))
fam_pic <- fam_pic[order(fam_pic$name),]
fam_col <- c(my_color[1:(length(unique(fam_pic$name))-1)],'gray')
fam_pic$mean <- fam_pic$mean*100

# x_axis info
x_axis_Mock <- scale_x_discrete(
  labels = c("3", "5", "7", "14", "21")
)
x_axis_EM <- scale_x_discrete(
  labels = c("3", "5", "7", "14", "21dpi")
)

# plot
p2 <- ggplot(fam_pic,aes(x=day,y=mean,fill=name))+
  facet_grid(.~treatment,scales = "free_x")+
  geom_col()+
  facetted_pos_scales(
    x = list(
      treatment == "Mock" ~ x_axis_Mock,
      treatment == "EM" ~ x_axis_EM
    )) +
  scale_x_discrete(labels = c(3,5,7,14,21),expand = c(0, 0))+
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = 'Relative Abundance (%)',fill='Family') +
  scale_fill_manual(values = rev(fam_col),label=unique(fam_pic$font_name)) +  
  scale_color_manual(values = rev(fam_col),label=unique(fam_pic$font_name)) +
  theme_bw()+
  guides(color = guide_legend(ncol = 1), fill = guide_legend(ncol = 1))+
  theme(
        panel.background = element_rect(fill = "gray90"),
        panel.grid = element_blank(),
        panel.border = element_rect(linewidth = 0.5,color = 'white'),
        axis.text.x = element_text(size = 8,angle = 0,vjust = 1,
                                   hjust = 0.5,color = "black",
                                   margin = margin(t=1)),
        axis.text.y = element_text(size = 8,color = "black"),
        axis.title.y = element_text(margin = margin(r = 4), size = 10, color = "black"), 
        axis.ticks.x = element_blank(),  
        axis.ticks.y = element_blank(),  
        legend.key.size = unit(0.3, "cm"),  
        legend.text = element_markdown(size =7, color = "black",margin = margin(l=1)), 
        legend.margin = margin(l=-8),
        legend.spacing.x = unit(0.02, "cm"), 
        legend.spacing.y = unit(0.02, "cm"), 
        legend.title = element_text(size = 8, margin = margin(b=1)),
        strip.background = element_rect(linetype = "solid",
                                        color = "white", fill = "gray90",
                                        linewidth = 0.5),
        strip.text = element_text(size = 9),
        panel.spacing.x= unit(0.15, "lines")
  )
p2

library(ggh4x)
p2 <- p2 + force_panelsizes(total_width=unit(2.5, "in"),total_height=unit(3, "in"))
p2 
ggsave("composition/family.pdf",p2,width = 2.5*1.9,height = 3*1.2)
ggsave("composition/family.png",p2,width = 2.5*1.9,height = 3*1.2,dpi = 600)

###  Genus  ------------------
genus <- read.csv('composition/genus_top20.csv',header = T,row.names = 1)
gen_name <- rownames(genus)
genus <- as.data.frame(t(genus))
genus$rownames <- rownames(genus)
gen_pic <- merge(genus, sample_data, by='rownames')
gen_pic <- pivot_longer(gen_pic,cols = colnames(genus)[-ncol(genus)],
                        names_to = 'name',values_to = 'abundance')

gen_pic <- gen_pic %>%
  group_by(treatment,name,day) %>%
  summarise(mean=mean(abundance))

### Make legend name ###
name <- as.data.frame(t(as.data.frame(strsplit(gen_pic$name,"_"))))
for (i in 1:nrow(name)){
  if (name$V1[i] == name$V2[i]) {
    name$V2[i] = ""
  }
}
name$V2 <- gsub("Unclassified","_Unclassified",name$V2)

gen_pic$pre <- name$V1
gen_pic$post <- name$V2
gen_pic <- mutate(gen_pic, font_name=glue("<i>{pre}</i>{post}"))
gen_pic$font_name <- gsub("<i>Others</i>","Others",gen_pic$font_name)

gen_pic$font_name <- gsub("<i>Bacillota</i>_Unclassified","Bacillota_Unclassified",gen_pic$font_name)
gen_pic$font_name <- gsub("<i>Eubacteriales</i>_Unclassified","Eubacteriales_Unclassified",gen_pic$font_name)
### Over ###

gen_pic$treatment <- factor(gen_pic$treatment,levels = c("Mock","EM"))
gen_pic$day <- factor(gen_pic$day,levels = c("3 dpi","5 dpi","7 dpi","14 dpi","21 dpi"))
gen_pic$name <- factor(gen_pic$name,levels = rev(gen_name))
gen_pic <- gen_pic[order(gen_pic$name),]
gen_col <- c(my_color[1:(length(unique(gen_pic$name))-1)],'gray')
gen_pic$mean <- gen_pic$mean*100

# x_axis info
x_axis_Mock <- scale_x_discrete(
  labels = c("3", "5", "7", "14", "21")
)
x_axis_EM <- scale_x_discrete(
  labels = c("3", "5", "7", "14", "21dpi")
)

# plot
p3 <- ggplot(gen_pic,aes(x=day,y=mean,fill=name))+
  facet_grid(.~treatment,scales = "free_x")+
  geom_col()+
  facetted_pos_scales(
    x = list(
      treatment == "Mock" ~ x_axis_Mock,
      treatment == "EM" ~ x_axis_EM
    )) +
  scale_x_discrete(labels = c(3,5,7,14,21),expand = c(0, 0))+
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = 'Relative Abundance (%)',fill='Genus') +
  scale_fill_manual(values = rev(gen_col),label=unique(gen_pic$font_name)) +  
  scale_color_manual(values = rev(gen_col),label=unique(gen_pic$font_name)) +
  theme_bw()+
  guides(color = guide_legend(ncol = 1), fill = guide_legend(ncol = 1))+
  theme(
    panel.background = element_rect(fill = "gray90"),
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth = 0.5,color = 'white'),
    axis.text.x = element_text(size = 8,angle = 0,vjust = 1,
                               hjust = 0.5,color = "black",
                               margin = margin(t=1)),
    axis.text.y = element_text(size = 8,color = "black"),
    axis.title.y = element_text(margin = margin(r = 4), size = 10, color = "black"), 
    axis.ticks.x = element_blank(),  
    axis.ticks.y = element_blank(),  
    legend.key.size = unit(0.3, "cm"),  
    legend.text = element_markdown(size =7, color = "black",margin = margin(l=1)), 
    legend.margin = margin(l=-8),
    legend.spacing.x = unit(0.02, "cm"), 
    legend.spacing.y = unit(0.02, "cm"), 
    legend.title = element_text(size = 8, margin = margin(b=1)),
    strip.background = element_rect(linetype = "solid",
                                    color = "white", fill = "gray90",
                                    linewidth = 0.5),
    strip.text = element_text(size = 9),
    panel.spacing.x= unit(0.15, "lines")
  )
p3

library(ggh4x)
p3 <- p3 + force_panelsizes(total_width=unit(2.5, "in"),total_height=unit(3, "in"))
p3 
ggsave("composition/genus.pdf",p3,width = 2.5*1.9,height = 3*1.2)
ggsave("composition/genus.png",p3,width = 2.5*1.9,height = 3*1.2,dpi = 600)

###  Feature  -----------------------------
feature <- read.csv('composition/feature_top25.csv',header = T,row.names = 1)
#re-order the feature name
feature$name <- rownames(feature)
feature <- separate(feature,name,into = c("species","order_num"),sep = "_F")
feature <- feature[order(as.numeric(feature$order_num)),]
feature <- feature[,-c(ncol(feature)-1,ncol(feature))]

#modify data format
fea_name <- rownames(feature)  #using for setting order
feature <- as.data.frame(t(feature))
feature$rownames <- rownames(feature)
fea_pic <- merge(feature, sample_data, by='rownames')
fea_pic <- pivot_longer(fea_pic,cols = colnames(feature)[-ncol(feature)],
                        names_to = 'name',values_to = 'abundance')

fea_pic <- fea_pic %>%
  group_by(treatment,day,name) %>%
  summarise(mean=mean(abundance))


### Make legend name ###
# get the name
name <- as.data.frame(t(as.data.frame(strsplit(fea_pic$name,"_F"))))
for (i in 1:nrow(name)){
  if (name$V1[i] == name$V2[i]) {
    name$V2[i] = ""
  }
}
# split strings according different situations
name$V3 <- NA
name$V4 <- NA
for (i in 1:nrow(name)){
  if (grepl("Unclassified", name$V1[i])){
    name$V3[i] <- str_split_1(name$V1[i],"_U")[1]
    name$V4[i] <- str_split_1(name$V1[i],"_U")[2]
  }
}
# make font name
name <- name %>%
  mutate(
    font_name = case_when(
      V4 == "nclassified" ~ glue("<i>{V3}</i>_U{V4}_F{V2}"),
      TRUE ~ glue("<i>{V1}</i>_F{V2}")
    )
  )

fea_pic$font_name <- name$font_name
fea_pic$font_name <- gsub("<i>Others</i>_F","Others",fea_pic$font_name)

### Over ###

fea_pic$treatment <- factor(fea_pic$treatment,levels = c("Mock","EM"))
fea_pic$day <- factor(fea_pic$day,levels = c("3 dpi","5 dpi","7 dpi","14 dpi","21 dpi"))
fea_pic$name <- factor(fea_pic$name,levels = rev(fea_name))
fea_pic <- fea_pic[order(fea_pic$name),]
fea_col <- c(my_color[1:(length(unique(fea_pic$name))-1)],'gray')
fea_pic$mean <- fea_pic$mean*100


# x_axis info
x_axis_Mock <- scale_x_discrete(
  labels = c("3", "5", "7", "14", "21")
)
x_axis_EM <- scale_x_discrete(
  labels = c("3", "5", "7", "14", "21dpi")
)

# plot
p4 <- ggplot(fea_pic,aes(x=day,y=mean,fill=name))+
  facet_grid(.~treatment,scales = "free_x")+
  geom_col()+
  facetted_pos_scales(
    x = list(
      treatment == "Mock" ~ x_axis_Mock,
      treatment == "EM" ~ x_axis_EM
    )) +
  scale_x_discrete(labels = c(3,5,7,14,21),expand = c(0, 0))+
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = 'Relative Abundance (%)',fill='ASVs') +
  scale_fill_manual(values = rev(fea_col),label=unique(fea_pic$font_name)) +  
  scale_color_manual(values = rev(fea_col),label=unique(fea_pic$font_name)) +
  theme_bw()+
  guides(color = guide_legend(ncol = 1), fill = guide_legend(ncol = 1))+
  theme(
    panel.background = element_rect(fill = "gray90"),
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth = 0.5,color = 'white'),
    axis.text.x = element_text(size = 8,angle = 0,vjust = 1,
                               hjust = 0.5,color = "black",
                               margin = margin(t=1)),
    axis.text.y = element_text(size = 8,color = "black"),
    axis.title.y = element_text(margin = margin(r = 4), size = 10, color = "black"), 
    axis.ticks.x = element_blank(),  
    axis.ticks.y = element_blank(),  
    legend.key.size = unit(0.3, "cm"),  
    legend.text = element_markdown(size =7, color = "black",margin = margin(l=1)), 
    legend.margin = margin(l=-8),
    legend.spacing.x = unit(0.02, "cm"), 
    legend.spacing.y = unit(0.02, "cm"), 
    legend.title = element_text(size = 8, margin = margin(b=1)),
    strip.background = element_rect(linetype = "solid",
                                    color = "white", fill = "gray90",
                                    linewidth = 0.5),
    strip.text = element_text(size = 9),
    panel.spacing.x= unit(0.15, "lines")
  )
p4

library(ggh4x)
p4 <- p4 + force_panelsizes(total_width=unit(2.5, "in"),total_height=unit(3, "in"))
p4 
ggsave("composition/feature.pdf",p4,width = 2.5*2.0,height = 3*1.2)
ggsave("composition/feature.png",p4,width = 2.5*2.0,height = 3*1.2,dpi = 600)


# Differential enrichment analysis-----------
library(phyloseq)
library(ANCOMBC)
library(tidyverse)
dir.create("deg")
f_table <- "feature_table_filter.csv"
s_table <- "metadata.csv"
t_table <- "taxonomy.csv"
feature_table <- read.csv(f_table,header = T,row.names = 1)
taxonomy <- read.csv(t_table,header = T,row.names = 1)
sample_data <- read.csv(s_table,header = T,row.names = 1)
#Transform feature table and taxonomy table to matrix
sample_data$day <- factor(sample_data$day,levels = c("3 dpi","5 dpi","7 dpi","14 dpi","21 dpi"))
sample_data$treatment <- factor(sample_data$treatment,levels = c("Mock","EM"))
OTU<- otu_table(as.matrix(feature_table),taxa_are_rows = T)
TAX<- tax_table(as.matrix(taxonomy))
SAM<-sample_data(sample_data)
physeq<- phyloseq(OTU,TAX,SAM)

# subset data 
phy_3 <- subset_samples(physeq,day=="3 dpi")
phy_5 <- subset_samples(physeq,day=="5 dpi")
phy_7 <- subset_samples(physeq,day=="7 dpi")
phy_14 <- subset_samples(physeq,day=="14 dpi")
phy_21 <- subset_samples(physeq,day=="21 dpi")
phy_mock <- subset_samples(physeq, treatment=="Mock")
phy_em <- subset_samples(physeq, treatment=="EM"|(treatment=="Mock"&day=="3 dpi"))
x1 <- as.data.frame(phy_em@otu_table)
x2 <- as.data.frame(phy_em@tax_table)
x3 <- as.data.frame(phy_em@sam_data)
x3$group <- paste(x3$day,x3$treatment,sep = "_")
x3$group <- factor(x3$group,levels=c("3 dpi_Mock","3 dpi_EM","5 dpi_EM", "7 dpi_EM",   "14 dpi_EM","21 dpi_EM"))
OTU<- otu_table(as.matrix(x1),taxa_are_rows = T)
TAX<- tax_table(as.matrix(x2))
SAM<-sample_data(x3)
phy_em<- phyloseq(OTU,TAX,SAM)
# modify tax table
taxonomy <- taxonomy %>% 
  rownames_to_column(var = "taxon")
## For ASVs--------------------
library(crayon)
# for different dpi datasets
for (i in c(3,5,7,14,21)){
  if (i == 3){
    name = phy_3
  } else if (i == 5){
    name = phy_5
  } else if (i == 7){
    name = phy_7
  } else if (i == 14){
    name = phy_14
  } else if (i == 21){
    name = phy_21
  }
  file_name=paste("deg/res_prim_",i,".csv",sep = "")
  message(bold(blue(paste0("--- Run analysis for ", i, " dpi dataset ---"))))
  set.seed(123)
  output <- ancombc2(data=name,
                     fix_formula = "treatment",group = "treatment",tax_level="Species",
                     prv_cut = 0,lib_cut = 0, alpha=0.05,p_adj_method="fdr",
                     struc_zero = F,neg_lb = F,pseudo = 1,pseudo_sens = FALSE,
                     global = FALSE,pairwise = TRUE,dunnet = FALSE,trend = FALSE,
                     n_cl = 3
  )
  # export results
  # select top ASVs
  feature_table2 <- feature_table %>% 
    mutate(across(everything(), ~./sum(.))) %>% 
    arrange(desc(rowSums(.))) 
  feature_table2 <- feature_table2[1:50,]
  # get results
  res_prim = output$res %>% 
    select(c(1,3,q_treatmentEM)) %>% 
    filter(taxon %in% rownames(feature_table2)) 
  res_prim <- left_join(res_prim,taxonomy[,c(1,6,8)],by="taxon")
  write.csv(res_prim,file = file_name)
}

# for mock dataset
set.seed(123)
output <- ancombc2(data=phy_mock,
                   fix_formula = "day",group = "day",tax_level="Species",
                   prv_cut = 0,lib_cut = 0, alpha=0.05,p_adj_method="fdr",
                   struc_zero = F,neg_lb = F,pseudo = 1,pseudo_sens = FALSE,
                   global = FALSE,pairwise = FALSE,dunnet = FALSE,trend = FALSE,
                   n_cl = 3
)
res_prim = output$res %>% 
  select(1, starts_with("lfc_d"), starts_with("q_d")) %>% 
  filter(taxon %in% rownames(feature_table2)) 
res_prim <- left_join(res_prim,taxonomy[,c(1,6,8)],by="taxon")
write.csv(res_prim,file = "deg/res_prim_mock.csv")

# for em dataset
set.seed(123)
output <- ancombc2(data=phy_em,
                   fix_formula = "group",group = "group",tax_level="Species",
                   prv_cut = 0,lib_cut = 0, alpha=0.05,p_adj_method="fdr",
                   struc_zero = F,neg_lb = F,pseudo = 1,pseudo_sens = FALSE,
                   global = FALSE,pairwise = FALSE,dunnet = FALSE,trend = FALSE,
                   n_cl = 3
)
res_prim = output$res %>% 
  select(1, starts_with("lfc_g"), starts_with("q_g")) %>% 
  filter(taxon %in% rownames(feature_table2)) 
res_prim <- left_join(res_prim,taxonomy[,c(1,6,8)],by="taxon")
write.csv(res_prim,file = "deg/res_prim_em.csv")

#combine results
# import
r3 <- read.csv("deg/res_prim_3.csv",header = T,row.names = 1) %>% 
  rename("lfc_3_EM/Mock" = lfc_treatmentEM, "q_3_EM/Mock" = q_treatmentEM)
r5 <- read.csv("deg/res_prim_5.csv",header = T,row.names = 1) %>% 
  rename("lfc_5_EM/Mock" = lfc_treatmentEM, "q_5_EM/Mock" = q_treatmentEM) %>% 
  select(-4,-5)
r7 <- read.csv("deg/res_prim_7.csv",header = T,row.names = 1) %>% 
  rename("lfc_7_EM/Mock" = lfc_treatmentEM, "q_7_EM/Mock" = q_treatmentEM) %>% 
  select(-4,-5)
r14 <- read.csv("deg/res_prim_14.csv",header = T,row.names = 1) %>% 
  rename("lfc_14_EM/Mock" = lfc_treatmentEM, "q_14_EM/Mock" = q_treatmentEM) %>% 
  select(-4,-5)
r21 <- read.csv("deg/res_prim_21.csv",header = T,row.names = 1) %>% 
  rename("lfc_21_EM/Mock" = lfc_treatmentEM, "q_21_EM/Mock" = q_treatmentEM) %>% 
  select(-4,-5)
mock <- read.csv("deg/res_prim_mock.csv",header = T,row.names = 1) %>% 
  rename("lfc_5_Mock" = lfc_day5.dpi, "lfc_7_Mock" = lfc_day7.dpi,
         "lfc_14_Mock" = lfc_day14.dpi, "lfc_21_Mock" = lfc_day21.dpi,
         "q_5_Mock" = q_day5.dpi, "q_7_Mock" = q_day7.dpi,
         "q_14_Mock" = q_day14.dpi, "q_21_Mock" = q_day21.dpi) %>% 
  select(-10,-11) %>% 
  mutate("lfc_3_Mock"=0,"q_3_Mock"=1)
em <- read.csv("deg/res_prim_em.csv",header = T,row.names = 1) %>% 
  rename("lfc_3_EM" = lfc_group3.dpi_EM, "lfc_5_EM" = lfc_group5.dpi_EM,
         "lfc_7_EM" = lfc_group7.dpi_EM, "lfc_14_EM" = lfc_group14.dpi_EM,
         "lfc_21_EM" = lfc_group21.dpi_EM,
         "q_3_EM" = q_group3.dpi_EM, "q_5_EM" = q_group5.dpi_EM,
         "q_7_EM" = q_group7.dpi_EM, "q_14_EM" = q_group14.dpi_EM,
         "q_21_EM" = q_group21.dpi_EM,) %>% 
  select(-12,-13)

#combine
results <- left_join(r3,r5,by="taxon") %>% 
  left_join(.,r7,by="taxon") %>% 
  left_join(.,r14,by="taxon") %>% 
  left_join(.,r21,by="taxon") %>% 
  left_join(.,mock,by="taxon") %>% 
  left_join(.,em,by="taxon") %>% 
  select(1,4,5,everything())

write.csv(results,file = "deg/ancombc_species.csv")


## fig ---------------------
library(tidyverse)
library(glue)
library(ggtext)
library(ggh4x)

# define Family color
my_color <- c("Oscillospiraceae"="#1F4880","Lactobacillaceae"="#004300", 
              "Lachnospiraceae"="#009100","Peptostreptococcaceae"="#91B585",
              "Erysipelotrichaceae"="#7B24B5","Enterobacteriaceae"="#D13D20",
              "Coprobacillaceae"="#D2A517","Clostridiaceae"="#DE7C68",
              "Bacillota_Unclassified"="#7E8F8F"
              )

# import data
df_wide <- read.csv("deg/ancombc_species.csv",header = T,row.names = 1)
df <- pivot_longer(df_wide,
                   cols = starts_with(c("lfc_", "q_")),
                   names_sep = "_",
                   names_to = c(".value", "day", "cate"))

### Make legend name for Species -----------------------------
# get the name
name <- df %>% 
  select(c(Family,Species)) %>% 
  separate(Species,c("V1","V2"),sep = "_F")
name$color <- NA
# define color for label
for (i in 1:nrow(name)) {
  if (name$Family[i] == "Lactobacillaceae"){
    name$color[i] = "#004300"
  } else if (name$Family[i] == "Lachnospiraceae"){
    name$color[i] = "#009100"
  } else if (name$Family[i] == "Oscillospiraceae"){
    name$color[i] = "#1F4880"
  } else if (name$Family[i] == "Coprobacillaceae"){
    name$color[i] = "#D2A517"
  } else if (name$Family[i] == "Enterobacteriaceae"){
    name$color[i] = "#D13D20"
  } else if (name$Family[i] == "Peptostreptococcaceae"){
    name$color[i] = "#91B585"
  } else if (name$Family[i] == "Erysipelotrichaceae"){
    name$color[i] = "#7B24B5"
  } else if (name$Family[i] == "Clostridiaceae"){
    name$color[i] = "#DE7C68"
  } else if (name$Family[i] == "Bacillota_Unclassified"){
    name$color[i] = "#7E8F8F"
  }
}
# split strings according different situations
name$V3 <- NA
name$V4 <- NA
for (i in 1:nrow(name)){
  if (grepl("Unclassified", name$V1[i])){
    name$V3[i] <- str_split_1(name$V1[i],"_U")[1]
    name$V4[i] <- str_split_1(name$V1[i],"_U")[2]
  }
}
# make font name
name <- name %>%
  mutate(
    font_name = case_when(
      V4 == "nclassified" ~ glue("<span style='color:{color}'><i>{V3}</i>_U{V4}_F{V2}"),
      TRUE ~ glue("<span style='color:{color}'><i>{V1}</i>_F{V2}")
    )
  )
# modify specific name
name$font_name <-  gsub(
  pattern ="<span style='color:{color}'><i>Bacillota</i>_Unclassified_F39",
  replacement = "<span style='color:{color}'>Bacillota_Unclassified_F39",
  x = name$font_name,
  fixed = TRUE)

# add to fig data and check name
df$font_name <- name$font_name
name <- name %>% 
  unite("Species",c(V1,V2),sep = "_F")
identical(name$Species, df$Species)   # MUST output TRUE

### make legend name for family -----------------
name <- df %>% 
  select(Family) %>%
  mutate(
    font_name = case_when(
      grepl("Unclassified", Family) ~ glue("{Family}"),
      TRUE ~ glue("<i>{Family}</i>")
    )
  )
# add to fig data and check name
df$font_Family <- name$font_name
identical(name$Family, df$Family)   # MUST output TRUE

# re-define Family color
my_color2 <- c("<i>Oscillospiraceae</i>"="#1F4880","<i>Lactobacillaceae</i>"="#004300", 
               "<i>Lachnospiraceae</i>"="#009100","<i>Peptostreptococcaceae</i>"="#91B585",
               "<i>Erysipelotrichaceae</i>"="#7B24B5","<i>Enterobacteriaceae</i>"="#D13D20",
               "<i>Coprobacillaceae</i>"="#D2A517","<i>Clostridiaceae</i>"="#DE7C68",
               "Bacillota_Unclassified"="#7E8F8F"
)

# define Family legend order
df$font_Family <- factor(df$font_Family,levels = c(
  "<i>Enterobacteriaceae</i>","<i>Clostridiaceae</i>",
  "<i>Coprobacillaceae</i>","<i>Erysipelotrichaceae</i>",
  "Bacillota_Unclassified","<i>Oscillospiraceae</i>",
  "<i>Peptostreptococcaceae</i>","<i>Lachnospiraceae</i>",
  "<i>Lactobacillaceae</i>"
))

### translate q value to asterisk -----------------
df <- df %>%
  mutate(
    sig = case_when(
      q < 0.001 ~ "***",
      q < 0.01 ~ "**",
      q < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

# calculate log2FC
df$lfc <- df$lfc / log(2)

#adjust order in the fig
df$font_name <- factor(df$font_name,levels = rev(unique(df$font_name)))
df$day <- factor(df$day ,levels = c("3",  "5",  "7",  "14", "21"))
df$cate <- factor(df$cate,levels = c("Mock","EM","EM.Mock"),
                  labels = c("Mock","EM","EM/Mock"))

# x_axis info
x_axis_Mock <- scale_x_discrete(
  labels = c("3", "5", "7", "14", "21")
)
x_axis_EM <- scale_x_discrete(
  labels = c("3", "5", "7", "14", "21")
)

x_axis_Mock_EM <- scale_x_discrete(
  labels = c("3", "5", "7", "14", "21dpi")
)


### draw the fig -------------
range(df$lfc)
p1 <- ggplot(df,aes(y=font_name,x=day))+
  geom_tile(aes(fill=lfc),color="gray96")+
  facet_grid(.~cate,scales = "free_x")+
  facetted_pos_scales(
    x = list(
      cate == "Mock" ~ x_axis_Mock,
      cate == "EM" ~ x_axis_EM,
      cate == "EM/Mock" ~ x_axis_Mock_EM
    )) +
  scale_fill_gradient2(limits=c(-8,10),
                       breaks = c(-7,0,9),low="#0000C7", high="#FF000C", mid="white",midpoint = 0)+
  guides(fill=guide_colorbar(title = "log2FC",title.hjust = 0,title.vjust = 0))+
  theme(axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "gray90"),
        #panel.border = element_rect(colour = "black"),
        axis.text.x = element_text(size = 7,colour='black'),
        #axis.text.y = element_markdown(size = 7),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        #legend.key.size = unit(0.6, "cm"),  #key size
        legend.key.height = unit(0.2, "cm"),
        legend.key.width = unit(0.3, "cm"),
        legend.text = element_text(size = 7, color = "black"), #legend label
        legend.spacing.x = unit(0.2, "cm"),#distance between key and label)
        legend.position = "right",
        legend.direction = "horizontal",
        legend.margin = margin(0,0,0,0),
        legend.title = element_text(size = 7,color = "black"),
        legend.title.position = "top",
        plot.margin = margin(0,0,0,0),
        strip.text = element_text(size = 8,colour = "black"),
        strip.background = element_rect(linetype = "solid",
                                        color = "white", fill = "gray90",
                                        linewidth = 0.5),
        panel.spacing.x= unit(0.15, "lines")
  ) +
  geom_text(aes(label = sig),size=2)


# Family legend part
df$x_name <-"x_name"

p2 <- ggplot(df,aes(x=x_name,y=font_name))+
  geom_tile(aes(fill=font_Family),color="black")+
  scale_fill_manual(values = my_color2)+
  scale_y_discrete(position = "right",
                   sec.axis = dup_axis(name = waiver(), labels = waiver()))+
  labs(x="",y="",fill="Family")+
  theme(
        axis.text.y.right = element_markdown(margin = margin(l = 2),size = 7,hjust = 0,vjust = 0.5),
        axis.text.y.left = element_blank(),
        axis.text.x =element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),,
        legend.key.size = unit(0.8, "cm"),  #key size
        legend.key.height = unit(0.2, "cm"),
        legend.key.width = unit(0.3, "cm"),
        legend.text = element_markdown(size = 7, color = "black",margin = margin(l=4)), #legend label
        legend.spacing.x = unit(1, "cm"),
        legend.position = "right",
        legend.direction = "vertical",
        legend.title = element_text(size = 7,color = "black"),
        legend.margin = margin(0,0,0,0),
        legend.key = element_rect(color = NA),
        plot.margin = margin(0,0,0,0),
  )


library(aplot)
p3 <- p1 %>%  insert_right(p2, width = 0.02)
p3

png("deg/ancombc_species.png", width = 7, height = 6,res = 600, unit="in")
p3
dev.off()



pdf("deg/ancombc_species.pdf", width = 7, height = 6)
p3
dev.off()


############################################
#####                                  #####
#####            STOP HERE             #####
#####                                  #####
############################################



# Generate table--------------
## For ASVs--------------------
library(phyloseq)
library(ANCOMBC)
library(dplyr)
f_table <- "feature_table_filter_ileum.csv"
s_table <- "sample_data_ileum.csv"
t_table <- "taxonomy.csv"
feature_table <- read.csv(f_table,header = T,row.names = 1)
taxonomy <- read.csv(t_table,header = T,row.names = 1)
sample_data <- read.csv(s_table,header = T,row.names = 1)
#Transform feature table and taxonomy table to matrix
sample_data$day <- factor(sample_data$day,levels = c("3 dpi","5 dpi","7 dpi","14 dpi","21 dpi"))
sample_data$treatment <- factor(sample_data$treatment,levels = c("Mock","EM"))
OTU<- otu_table(as.matrix(feature_table),taxa_are_rows = T)
TAX<- tax_table(as.matrix(taxonomy))
SAM<-sample_data(sample_data)
physeq<- phyloseq(OTU,TAX,SAM)

phy_3 <- subset_samples(physeq,day=="3 dpi")
phy_5 <- subset_samples(physeq,day=="5 dpi")
phy_7 <- subset_samples(physeq,day=="7 dpi")
phy_14 <- subset_samples(physeq,day=="14 dpi")
phy_21 <- subset_samples(physeq,day=="21 dpi")


my_daa <- function(data,tax_level,filename){
  set.seed(10022024)
  output <- ancombc(data=data,assay_name =NULL,
                    formula = "treatment",group = "treatment",tax_level=tax_level,
                    prv_cut = 0.0,lib_cut = 0, struc_zero = T,neg_lb = F,alpha=0.05,
                    n_cl = 3,verbose = T)
  #primary results
  res_prim = output$res
  #export
  write.csv(res_prim,file = filename)
}

my_daa(phy_3,"Species","res_prim_3.csv")
my_daa(phy_5,"Species","res_prim_5.csv")
my_daa(phy_7,"Species","res_prim_7.csv")
my_daa(phy_14,"Species","res_prim_14.csv")
my_daa(phy_21,"Species","res_prim_21.csv")



#combine results
r3 <- read.csv("res_prim_3.csv",header = T,row.names = 1)
r5 <- read.csv("res_prim_5.csv",header = T,row.names = 1)
r7 <- read.csv("res_prim_7.csv",header = T,row.names = 1)
r14 <- read.csv("res_prim_14.csv",header = T,row.names = 1)
r21 <- read.csv("res_prim_21.csv",header = T,row.names = 1)

r3$zero <- ifelse(r3$se.treatmentEM == 0 & r3$q_val.treatmentEM ==0, "szero","notzero")
r5$zero <- ifelse(r5$se.treatmentEM == 0 & r5$q_val.treatmentEM ==0, "szero","notzero")
r7$zero <- ifelse(r7$se.treatmentEM == 0 & r7$q_val.treatmentEM ==0, "szero","notzero")
r14$zero <- ifelse(r14$se.treatmentEM == 0 & r14$q_val.treatmentEM ==0, "szero","notzero")
r21$zero <- ifelse(r21$se.treatmentEM == 0 & r21$q_val.treatmentEM ==0, "szero","notzero")


#select top 50 ASV
f_table <- "feature_table_relative.csv"
s_table <- "sample_data_ileum.csv"
t_table <- "taxonomy.csv"
otu <-read.csv(f_table,header = T,row.names = 1)
sample_data <-read.csv(s_table,header = T,row.names = 1)
tax <- read.csv(t_table,header = T,row.names = 1)
otu$rownames <- rownames(otu)
tax$rownames <- rownames(tax)

number=50
feature_otu <- merge(tax, otu, by='rownames')
rownames(feature_otu) <-feature_otu$Species
feature_otu <- feature_otu[,-c(1:9)]
feature_otu <- feature_otu[order(rowSums(feature_otu),decreasing = T),]
feature_otu <- feature_otu[1:number,]

#calculate mean and sd
feature_otu <- as.data.frame(t(feature_otu))
identical(rownames(feature_otu),rownames(sample_data))

feature_otu <- feature_otu*100
#v_mean <- aggregate(feature_otu,by=list(trt=sample_data$trt),mean)
#v_sd <- aggregate(feature_otu,by=list(trt=sample_data$trt),sd)
sample_data$trt <- paste(sample_data$day,sample_data$treatment,sep = "_")
v_median <- aggregate(feature_otu,by=list(trt=sample_data$trt),mean)
#modify format of mean and sd
rownames(v_median) <- v_median$trt
v_median <- v_median[,-1]
v_median <- as.data.frame(t(v_median))
v_median <- apply(v_median,c(1,2),FUN=function(x){
  sprintf("%0.2f",x)
})

#rownames(v_mean) <- v_mean$trt
#v_mean <- v_mean[,-1]
#v_mean <- as.data.frame(t(v_mean))
#v_mean <- apply(v_mean,c(1,2),FUN=function(x){
#  sprintf("%0.3f",x)
#})

#rownames(v_sd) <- v_sd$trt
#v_sd <- v_sd[,-1]
#v_sd <- as.data.frame(t(v_sd))
#v_sd <- apply(v_sd,c(1,2),FUN=function(x){
#  sprintf("%0.2f",x)
#})

#combine mean and sd
#identical(rownames(v_mean),rownames(v_sd))

#result_table <- data.frame(matrix(ncol = ncol(v_mean), nrow = nrow(v_mean)))
#rownames(result_table) <- rownames(v_mean)
#colnames(result_table) <- colnames(v_mean)

#for (i in 1:nrow(result_table)){
#  for (j in 1:ncol(result_table)){
#    result_table[i,j] <- paste(v_mean[i,j],v_sd[i,j],sep = " ± ")
#  }
#}

result_table <- as.data.frame(v_median)
result_table <- result_table[,c(6,5,8,7,10,9,2,1,4,3)]

#add logFC info
my_info <- function(data){
  df <- data
  info <- data.frame(row.names = rownames(result_table))
  info$lfc <- df$lfc.treatmentEM[match(rownames(info),df$lfc.taxon)]
  info$qval <- df$q_val.treatmentEM[match(rownames(info),df$lfc.taxon)]
  info$zero <- df$zero[match(rownames(info),df$lfc.taxon)]
  identical(rownames(info),rownames(result_table))
  
  info$sig <-  case_when(
    info$zero == "szero" ~ "***",
    info$qval >= 0.05 ~ "",
    info$qval >= 0.01 & info$qval < 0.05 ~ "*",
    info$qval >= 0.001 & info$qval < 0.01 ~ "**",
    info$qval < 0.001 ~ "***"
  )
  return(info)
}



my_lfc <- function(fcname){
  result_table[[fcname]] <-  sprintf("%0.2f",info$lfc)
  for (i in 1:nrow(result_table)){
    if (info$sig[i] != ""){
      result_table[[fcname]][i] <- paste(result_table[[fcname]][i],info$sig[i],sep = "")
    }
  }
  return(result_table)
}

#3 dpi
info <- my_info(r3)
result_table <- my_lfc("lfc3")

#5 dpi
info <- my_info(r5)
result_table <- my_lfc("lfc5")

#7 dpi
info <- my_info(r7)
result_table <- my_lfc("lfc7")

#14 dpi
info <- my_info(r14)
result_table <- my_lfc("lfc14")

#21 dpi
info <- my_info(r21)
result_table <- my_lfc("lfc21")

#export
library(openxlsx)
write.xlsx(result_table, file = "diff_table_ASV.xlsx",rowNames = T)


# For Genus-----------------
library(phyloseq)
library(ANCOMBC)
library(dplyr)
f_table <- "feature_table_filter_ileum.csv"
s_table <- "sample_data_ileum.csv"
t_table <- "taxonomy.csv"
feature_table <- read.csv(f_table,header = T,row.names = 1)
taxonomy <- read.csv(t_table,header = T,row.names = 1)
taxonomy$Order <- gsub("Eubacteriales","Clostridiales",taxonomy$Order)
taxonomy$Class <- gsub("Bacillota_incertae_sedis","Clostridia",taxonomy$Class)
taxonomy$Order <- gsub("Bacillota_incertae_sedis","Clostridiales",taxonomy$Order)
taxonomy$Family <- gsub("Bacillota_incertae_sedis","Oscillospiraceae",taxonomy$Family)
taxonomy$Class <- gsub("Actinobacteria","Actinomycetia",taxonomy$Class)
taxonomy$Order <- gsub("Bacillales","Caryophanales",taxonomy$Order)
taxonomy$Order <- gsub("Dermabacterales","Micrococcales",taxonomy$Order)
taxonomy$Family <- gsub("Mycobacteriaceae","Corynebacteriaceae",taxonomy$Family)

sample_data <- read.csv(s_table,header = T,row.names = 1)
#Transform feature table and taxonomy table to matrix
sample_data$day <- factor(sample_data$day,levels = c("3 dpi","5 dpi","7 dpi","14 dpi","21 dpi"))
sample_data$treatment <- factor(sample_data$treatment,levels = c("Mock","EM"))
OTU<- otu_table(as.matrix(feature_table),taxa_are_rows = T)
TAX<- tax_table(as.matrix(taxonomy))
SAM<-sample_data(sample_data)
physeq<- phyloseq(OTU,TAX,SAM)

phy_3 <- subset_samples(physeq,day=="3 dpi")
phy_5 <- subset_samples(physeq,day=="5 dpi")
phy_7 <- subset_samples(physeq,day=="7 dpi")
phy_14 <- subset_samples(physeq,day=="14 dpi")
phy_21 <- subset_samples(physeq,day=="21 dpi")


my_daa <- function(data,tax_level,filename){
  set.seed(10022024)
  output <- ancombc(data=data,assay_name =NULL,
                    formula = "treatment",group = "treatment",tax_level=tax_level,
                    prv_cut = 0.0,lib_cut = 0, struc_zero = T,neg_lb = F,alpha=0.05,
                    n_cl = 3,verbose = T)
  #primary results
  res_prim = output$res
  #export
  write.csv(res_prim,file = filename)
}

my_daa(phy_3,"Genus","res_prim_3.csv")
my_daa(phy_5,"Genus","res_prim_5.csv")
my_daa(phy_7,"Genus","res_prim_7.csv")
my_daa(phy_14,"Genus","res_prim_14.csv")
my_daa(phy_21,"Genus","res_prim_21.csv")



#combine results
r3 <- read.csv("res_prim_3.csv",header = T,row.names = 1)
r5 <- read.csv("res_prim_5.csv",header = T,row.names = 1)
r7 <- read.csv("res_prim_7.csv",header = T,row.names = 1)
r14 <- read.csv("res_prim_14.csv",header = T,row.names = 1)
r21 <- read.csv("res_prim_21.csv",header = T,row.names = 1)

r3$zero <- ifelse(r3$se.treatmentEM == 0 & r3$q_val.treatmentEM ==0, "szero","notzero")
r5$zero <- ifelse(r5$se.treatmentEM == 0 & r5$q_val.treatmentEM ==0, "szero","notzero")
r7$zero <- ifelse(r7$se.treatmentEM == 0 & r7$q_val.treatmentEM ==0, "szero","notzero")
r14$zero <- ifelse(r14$se.treatmentEM == 0 & r14$q_val.treatmentEM ==0, "szero","notzero")
r21$zero <- ifelse(r21$se.treatmentEM == 0 & r21$q_val.treatmentEM ==0, "szero","notzero")


#select top 25 genus
f_table <- "feature_table_relative.csv"
s_table <- "sample_data_ileum.csv"
t_table <- "taxonomy.csv"
otu <-read.csv(f_table,header = T,row.names = 1)
sample_data <-read.csv(s_table,header = T,row.names = 1)
tax <- read.csv(t_table,header = T,row.names = 1)
otu$rownames <- rownames(otu)
tax$rownames <- rownames(tax)

number=25
feature_otu <- merge(tax, otu, by='rownames')
feature_otu <- feature_otu[,-c(1:6,8,9)]
feature_otu <- aggregate(feature_otu[,-1],by=list(Genus=feature_otu$Genus),FUN = "sum")
rownames(feature_otu) <- feature_otu$Genus
feature_otu <- feature_otu[,-1]
feature_otu <- feature_otu[order(rowSums(feature_otu),decreasing = T),]
feature_otu <- feature_otu[1:number,]

#calculate mean and sd
feature_otu <- as.data.frame(t(feature_otu))
identical(rownames(feature_otu),rownames(sample_data))

feature_otu <- feature_otu*100
sample_data$trt <- paste(sample_data$day,sample_data$treatment,sep = "_")
v_median <- aggregate(feature_otu,by=list(trt=sample_data$trt),mean)
#modify format of mean and sd
rownames(v_median) <- v_median$trt
v_median <- v_median[,-1]
v_median <- as.data.frame(t(v_median))
v_median <- apply(v_median,c(1,2),FUN=function(x){
  sprintf("%0.2f",x)
})

result_table <- as.data.frame(v_median)
result_table <- result_table[,c(6,5,8,7,10,9,2,1,4,3)]

#add logFC info
my_info <- function(data){
  df <- data
  info <- data.frame(row.names = rownames(result_table))
  info$lfc <- df$lfc.treatmentEM[match(rownames(info),df$lfc.taxon)]
  info$qval <- df$q_val.treatmentEM[match(rownames(info),df$lfc.taxon)]
  info$zero <- df$zero[match(rownames(info),df$lfc.taxon)]
  identical(rownames(info),rownames(result_table))
  
  info$sig <-  case_when(
    info$zero == "szero" ~ "***",
    info$qval >= 0.05 ~ "",
    info$qval >= 0.01 & info$qval < 0.05 ~ "*",
    info$qval >= 0.001 & info$qval < 0.01 ~ "**",
    info$qval < 0.001 ~ "***"
  )
  return(info)
}



my_lfc <- function(fcname){
  result_table[[fcname]] <-  sprintf("%0.2f",info$lfc)
  for (i in 1:nrow(result_table)){
    if (info$sig[i] != ""){
      result_table[[fcname]][i] <- paste(result_table[[fcname]][i],info$sig[i],sep = "")
    }
  }
  return(result_table)
}

#3 dpi
info <- my_info(r3)
result_table <- my_lfc("lfc3")

#5 dpi
info <- my_info(r5)
result_table <- my_lfc("lfc5")

#7 dpi
info <- my_info(r7)
result_table <- my_lfc("lfc7")

#14 dpi
info <- my_info(r14)
result_table <- my_lfc("lfc14")

#21 dpi
info <- my_info(r21)
result_table <- my_lfc("lfc21")

#export
library(openxlsx)
write.xlsx(result_table, file = "diff_table_genus.xlsx",rowNames = T)


# For Family-----------------
library(phyloseq)
library(ANCOMBC)
library(dplyr)
f_table <- "feature_table_filter_ileum.csv"
s_table <- "sample_data_ileum.csv"
t_table <- "taxonomy.csv"
feature_table <- read.csv(f_table,header = T,row.names = 1)
taxonomy <- read.csv(t_table,header = T,row.names = 1)

taxonomy$Order <- gsub("Eubacteriales","Clostridiales",taxonomy$Order)
taxonomy$Class <- gsub("Bacillota_incertae_sedis","Clostridia",taxonomy$Class)
taxonomy$Order <- gsub("Bacillota_incertae_sedis","Clostridiales",taxonomy$Order)
taxonomy$Class <- gsub("Actinobacteria","Actinomycetia",taxonomy$Class)
taxonomy$Order <- gsub("Bacillales","Caryophanales",taxonomy$Order)
taxonomy$Order <- gsub("Dermabacterales","Micrococcales",taxonomy$Order)
taxonomy$Order <- gsub("Hyphomicrobiales","Rhizobiales",taxonomy$Order)

sample_data <- read.csv(s_table,header = T,row.names = 1)
#Transform feature table and taxonomy table to matrix
sample_data$day <- factor(sample_data$day,levels = c("3 dpi","5 dpi","7 dpi","14 dpi","21 dpi"))
sample_data$treatment <- factor(sample_data$treatment,levels = c("Mock","EM"))
OTU<- otu_table(as.matrix(feature_table),taxa_are_rows = T)
TAX<- tax_table(as.matrix(taxonomy))
SAM<-sample_data(sample_data)
physeq<- phyloseq(OTU,TAX,SAM)

phy_3 <- subset_samples(physeq,day=="3 dpi")
phy_5 <- subset_samples(physeq,day=="5 dpi")
phy_7 <- subset_samples(physeq,day=="7 dpi")
phy_14 <- subset_samples(physeq,day=="14 dpi")
phy_21 <- subset_samples(physeq,day=="21 dpi")



my_daa <- function(data,tax_level,filename){
  set.seed(10022024)
  output <- ancombc(data=data,assay_name =NULL,
                    formula = "treatment",group = "treatment",tax_level=tax_level,
                    prv_cut = 0.0,lib_cut = 0, struc_zero = T,neg_lb = F,alpha=0.05,
                    n_cl = 3,verbose = T)
  #primary results
  res_prim = output$res
  #export
  write.csv(res_prim,file = filename)
}

my_daa(phy_3,"Family","res_prim_3.csv")
my_daa(phy_5,"Family","res_prim_5.csv")
my_daa(phy_7,"Family","res_prim_7.csv")
my_daa(phy_14,"Family","res_prim_14.csv")
my_daa(phy_21,"Family","res_prim_21.csv")


#combine results
r3 <- read.csv("res_prim_3.csv",header = T,row.names = 1)
r5 <- read.csv("res_prim_5.csv",header = T,row.names = 1)
r7 <- read.csv("res_prim_7.csv",header = T,row.names = 1)
r14 <- read.csv("res_prim_14.csv",header = T,row.names = 1)
r21 <- read.csv("res_prim_21.csv",header = T,row.names = 1)

r3$zero <- ifelse(r3$se.treatmentEM == 0 & r3$q_val.treatmentEM ==0, "szero","notzero")
r5$zero <- ifelse(r5$se.treatmentEM == 0 & r5$q_val.treatmentEM ==0, "szero","notzero")
r7$zero <- ifelse(r7$se.treatmentEM == 0 & r7$q_val.treatmentEM ==0, "szero","notzero")
r14$zero <- ifelse(r14$se.treatmentEM == 0 & r14$q_val.treatmentEM ==0, "szero","notzero")
r21$zero <- ifelse(r21$se.treatmentEM == 0 & r21$q_val.treatmentEM ==0, "szero","notzero")


#select top 10 Family
f_table <- "feature_table_relative.csv"
s_table <- "sample_data_ileum.csv"
t_table <- "taxonomy.csv"
otu <-read.csv(f_table,header = T,row.names = 1)
sample_data <-read.csv(s_table,header = T,row.names = 1)
tax <- read.csv(t_table,header = T,row.names = 1)
otu$rownames <- rownames(otu)
tax$rownames <- rownames(tax)

number=10
feature_otu <- merge(tax, otu, by='rownames')
feature_otu <- feature_otu[,-c(1:5,7:9)]
feature_otu <- aggregate(feature_otu[,-1],by=list(Genus=feature_otu$Family),FUN = "sum")
rownames(feature_otu) <- feature_otu$Genus
feature_otu <- feature_otu[,-1]
feature_otu <- feature_otu[order(rowSums(feature_otu),decreasing = T),]
feature_otu <- feature_otu[1:number,]

#calculate mean and sd
feature_otu <- as.data.frame(t(feature_otu))
identical(rownames(feature_otu),rownames(sample_data))

feature_otu <- feature_otu*100
sample_data$trt <- paste(sample_data$day,sample_data$treatment,sep = "_")
v_median <- aggregate(feature_otu,by=list(trt=sample_data$trt),mean)
#modify format of mean and sd
rownames(v_median) <- v_median$trt
v_median <- v_median[,-1]
v_median <- as.data.frame(t(v_median))
v_median <- apply(v_median,c(1,2),FUN=function(x){
  sprintf("%0.2f",x)
})

result_table <- as.data.frame(v_median)
result_table <- result_table[,c(6,5,8,7,10,9,2,1,4,3)]

#add logFC info
my_info <- function(data){
  df <- data
  info <- data.frame(row.names = rownames(result_table))
  info$lfc <- df$lfc.treatmentEM[match(rownames(info),df$lfc.taxon)]
  info$qval <- df$q_val.treatmentEM[match(rownames(info),df$lfc.taxon)]
  info$zero <- df$zero[match(rownames(info),df$lfc.taxon)]
  identical(rownames(info),rownames(result_table))
  
  info$sig <-  case_when(
    info$zero == "szero" ~ "***",
    info$qval >= 0.05 ~ "",
    info$qval >= 0.01 & info$qval < 0.05 ~ "*",
    info$qval >= 0.001 & info$qval < 0.01 ~ "**",
    info$qval < 0.001 ~ "***"
  )
  return(info)
}



my_lfc <- function(fcname){
  result_table[[fcname]] <-  sprintf("%0.2f",info$lfc)
  for (i in 1:nrow(result_table)){
    if (info$sig[i] != ""){
      result_table[[fcname]][i] <- paste(result_table[[fcname]][i],info$sig[i],sep = "")
    }
  }
  return(result_table)
}

#3 dpi
info <- my_info(r3)
result_table <- my_lfc("lfc3")

#5 dpi
info <- my_info(r5)
result_table <- my_lfc("lfc5")

#7 dpi
info <- my_info(r7)
result_table <- my_lfc("lfc7")

#14 dpi
info <- my_info(r14)
result_table <- my_lfc("lfc14")

#21 dpi
info <- my_info(r21)
result_table <- my_lfc("lfc21")

#export
library(openxlsx)
write.xlsx(result_table, file = "diff_table_family.xlsx",rowNames = T)


# For Phylum-----------------
library(phyloseq)
library(ANCOMBC)
library(dplyr)
f_table <- "feature_table_filter_ileum.csv"
s_table <- "sample_data_ileum.csv"
t_table <- "taxonomy.csv"
feature_table <- read.csv(f_table,header = T,row.names = 1)
taxonomy <- read.csv(t_table,header = T,row.names = 1)
sample_data <- read.csv(s_table,header = T,row.names = 1)
#Transform feature table and taxonomy table to matrix
sample_data$day <- factor(sample_data$day,levels = c("3 dpi","5 dpi","7 dpi","14 dpi","21 dpi"))
sample_data$treatment <- factor(sample_data$treatment,levels = c("Mock","EM"))
OTU<- otu_table(as.matrix(feature_table),taxa_are_rows = T)
TAX<- tax_table(as.matrix(taxonomy))
SAM<-sample_data(sample_data)
physeq<- phyloseq(OTU,TAX,SAM)

phy_3 <- subset_samples(physeq,day=="3 dpi")
phy_5 <- subset_samples(physeq,day=="5 dpi")
phy_7 <- subset_samples(physeq,day=="7 dpi")
phy_14 <- subset_samples(physeq,day=="14 dpi")
phy_21 <- subset_samples(physeq,day=="21 dpi")


my_daa <- function(data,tax_level,filename){
  set.seed(10022024)
  output <- ancombc(data=data,assay_name =NULL,
                    formula = "treatment",group = "treatment",tax_level=tax_level,
                    prv_cut = 0.0,lib_cut = 0, struc_zero = T,neg_lb = F,alpha=0.05,
                    n_cl = 3,verbose = T)
  #primary results
  res_prim = output$res
  #export
  write.csv(res_prim,file = filename)
}

my_daa(phy_3,"Phylum","res_prim_3.csv")
my_daa(phy_5,"Phylum","res_prim_5.csv")
my_daa(phy_7,"Phylum","res_prim_7.csv")
my_daa(phy_14,"Phylum","res_prim_14.csv")
my_daa(phy_21,"Phylum","res_prim_21.csv")



#combine results
r3 <- read.csv("res_prim_3.csv",header = T,row.names = 1)
r5 <- read.csv("res_prim_5.csv",header = T,row.names = 1)
r7 <- read.csv("res_prim_7.csv",header = T,row.names = 1)
r14 <- read.csv("res_prim_14.csv",header = T,row.names = 1)
r21 <- read.csv("res_prim_21.csv",header = T,row.names = 1)

r3$zero <- ifelse(r3$se.treatmentEM == 0 & r3$q_val.treatmentEM ==0, "szero","notzero")
r5$zero <- ifelse(r5$se.treatmentEM == 0 & r5$q_val.treatmentEM ==0, "szero","notzero")
r7$zero <- ifelse(r7$se.treatmentEM == 0 & r7$q_val.treatmentEM ==0, "szero","notzero")
r14$zero <- ifelse(r14$se.treatmentEM == 0 & r14$q_val.treatmentEM ==0, "szero","notzero")
r21$zero <- ifelse(r21$se.treatmentEM == 0 & r21$q_val.treatmentEM ==0, "szero","notzero")


#select top phylum
f_table <- "feature_table_relative.csv"
s_table <- "sample_data_ileum.csv"
t_table <- "taxonomy.csv"
otu <-read.csv(f_table,header = T,row.names = 1)
sample_data <-read.csv(s_table,header = T,row.names = 1)
tax <- read.csv(t_table,header = T,row.names = 1)
otu$rownames <- rownames(otu)
tax$rownames <- rownames(tax)

number=3  #cecum 3
feature_otu <- merge(tax, otu, by='rownames')
feature_otu <- feature_otu[,-c(1:2,4:9)]
feature_otu <- aggregate(feature_otu[,-1],by=list(Genus=feature_otu$Phylum),FUN = "sum")
rownames(feature_otu) <- feature_otu$Genus
feature_otu <- feature_otu[,-1]
feature_otu <- feature_otu[order(rowSums(feature_otu),decreasing = T),]
feature_otu <- feature_otu[1:number,]

#calculate mean and sd
feature_otu <- as.data.frame(t(feature_otu))
identical(rownames(feature_otu),rownames(sample_data))

feature_otu <- feature_otu*100
sample_data$trt <- paste(sample_data$day,sample_data$treatment,sep = "_")
v_median <- aggregate(feature_otu,by=list(trt=sample_data$trt),mean)
#modify format of mean and sd
rownames(v_median) <- v_median$trt
v_median <- v_median[,-1]
v_median <- as.data.frame(t(v_median))
v_median <- apply(v_median,c(1,2),FUN=function(x){
  sprintf("%0.2f",x)
})

result_table <- as.data.frame(v_median)
result_table <- result_table[,c(6,5,8,7,10,9,2,1,4,3)]

#add logFC info
my_info <- function(data){
  df <- data
  info <- data.frame(row.names = rownames(result_table))
  info$lfc <- df$lfc.treatmentEM[match(rownames(info),df$lfc.taxon)]
  info$qval <- df$q_val.treatmentEM[match(rownames(info),df$lfc.taxon)]
  info$zero <- df$zero[match(rownames(info),df$lfc.taxon)]
  identical(rownames(info),rownames(result_table))
  
  info$sig <-  case_when(
    info$zero == "szero" ~ "***",
    info$qval >= 0.05 ~ "",
    info$qval >= 0.01 & info$qval < 0.05 ~ "*",
    info$qval >= 0.001 & info$qval < 0.01 ~ "**",
    info$qval < 0.001 ~ "***"
  )
  return(info)
}



my_lfc <- function(fcname){
  result_table[[fcname]] <-  sprintf("%0.2f",info$lfc)
  for (i in 1:nrow(result_table)){
    if (info$sig[i] != ""){
      result_table[[fcname]][i] <- paste(result_table[[fcname]][i],info$sig[i],sep = "")
    }
  }
  return(result_table)
}

#3 dpi
info <- my_info(r3)
result_table <- my_lfc("lfc3")

#5 dpi
info <- my_info(r5)
result_table <- my_lfc("lfc5")

#7 dpi
info <- my_info(r7)
result_table <- my_lfc("lfc7")

#14 dpi
info <- my_info(r14)
result_table <- my_lfc("lfc14")

#21 dpi
info <- my_info(r21)
result_table <- my_lfc("lfc21")

#export
library(openxlsx)
write.xlsx(result_table, file = "diff_table_phylum.xlsx",rowNames = T)
