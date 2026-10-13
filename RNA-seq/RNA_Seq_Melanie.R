setwd("D:\\OneDrive - Oklahoma A and M System\\PhD\\Project\\Melanie\\Transcriptome\\result")

setwd("C:\\Users\\jiaguo\\OneDrive - Oklahoma A and M System\\PhD\\Project\\Melanie\\Transcriptome\\result3")

#biomaRt

#modify raw counts file------------------
cts <- read.table("counts.txt",header = T,row.names=1,sep = "\t", skip = 1,check.names = F)
cts <- cts[,-c(1:5)]
colnames(cts) <- sub("/scratch/josh/Melanie/DCA_RNA_Seq/result/allign/allign_","",colnames(cts))
colnames(cts) <- sub(".sam","",colnames(cts))
write.csv(cts,file = "counts.csv")


#DEG analysis--------------------
library(DESeq2)

#import data as DESeq2 object
cts <- read.csv("counts.csv",header = T,row.names = 1)
sample_data <- read.csv("sample_data.csv",header = T,row.names = 1)
sample_data$group <- factor(sample_data$group)
dds <- DESeqDataSetFromMatrix(countData=cts, colData=sample_data, design= ~ group)

#filter genes 
#(keep only rows that have a count of at least 10 for a minimal number of samples(group size))
keep <- rowSums(counts(dds) >= 10) >= 6
dds <- dds[keep,]

#order levels
dds$group <- factor(dds$group,levels = c("Mock","DCA","NE","DCA_NE"))


#DEG analysis
dds2 <- DESeq(dds)

#check results
summary(results(dds2))


#get normalized counts
nor_counts <- as.data.frame(counts(dds2,normalized=T))
write.csv(nor_counts,file = "normalized_counts.csv")
#export normalized counts
library(dplyr)

for (i in 1:nrow(nor_counts)){
  nor_counts$mean_DCA[i] <- mean(as.numeric(nor_counts[i,1:6]))
  nor_counts$mean_Mock[i] <- mean(as.numeric(nor_counts[i,7:12]))
  nor_counts$mean_NE[i] <- mean(as.numeric(nor_counts[i,13:18]))
  nor_counts$"mean_DCA+NE"[i] <- mean(as.numeric(nor_counts[i,19:24]))
}

nor_counts$FC_DCAvsMock <- nor_counts$mean_DCA/nor_counts$mean_Mock 
nor_counts$FC_NEvsMock <- nor_counts$mean_NE/nor_counts$mean_Mock 
nor_counts$"FC_DCA+NEvsMock" <- nor_counts$"mean_DCA+NE"/nor_counts$mean_Mock  
nor_counts$"FC_DCA+NEvsNE" <- nor_counts$"mean_DCA+NE"/nor_counts$mean_NE

nor_counts$"manual_log2FC_DCAvsMock" <- log2(nor_counts$"FC_DCAvsMock")
nor_counts$"manual_log2FC_NEvsMock" <- log2(nor_counts$"FC_NEvsMock") 
nor_counts$"manual_log2FC_DCA+NEvsMock" <- log2(nor_counts$"FC_DCA+NEvsMock") 
nor_counts$"manual_log2FC_DCA+NEvsNE" <- log2(nor_counts$"FC_DCA+NEvsNE") 

write.csv(nor_counts,file = "2222.csv")

#get DEGs results--------------------
res_DCA_Mock <- results(dds2,contrast=c("group","DCA","Mock"))
res_NE_Mock <- results(dds2,contrast=c("group","NE","Mock"))
res_DCA_NE_Mock <- results(dds2,contrast=c("group","DCA_NE","Mock"))
res_DCA_NE_NE <- results(dds2,contrast=c("group","DCA_NE","NE"))

#export selected results
DCA_Mock <- as.data.frame(res_DCA_Mock)
NE_Mock <- as.data.frame(res_NE_Mock)
DCA_NE_Mock <- as.data.frame(res_DCA_NE_Mock)
DCA_NE_NE <- as.data.frame(res_DCA_NE_NE)
DCA_Mock[is.na(DCA_Mock)] <- 1
NE_Mock[is.na(NE_Mock)] <- 1
DCA_NE_Mock[is.na(DCA_NE_Mock)] <- 1
DCA_NE_NE[is.na(DCA_NE_NE)] <- 1

nor_counts$"p.adjust_DCAvsMock" <- DCA_Mock$padj
nor_counts$"DESeq2_log2FC_DCAvsMock" <- DCA_Mock$log2FoldChange

nor_counts$"p.adjust_NEvsMock" <- NE_Mock$padj
nor_counts$"DESeq2_log2FC_NEvsMock" <- NE_Mock$log2FoldChange

nor_counts$"p.adjust_DCA+NEvsDCA" <- DCA_NE_Mock$padj
nor_counts$"DESeq2_log2FC_DCA+NEvsMock" <- DCA_NE_Mock$log2FoldChange

nor_counts$"p.adjust_DCA+NEvsNE" <- DCA_NE_NE$padj
nor_counts$"DESeq2_log2FC_DCA+NEvsNE" <- DCA_NE_NE$log2FoldChange

write.csv(nor_counts,file = "normalized_counts_detail.csv")

write.csv(as.data.frame(res_DCA_Mock),file = "DCAvsMock.csv")
write.csv(as.data.frame(res_NE_Mock),file = "NEvsMock.csv")
write.csv(as.data.frame(res_DCA_NE_Mock),file = "DCA+NEvsMock.csv")
write.csv(as.data.frame(res_DCA_NE_NE),file = "DCA+NEvsNE.csv")



#select DEGs
res_DCA_Mock_sub <- subset(res_DCA_Mock, padj <= 0.001&abs(log2FoldChange) >= 1)
res_NE_Mock_sub <- subset(res_NE_Mock, padj <= 0.001&abs(log2FoldChange) >= 1)
res_DCA_NE_Mock_sub <- subset(res_DCA_NE_Mock, padj <= 0.001&abs(log2FoldChange) >= 1)
res_DCA_NE_NE_sub <- subset(res_DCA_NE_NE, padj <= 0.001&abs(log2FoldChange) >= 1)

res_DCA_Mock_sub@nrows
res_NE_Mock_sub@nrows
res_DCA_NE_Mock_sub@nrows
res_DCA_NE_NE_sub@nrows 


res_DCA_Mock_sub <- res_DCA_Mock_sub[,c(2,6)]
res_NE_Mock_sub <- res_NE_Mock_sub[,c(2,6)]
res_DCA_NE_Mock_sub <- res_DCA_NE_Mock_sub[,c(2,6)]

#Export DEGs
write.csv(res_DCA_Mock_sub,file = "res_DCA_Mock001.csv")
write.csv(res_NE_Mock_sub,file = "res_NE_Mock001.csv")
write.csv(res_DCA_NE_Mock_sub,file = "res_DCA_NE_Mock001.csv")

# Combine DEGs
library(tidyverse)
DCA_Mock <- read.csv("res_DCA_Mock.csv")
NE_Mock <- read.csv("res_NE_Mock.csv")
DCA_NE_Mock <- read.csv("res_DCA_NE_Mock.csv")

df2 <- DCA_Mock %>% 
  full_join(.,NE_Mock,by="X") %>% 
  full_join(.,DCA_NE_Mock,by="X")


df_final <- full_join(df,df2,by="X")
write.csv(df_final,"DEG.csv")


#export DEGS number of each comparsion
up_DCA_Mock <- sum(res_DCA_Mock_sub$log2FoldChange > 0)
down_DCA_Mock <- sum(res_DCA_Mock_sub$log2FoldChange < 0)
up_NE_Mock <- sum(res_NE_Mock_sub$log2FoldChange > 0)
down_NE_Mock <- sum(res_NE_Mock_sub$log2FoldChange < 0)
up_DCA_NE_Mock <- sum(res_DCA_NE_Mock_sub$log2FoldChange > 0)
down_DCA_NE_Mock <- sum(res_DCA_NE_Mock_sub$log2FoldChange < 0)


deg_num <- as.data.frame(rbind(up_DCA_Mock,down_DCA_Mock,up_NE_Mock,down_NE_Mock,
      up_DCA_NE_Mock,down_DCA_NE_Mock))

deg_num$group <- c(rep("DCA_Mock",2),rep( "NE_Mock",2),rep("DCA_NE_Mock",2))
deg_num$Regulation <- c(rep(c("Up","Down"),3))
colnames(deg_num) <-c("number","group","Regulation")
write.csv(deg_num,file = "DEGs_number001.csv")

#PCA analysis---------------------------------------
library(vegan)
library(tidyverse)

#data import
counts_norm <- read.csv("normalized_counts.csv",header = T,row.names = 1,check.names = F)
sample_data <- read.csv("sample_data.csv",header = T,row.names = 1,check.names = F)

#calculate dist matrix
dist <- vegdist(t(counts_norm),method = "bray")

#permanova analysis
set.seed(111)
permanova <- adonis2(dist ~ group,data=sample_data,permutations = 10000)

#Export statistical analysis result (R square and P value)
permanova$source <- 'Bray_Curtis'
write.csv(permanova,file = "permanova.csv")

#pcoa analysis
pcoa <- cmdscale(dist,k=2,eig = T,add = T)

pcoa_info <- as.data.frame(cbind(pcoa$points,pcoa$eig))
colnames(pcoa_info) <- c("axis1","axis2","eig")
sum_eig <- sum(pcoa_info$eig)
pcoa_info$eig <- pcoa_info$eig/sum_eig
pcoa_info <- pcoa_info[order(pcoa_info$eig,decreasing = T),]
#export result
write.csv(pcoa_info,file = "pcoa_info.csv")

#fig
library(dplyr)
library(ggplot2)
library(ggalt)
library(ggh4x)
library(ggpubr)


pcoa_info <- read.csv("pcoa_info.csv",header = T)
sample_data <- read.csv("sample_data.csv",header = T,row.names = 1)
sample_data$X <- rownames(sample_data)
pcoa_info <- left_join(pcoa_info,sample_data,by="X")
pcoa_info$group <- factor(pcoa_info$group,levels = c("Mock","DCA","NE","DCA_NE"),
                          labels = c("Mock","DCA","NE","DCA+NE"))

my_color <- c("Mock"='#0000C7',"NE"='#ED2024',"DCA"='#660D82',"DCA+NE"='#F091E6')


p1 <- ggplot(pcoa_info,aes(x=axis1,axis2,color=group)) +
  geom_point(size=1.8) + 
  scale_color_manual(values = my_color)+
  labs(x='PCoA1 (37.4%)',y='PCoA2 (16.1%)',color="Group")+
  annotate("text",x=Inf, y=Inf, size =3, vjust=1.5, hjust =1.2,
           label = expression(paste(italic(P), ' < 0.001')))+
  annotate("text",x=Inf, y=Inf, size = 3, vjust=2.5, hjust =1.2,
           label= expression(paste(R^2, ' = 0.416')))+
  scale_x_continuous(breaks = c(-0.1,0,0.1),limits = c(-0.15,0.15))+
  scale_y_continuous(breaks = c(-0.1,0,0.1))+
  #stat_centroid(aes(xend = axis1, yend = axis2, colour = group),
               # geom = "segment", crop_other = F,
              #  alpha=0.2,size = 1,show.legend = F)+
  #stat_chull(geom = "polygon", aes(group = group, fill = group), alpha = 0.1,show.legend = F)+
  geom_encircle(aes(fill=group), alpha = 0.1, show.legend = F,expand=0.03,spread=0.8,s_shape=0.5)+ 
  #stat_ellipse(aes(fill=group),level=0.95,type="norm",geom = "polygon",
   #            alpha=0.1,linetype='dashed',show.legend = FALSE)+
  scale_fill_manual(values = my_color)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 7,color = "black"),
        axis.text.x = element_text(size = 7,vjust = 0.9,color = "black"),
        axis.text.y = element_text(size = 7,color = "black"),
        panel.border = element_rect(linewidth = 0.5,color = 'black'),
        axis.ticks.x = element_line(linewidth = 0.3), 
        axis.ticks.y = element_line(linewidth  = 0.3),
        legend.text = element_text(size = 7,colour = "black"),
        legend.title = element_text(size = 8,colour = "black"),
        legend.spacing.x = unit(0.05, "cm"),
        plot.title = element_text(hjust = 0.5),
        legend.margin = margin(l=-5))
#geom_text(aes(label = SampleID), vjust = -0.5, hjust = 0.5)

p1 
library(ggh4x)
p1 <- p1 + force_panelsizes(total_width =unit(1.1, "in") ,total_height = unit(1.2, "in"))

p1 
#export fig
ggsave("pcoa.pdf",plot= p1,width = 5,height = 5)

ggsave("222.png",dpi=600,p1)

#figs of DEGs
#DEGs numbers---------------
library(ggplot2)
deg_num <- read.csv("DEGs_number.csv",header = T)
deg_num$group <- factor(deg_num$group,
              levels = c("DCA_Mock", "NE_Mock" ,"DCA_NE_Mock"),
              labels = c("DCA vs Mock", "NE vs Mock" ,"DCA+NE vs Mock"))

deg_num$Regulation <- factor(deg_num$Regulation,levels = c("Up","Down"))

mycolor <- c("Up"="#ED2024","Down"="#0000C7")

p <- ggplot(deg_num,aes(x=group,y=number))+
  geom_col(aes(fill=Regulation))+
  scale_fill_manual(values = mycolor)+
  labs(x="",y="")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.line.x = element_line(color = "black",linewidth = 0.5), 
        axis.line.y = element_line(color = "black",linewidth = 0.5),  
        axis.text.x = element_text(size = 8,angle = 45,hjust = 1,color = "black"),
        axis.text.y = element_text(size = 8,color = "black"),
        panel.border = element_rect(linewidth = 0.5,color = 'black'),  
        axis.ticks.x = element_line(linewidth = 0.3),  
        axis.ticks.y = element_line(linewidth = 0.3), 
  )
p
library(ggh4x)
p <- p + force_panelsizes(cols=unit(1.6, "in"),rows=unit(2, "in"))
p 
#export fig
ggsave("DEGs_number.pdf",plot= p,width = 5,height = 5)



#Venn plot----------------------
##Upset venn-----------
res_DCA_Mock <- read.csv("res_DCA_Mock.csv",header = T,row.names = 1)
res_NE_Mock <- read.csv("res_NE_Mock.csv",header = T,row.names = 1)
res_DCA_NE_Mock <- read.csv("res_DCA_NE_Mock.csv",header = T,row.names = 1)
res_ <- read.csv("res_.csv",header = T,row.names = 1)


#data modify
#res_DCA_Mock
res_DCA_Mock$fc <- abs(res_DCA_Mock$log2FoldChange)
res_DCA_Mock$gene <- rownames(res_DCA_Mock)
for (i in 1:nrow(res_DCA_Mock)){
  if (res_DCA_Mock$log2FoldChange[i] > 0){
    res_DCA_Mock$group_DCA_Mock[i] = "DCA"
  }
  else {
    res_DCA_Mock$group_DCA_Mock[i] = "Mock"
  }
}

res_DCA_Mock <- res_DCA_Mock[order(res_DCA_Mock$group,res_DCA_Mock$fc),]
res_DCA_Mock$"DCA vs Mock" <- 1
#res_DCA_Mock$gene <- factor(res_DCA_Mock$gene,levels = res_DCA_Mock$gene)
#res_DCA_Mock$group <- factor(res_DCA_Mock$group,levels = c("Mock","DCA"))


#res_NE_Mock
res_NE_Mock$fc <- abs(res_NE_Mock$log2FoldChange)
res_NE_Mock$gene <- rownames(res_NE_Mock)
for (i in 1:nrow(res_NE_Mock)){
  if (res_NE_Mock$log2FoldChange[i] > 0){
    res_NE_Mock$group_NE_Mock[i] = "NE"
  }
  else {
    res_NE_Mock$group_NE_Mock[i] = "Mock"
  }
}

res_NE_Mock <- res_NE_Mock[order(res_NE_Mock$group,res_NE_Mock$fc),]
res_NE_Mock$"NE vs Mock"<- 1
#res_NE_Mock$gene <- factor(res_NE_Mock$gene,levels = res_NE_Mock$gene)
#res_NE_Mock$group <- factor(res_NE_Mock$group,levels = c("Mock","NE"))


#res_DCA_NE_Mock
res_DCA_NE_Mock$fc <- abs(res_DCA_NE_Mock$log2FoldChange)
res_DCA_NE_Mock$gene <- rownames(res_DCA_NE_Mock)
for (i in 1:nrow(res_DCA_NE_Mock)){
  if (res_DCA_NE_Mock$log2FoldChange[i] > 0){
    res_DCA_NE_Mock$group_DCA_NE_Mock[i] = "DCA+NE"
  }
  else {
    res_DCA_NE_Mock$group_DCA_NE_Mock[i] = "DCA"
  }
}

res_DCA_NE_Mock <- res_DCA_NE_Mock[order(res_DCA_NE_Mock$group,res_DCA_NE_Mock$fc),]
res_DCA_NE_Mock$"DCA+NE vs DCA" <- 1
#res_DCA_NE_Mock$gene <- factor(res_DCA_NE_Mock$gene,levels = res_DCA_NE_Mock$gene)
#res_DCA_NE_Mock$group <- factor(res_DCA_NE_Mock$group,levels = c("DCA+NE","DCA"))


#res_
res_$fc <- abs(res_$log2FoldChange)
res_$gene <- rownames(res_)
for (i in 1:nrow(res_)){
  if (res_$log2FoldChange[i] > 0){
    res_$group_[i] = "DCA+NE"
  }
  else {
    res_$group_[i] = "NE"
  }
}

res_ <- res_[order(res_$group,res_$fc),]
res_$"DCA+NE vs NE"<- 1
#res_$gene <- factor(res_$gene,levels = res_$gene)
#res_$group <- factor(res_$group,levels = c("DCA+NE","NE"))


#modify data for Upset
library(dplyr)
data <- data.frame(gene=unique(c(res_DCA_Mock$gene,res_DCA_NE_Mock$gene,res_NE_Mock$gene,res_$gene)))
data <- left_join(data,res_DCA_Mock[,-c(1:3)],by="gene")
data <- left_join(data,res_NE_Mock[,-c(1:3)],by="gene")
data <- left_join(data,res_DCA_NE_Mock[,-c(1:3)],by="gene")
#data <- left_join(data,res_[,-c(1:3)],by="gene")
data$"DCA+NE vs NE" <- 0
compar <- c("DCA vs Mock","DCA+NE vs DCA","NE vs Mock","DCA+NE vs NE")
data[compar][is.na(data[compar])] <- 0

#fig for Upset
library(ggplot2)
library(ComplexUpset)

p1 <- upset(data, compar,
      width_ratio=0.35,
      set_sizes=(
        upset_set_size(
          geom=geom_bar(),
          position='right'
        )+
         geom_text(aes(label=after_stat(count)), hjust=1.1, stat='count')
        # you can also add annotations on top of bars:
        + annotate(geom='text', label='@', x='Drama', y=360, color='white', size=3)
        + expand_limits(y=360)
      ),
      matrix=(
        intersection_matrix(
          geom=geom_point(
            size=4
          ),
          segment=geom_segment(
            color='#660D82'
          ),
          outline_color=list(
            active='#660D82',
            inactive='grey70'
          )
        )+
          scale_color_manual(
            values=c('TRUE'='#660D82', 'FALSE'='grey'),
            labels=c('TRUE'='yes', 'FALSE'='no'),
            breaks=c('TRUE', 'FALSE'),
            name='Is intersection member?'
        )
      ),
      base_annotations=list(
        'Intersection size'=intersection_size(
          counts=TRUE,
          mapping=aes(fill='bars_color')
        ) + scale_fill_manual(values=c('bars_color'='#660D82'), guide='none')
      ),
    )

p1

ggsave("upset.pdf",plot = p1,width = 7,height = 4)


##ggVennDiagram Venn plot------------------------
library(ggVennDiagram)
library(ggplot2)
res_DCA_Mock <- read.csv("res_DCA_Mock001.csv",header = T)
res_NE_Mock <- read.csv("res_NE_Mock001.csv",header = T)
res_DCA_NE_Mock <- read.csv("res_DCA_NE_Mock001.csv",header = T)


DCA_Mock <- as.character(res_DCA_Mock$X)
NE_Mock <- as.character(res_NE_Mock$X)
DCA_NE_Mock <- as.character(res_DCA_NE_Mock$X)

gene_list <- list("DCA vs Mock"=DCA_Mock, "NE vs Mock"=NE_Mock,
                  "DCA+NE vs Mock"=DCA_NE_Mock,"DCA+NE vs NE"= as.character())


p <- ggVennDiagram(gene_list,label = "count")+
  scale_fill_gradient(low="white",high = "#ED2024")+
  scale_x_continuous(expand = expansion(mult = .08))+
  theme(legend.position = "none")
p
library(ggh4x)
p <- p + force_panelsizes(total_width=unit(4, "in"),total_height=unit(3, "in"))
p
#export fig
ggsave("venn_DEGs001.pdf",plot = p,width = 5,height = 5)

# Manhattan plot for DEGs-------------------------------
library(ggplot2)
library(dplyr)
library(ggrepel)

DCA_Mock <- read.csv("DCA_Mock.csv",header = T)
NE_Mock <- read.csv("NE_Mock.csv",header = T)
DCA_NE_Mock <- read.csv("DCA_NE_Mock.csv",header = T)


DCA_Mock$group <- "DCA_Mock"
NE_Mock$group <- "NE_Mock"
DCA_NE_Mock$group <- "DCA_NE_Mock"


colnames(DCA_Mock)[1] <- "gene"
colnames(NE_Mock)[1] <- "gene"
colnames(DCA_NE_Mock)[1] <- "gene"
colnames()[1] <- "gene"


#select gene will be shown on the plot (use yours criterion)
#top5sig1 <- filter(DCA_Mock,group=="DCA_Mock") %>% slice_max(abs(log2FoldChange),n=5)
#top5sig2 <- filter(NE_Mock,group=="NE_Mock") %>% slice_max(abs(log2FoldChange),n=5)
#top5sig3 <- filter(DCA_NE_Mock,group=="DCA_NE_Mock") %>% slice_max(abs(log2FoldChange),n=5)
#top5sig4 <- filter(,group=="") %>% slice_max(abs(log2FoldChange),n=5)


#add info about label gene
DCA_Mock$size <- case_when(!(DCA_Mock$gene %in% top5sig1$gene)~ 1,
                           DCA_Mock$gene %in% top5sig1$gene ~ 2)
NE_Mock$size <- case_when(!(NE_Mock$gene %in% top5sig1$gene)~ 1,
                           NE_Mock$gene %in% top5sig1$gene ~ 2)
DCA_NE_Mock$size <- case_when(!(DCA_NE_Mock$gene %in% top5sig1$gene)~ 1,
                             DCA_NE_Mock$gene %in% top5sig1$gene ~ 2)

#combine data
mdata <- rbind(DCA_Mock,NE_Mock,DCA_NE_Mock,)
mdata[is.na(mdata)] <- 1
mdata$group <- factor(mdata$group,levels = c("DCA_Mock","NE_Mock","DCA_NE_Mock",""),
                      labels = c("DCA vs Mock", "NE vs Mock" ,"DCA+NE vs DCA","DCA+NE vs NE"))
#mdata <- mdata[abs(mdata$log2FoldChange)>=1,]


#add p value label
#mdata$ptype<- ifelse(mdata$padj<0.01,"P-adj<0.01","P-adj>=0.01")
for (i in 1:nrow(mdata)) {
  if (mdata$padj[i] > 0.01 | abs(mdata$log2FoldChange[i])<2){
    mdata$ptype[i] = "p.adj > 0.01"
  }
  #else if (mdata$padj[i] >= 0.05 & mdata$log2FoldChange[i]<2){
   # mdata$ptype[i] = "p.adj<=0.01"
  #}
  else {
    mdata$ptype[i] = "p.adj <= 0.01"
  }
}

mdata$ptype <- factor(mdata$ptype,levels = c("p.adj > 0.01","p.adj <= 0.01"))
mycol1 <- c("gray","#ED2024")


#select genes with or without label
mdata_nolabel <- filter(mdata,size==1)
mdata_label <- filter(mdata,size==2)


#make data for background bar
range(DCA_Mock$log2FoldChange)
range(NE_Mock$log2FoldChange)
range(DCA_NE_Mock$log2FoldChange)

dfbar<-data.frame(x=c("DCA_Mock","NE_Mock","DCA_NE_Mock","" ),
                  y=c(3.822238,7.630585,5.502281,3.511136))
dfbar1<-data.frame(x=c("DCA_Mock","NE_Mock","DCA_NE_Mock","" ),
                   y=c(-3.388952,-8.133910,-5.373394,-4.315670))

#make data for color block
dfcol<-data.frame(x=c("DCA_Mock","NE_Mock","DCA_NE_Mock","" ),
                  y=0,
                  label=c("DCA\nvs\nMock","NE\nvs\nMock","DCA+NE\nvs\nDCA","DCA+NE\nvs\nNE" ))
mycol2 <- c('#0000C7','#66ff00','#660D82',"#FF000C")


#fig:
p <- ggplot()+
  geom_jitter(data = mdata,
                          aes(x = group, y = log2FoldChange, color = ptype),
                          size = 0.85,
                          width =0.4)+
  scale_color_manual(values = mycol1)+
  scale_y_continuous(breaks = c(-6,-2,2,6),limits = c(-8.2,8.2))+
  geom_hline(yintercept = c(-2,2),linetype = "dashed",linewidth=0.5,color="black")+
  geom_hline(yintercept = c(-6,6),linetype = "dashed",linewidth=0.5,color="#dcdcdc")+
  labs(x="",y="Log2FoldChange")+
  theme_bw()+
  theme(
    axis.title = element_text(size = 10,color = "black"),
    axis.line.y = element_line(color = "black"),
    axis.line.x = element_line(color = "black"),
    axis.text.x = element_text(color = "black"),
    panel.grid = element_blank(),
    legend.text = element_text(size = 9),
  )

p
ggsave("Manhattan.pdf",plot = p,width = 7,height = 5)


p <- ggplot()+
  #background bar
  geom_col(data = dfbar,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)+
  scale_y_continuous(breaks = c(-8,-6,-4,-2,0,2,4,6,8))+
  geom_hline(yintercept = c(-6,-2,2,6),linetype = "dashed",linewidth=0.5,color="#dcdcdc")+
  geom_col(data = dfbar1,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)+
  #point plot
  geom_jitter(data = mdata_nolabel,
              aes(x = group, y = log2FoldChange, color = ptype),
              size = 0.85,
              width =0.4)+
  scale_color_manual(values = mycol1)+
  geom_jitter(data = mdata_label,
              aes(x = group, y = log2FoldChange),
              color = "#ffb61e",
              size = 3,
              width =0.4)+
  #add color block
  geom_tile(data = dfcol,
            aes(x=x,y=y),
            height=1.8,
            color = "black",
            fill = mycol2,
            alpha = 0.8,
            show.legend = F)+
  geom_text(data=dfcol,
            aes(x=x,y=y,label=label),
            size = 3,
            color ="white"
  )+
  #add label
  geom_text_repel(
    data=mdata_label,
    aes(x = group, y = log2FoldChange,label=gene),
  ) +
  #legend customize
  labs(x="",y="Log2 Fold Change")+
  theme_minimal()+
  theme(
    axis.title = element_text(size = 10,
                              color = "black"),
    axis.line.y = element_line(color = "black",
                               linewidth = 1),
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "top",
    legend.direction = "vertical",
    legend.justification = c(1,0),
    legend.text = element_text(size = 9),
  )
p

library(plotly)
ggplotly(p)

ggsave("Manhattan.pdf",plot = p,width = 5,height = 6)

#volcano plot for DEGs---------------
library(dplyr)
library(ggplot2)

#import data
DCA_Mock <- read.csv("DCAvsMock.csv",header = T)
NE_Mock <- read.csv("NEvsMock.csv",header = T)
DCA_NE_Mock <- read.csv("DCA+NEvsMock.csv",header = T)
DCA_NE_NE <- read.csv("DCA+NEvsNE.csv",header = T)

DCA_Mock[is.na(DCA_Mock)] <- 1
NE_Mock[is.na(NE_Mock)] <- 1
DCA_NE_Mock[is.na(DCA_NE_Mock)] <- 1
DCA_NE_NE[is.na(DCA_NE_NE)] <- 1

DCA_Mock$'log10.p.adj' <- -log10(DCA_Mock$padj)
NE_Mock$'log10.p.adj' <- -log10(NE_Mock$padj)
DCA_NE_Mock$'log10.p.adj' <- -log10(DCA_NE_Mock$padj)
DCA_NE_NE$'log10.p.adj' <- -log10(DCA_NE_NE$padj)

DCA_Mock <- mutate(DCA_Mock,change = ifelse(padj <= 0.01&abs(log2FoldChange)>=1,
                          ifelse(log2FoldChange>0,"Up","Down"),"No change"))
NE_Mock <- mutate(NE_Mock,change = ifelse(padj <= 0.01&abs(log2FoldChange)>=1,
                           ifelse(log2FoldChange>0,"Up","Down"),"No change"))
DCA_NE_Mock <- mutate(DCA_NE_Mock,change = ifelse(padj <= 0.01&abs(log2FoldChange)>=1,
                           ifelse(log2FoldChange>0,"Up","Down"),"No change"))
DCA_NE_NE <- mutate(DCA_NE_NE,change = ifelse(padj <= 0.01&abs(log2FoldChange)>=1,
                                                  ifelse(log2FoldChange>0,"Up","Down"),"No change"))


DCA_Mock$compar <- "DCA vs Mock"
NE_Mock$compar <- "NE vs Mock"
DCA_NE_Mock$compar <- "DCA+NE vs Mock"
DCA_NE_NE$compar <- "DCA+NE vs NE"

data <- data.frame(rbind(DCA_Mock, NE_Mock, DCA_NE_Mock, DCA_NE_NE))

data$compar <- factor(data$compar,levels = c("DCA vs Mock", "NE vs Mock" ,"DCA+NE vs Mock", "DCA+NE vs NE"))

#fig
my_color <- c("Down"='#0000C7',"Up"='#ED2024',"No change"='gray')


p1 <- ggplot(data,aes(x=log2FoldChange,y=log10.p.adj))+
  geom_point(size=1,aes(color=change)) + 
  scale_color_manual(values = my_color)+
  labs(color="Regulation",y="-log10 p.adj")+
  facet_grid(~compar)+
  theme_bw()+
  scale_x_continuous(breaks = c(-5,-1,0,1,5),expand = c(0.05,0))+
  scale_y_continuous(breaks = c(0,2,20,40),expand = c(0.01,0))+
  geom_vline(xintercept=c(-1,1),color = "black", linetype = "dashed")+
  geom_hline(yintercept=c(-log10(0.01)),color = "black", linetype = "dashed")+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 8),
        legend.title = element_text(size = 8),
        axis.text.x = element_text(size = 7,vjust = 0.9,color = "black"),
        axis.text.y = element_text(size = 7,color = "black"),
        panel.border = element_rect(linewidth = 0.5,color = 'black'),
        axis.ticks.x = element_line(linewidth = 0.3), 
        axis.ticks.y = element_line(linewidth  = 0.3),
        legend.text = element_text(size = 7,colour = "black"),
        legend.spacing.x = unit(0.05, "cm"),
        plot.title = element_text(hjust = 0.5),
        strip.text = element_text(size = 8,color = "black"),
        legend.margin = margin(l=-5))
#geom_text(aes(label = SampleID), vjust = -0.5, hjust = 0.5)

p1 
library(ggh4x)
p1 <- p1 + force_panelsizes(total_width =unit(4.3, "in") ,total_height = unit(1.2, "in"))
p1 
#export fig
ggsave("volcano.pdf",plot= p1,width = 5,height = 5)

ggsave("111.png",dpi = 600,p1)


#heatmap for DEGs------------------------
library(ggplot2)
library(ggtext)
library(glue)
library(dplyr)
library(ggtree)
library(tidyr)
library(patchwork)

# prepare DEGs table-------------------------
res_DCA_Mock <- read.csv("res_DCA_Mock.csv",header = T)
res_NE_Mock <- read.csv("res_NE_Mock.csv",header = T)
res_DCA_NE_Mock <- read.csv("res_DCA_NE_Mock.csv",header = T)

gene_count <- read.csv("normalized_counts.csv",header = T)

# extract gene ID
gene_DCA_Mock <- data.frame("X"=res_DCA_Mock$X)
gene_NE_Mock <- data.frame("X"=res_NE_Mock$X)
gene_DCA_NE_Mock <- data.frame("X"=res_DCA_NE_Mock$X)

gene_id <- data.frame(unique(rbind(gene_DCA_Mock,gene_NE_Mock,
                               gene_DCA_NE_Mock)))


gene_count_sub <- left_join(gene_id,gene_count,by="X")
rownames(gene_count_sub) <- gene_count_sub$X

gene_count_sub <- as.data.frame(t(gene_count_sub[,-1]))

gene_count_sub$group <- rep(c("DCA","Mock","NE","DCA+NE"),each=6)

gene_count_sub <- aggregate(gene_count_sub[,-ncol(gene_count_sub)],
                            by=list(gene_count_sub$group),FUN="median")

rownames(gene_count_sub) <- gene_count_sub$Group.1

gene_count_sub <- as.data.frame(t(gene_count_sub[,-1]))

gene_count_sub <- gene_count_sub+1
#export
write.csv(gene_count_sub,file = "DEGs_heatmap.csv")

#import DEGs heatmap data
data <- read.csv("DEGs_heatmap.csv",header = T,row.names = 1)

#calculate Log2FC
data_empty <- data.frame(matrix(0,nrow = nrow(data),ncol = ncol(data)))
rownames(data_empty) <- rownames(data)
colnames(data_empty) <- colnames(data)
data_empty <- log2(data/data[,3])


##fig heatmap----------------------
#tree plot
tree_y <- hclust(dist(data_empty)) 

p1 <- ggtree(tree_y,layout="rectangular",branch.length="none",size=0.2)+
  theme(plot.margin = margin(c(0,0,0,0))) +
geom_tiplab()+ 
theme_tree2()+
scale_x_continuous(limits = c(0,25))
p1
#get label order
tree_d <- fortify(tree_y)
tree_d <- subset(tree_d, isTip)
lab_order <- tree_d$label[order(tree_d$y,decreasing=F)]


#data modify
data_empty$gene <- row.names(data_empty)
data <- pivot_longer(data_empty,cols = colnames(data_empty[,-ncol(data_empty)]),
                     names_to = "group",values_to = "log2fc")
data$group <- factor(data$group,levels = c("Mock","DCA","NE","DCA.NE"),
                     labels = c("Mock","DCA","NE","DCA+NE"))
data$gene <- factor(data$gene,levels = lab_order)

#fig
range(data$log2fc)

p2 <- ggplot(data,aes(x=group,y=gene))+
  geom_tile(aes(fill=log2fc),color="black")+
  scale_x_discrete(expand=c(0,0),position = "top")+
  scale_y_discrete(expand=c(0,0),position = "left")+
  scale_fill_gradient2(limits=c(-7,7),
                       breaks = c(-7,0,7),low="#0000C7", high="#ED2024", mid="white",midpoint = 0)+
  guides(fill=guide_colorbar(title = "log2FC",title.hjust = 0.5,title.vjust = 1))+
  theme(axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_text(margin = margin(r = 0),colour='black',size = 8,angle = 90,hjust = 1,vjust =0.5 ),
        axis.text.y = element_text(margin = margin(r = 0),colour='black',size = 8),
        axis.title = element_blank(),
        #legend.key.size = unit(0.6, "cm"),  #key size
        legend.key.height = unit(0.2, "cm"),
        legend.key.width = unit(0.3, "cm"),
        legend.text = element_text(size = 8, color = "black"), #legend label
        legend.spacing.x = unit(0.2, "cm"),#distance between key and label)
        legend.position = "bottom",
        #legend.margin = margin(b=-5),
        legend.title = element_text(size = 8,color = "black"),
        plot.margin = margin(0,0,0,0)
  ) +
  coord_flip()
p2

ggsave("DEGs_heatmap.pdf",plot = p2,width = 150,height = 6,limitsize = FALSE)
       
ggsave("DEGs_heatmap.png",plot = p2,width = 100,height = 6,dpi = 300,limitsize = FALSE)

##circlize heatmap--------------
library(circlize)
library(ComplexHeatmap)
col_fun = colorRamp2(c(-7, 0, 7), c("#0000C7", "white", "#ED2024"))
#import DEGs heatmap data
data <- read.csv("DEGs_heatmap.csv",header = T,row.names = 1)
#calculate Log2FC
data_empty <- data.frame(matrix(0,nrow = nrow(data),ncol = ncol(data)))
rownames(data_empty) <- rownames(data)
colnames(data_empty) <- colnames(data)
data_empty <- log2(data/data[,3])
data_empty <- as.matrix(data_empty)
mat <- as.matrix(data_empty[,c(2,4,1,3)])
colnames(mat) <- c("DCA+NE","NE","DCA","Mock")

#fig
circos.clear()
circos.heatmap(mat, col = col_fun,track.height = 0.2,cluster = TRUE,
               rownames.side = "inside",cell.border ="black"
               )

#legend
lgd = Legend(title = "log2FC", col_fun = col_fun)
circle_size = unit(10, "inch")
draw(lgd, x = circle_size,just="right")

#Enrichment analysis----------------------------
library(clusterProfiler)
library(org.Gg.eg.db)  #genome annotation of chicken

#import data
res_DCA_Mock <- read.csv("res_DCA_Mock001.csv",header = T)
res_NE_Mock <- read.csv("res_NE_Mock001.csv",header = T)
res_DCA_NE_Mock <- read.csv("res_DCA_NE_Mock001.csv",header = T)


# extract gene ID
gene_DCA_Mock <- res_DCA_Mock$X
gene_NE_Mock <- res_NE_Mock$X
gene_DCA_NE_Mock <- res_DCA_NE_Mock$X



# GO Enrichment Analysis---------------------------
ego_DCA_Mock <- enrichGO(gene = gene_DCA_Mock,
                OrgDb         = org.Gg.eg.db,
                keyType       = "SYMBOL",
                ont           = "ALL",
                pAdjustMethod = "fdr",
                pvalueCutoff  = 1,
                qvalueCutoff  = 1,
                readable      = TRUE)

ego_NE_Mock <- enrichGO(gene = gene_NE_Mock,
                        OrgDb         = org.Gg.eg.db,
                        keyType       = "SYMBOL",
                        ont           = "ALL",
                        pAdjustMethod = "fdr",
                        pvalueCutoff  = 1,
                        qvalueCutoff  = 1,
                        readable      = TRUE)

ego_DCA_NE_Mock <- enrichGO(gene = gene_DCA_NE_Mock,
                           OrgDb         = org.Gg.eg.db,
                           keyType       = "SYMBOL",
                           ont           = "ALL",
                           pAdjustMethod = "fdr",
                           pvalueCutoff  = 1,
                           qvalueCutoff  = 1,
                           readable      = TRUE)



#calculate richfactor
ego_DCA_Mock <- mutate(ego_DCA_Mock, richfactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))
ego_NE_Mock <- mutate(ego_NE_Mock, richfactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))
ego_DCA_NE_Mock <- mutate(ego_DCA_NE_Mock, richfactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))


#extract and export go enrichment results
go_result_DCA_Mock <- ego_DCA_Mock@result 
go_result_NE_Mock <- ego_NE_Mock@result
go_result_DCA_NE_Mock <- ego_DCA_NE_Mock@result

go_result_DCA_Mock <- go_result_DCA_Mock[which(go_result_DCA_Mock$p.adjust<0.05),]
go_result_NE_Mock <- go_result_NE_Mock[which(go_result_NE_Mock$p.adjust<0.05),]
go_result_DCA_NE_Mock <- go_result_DCA_NE_Mock[which(go_result_DCA_NE_Mock$p.adjust<0.05),]

write.csv(go_result_DCA_Mock,file = "go_DCA_Mock.csv")
write.csv(go_result_NE_Mock,file = "go_NE_Mock.csv")
write.csv(go_result_DCA_NE_Mock,file = "go_DCA_NE_Mock.csv")


#KEGG enrichment analysis----------------------
#translate gene ID
keytypes(org.Gg.eg.db)  #check key type (gene ID type)
gene_DCA_Mock_tr <- bitr(gene_DCA_Mock, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Gg.eg.db")
gene_DCA_Mock_tr <- gene_DCA_Mock_tr[!duplicated(gene_DCA_Mock_tr$SYMBOL),]

gene_NE_Mock_tr <- bitr(gene_NE_Mock, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Gg.eg.db")
gene_NE_Mock_tr <- gene_NE_Mock_tr[!duplicated(gene_NE_Mock_tr$SYMBOL),]

gene_DCA_NE_Mock_tr <- bitr(gene_DCA_NE_Mock, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Gg.eg.db")
gene_DCA_NE_Mock_tr <- gene_DCA_NE_Mock_tr[!duplicated(gene_DCA_NE_Mock_tr$SYMBOL),]


#richment analysis
kk_DCA_Mock <- enrichKEGG(gene = gene_DCA_Mock_tr$ENTREZID,
                         organism = "gga",
                         keyType = "ncbi-geneid",
                         pAdjustMethod = "fdr",
                         pvalueCutoff  = 1,
                         qvalueCutoff  = 1)

kk_NE_Mock <- enrichKEGG(gene = gene_NE_Mock_tr$ENTREZID,
                         organism = "gga",
                         keyType = "ncbi-geneid",
                         pAdjustMethod = "fdr",
                         pvalueCutoff  = 1,
                         qvalueCutoff  = 1)

kk_DCA_NE_Mock <- enrichKEGG(gene = gene_DCA_NE_Mock_tr$ENTREZID,
                            organism = "gga",
                            keyType = "ncbi-geneid",
                            pAdjustMethod = "fdr",
                            pvalueCutoff  = 1,
                            qvalueCutoff  = 1)


#transform gene ID to human readable
kk_DCA_Mock <- setReadable(kk_DCA_Mock, OrgDb = org.Gg.eg.db, keyType="ENTREZID")
kk_NE_Mock <- setReadable(kk_NE_Mock, OrgDb = org.Gg.eg.db, keyType="ENTREZID")
kk_DCA_NE_Mock <- setReadable(kk_DCA_NE_Mock, OrgDb = org.Gg.eg.db, keyType="ENTREZID")


#extract and export KEGG enrichment results
kegg_result_DCA_Mock <- kk_DCA_Mock@result 
kegg_result_NE_Mock <- kk_NE_Mock@result
kegg_result_DCA_NE_Mock <- kk_DCA_NE_Mock@result

kegg_result_DCA_Mock <- kegg_result_DCA_Mock[which(kegg_result_DCA_Mock$p.adjust<0.05),]
kegg_result_NE_Mock <- kegg_result_NE_Mock[which(kegg_result_NE_Mock$p.adjust<0.05),]
kegg_result_DCA_NE_Mock <- kegg_result_DCA_NE_Mock[which(kegg_result_DCA_NE_Mock$p.adjust<0.05),]

#export
write.csv(kegg_result_DCA_Mock,file = "kegg_DCA_Mock001.csv")
write.csv(kegg_result_NE_Mock,file = "kegg_NE_Mock001.csv")
write.csv(kegg_result_DCA_NE_Mock,file = "kegg_DCA_NE_Mock001.csv")

#figs for GO -----------------------
go_DCA_Mock <-read.csv("go_DCA_Mock.csv",header = T,row.names = 1)
go_NE_Mock <-read.csv("go_NE_Mock.csv",header = T,row.names = 1)
go_DCA_NE_Mock <-read.csv("go_DCA_NE_Mock.csv",header = T,row.names = 1)
# add comparison info
go_DCA_NE_Mock$compar <- "DCA+NE vs Mock"
go_NE_Mock$compar <- "NE vs Mock"

go_DCA_Mock[1,] <- NA
go_DCA_Mock$compar <- "DCA vs Mock"

#combine comparisons
total <- rbind(go_DCA_Mock,go_NE_Mock,go_DCA_NE_Mock)
total$compar <- factor(total$compar,levels = c("DCA vs Mock","NE vs Mock","DCA+NE vs Mock"))

#make colorful name
library(ggtext)
library(glue)
library(dplyr)

for (i in 1:nrow(total)) {
  if (is.na(total$ONTOLOGY[i])) {
    total$color[i] = "#000000"  
  } else if (total$ONTOLOGY[i] == "BP") {
    total$color[i] = "#085f63"
  } else if (total$ONTOLOGY[i] == "MF") {
    total$color[i] = "#FFAB00"
  } else {
    total$color[i] = "#ff5959"
  }
}

total <- mutate(total,font_name=glue(("<span style='color:{color}'>{Description}")))

total <- total[order(total$ONTOLOGY,total$Description,decreasing = T),]

total$font_name <- factor(total$font_name,levels = unique(total$font_name))

#modify p  
total$p_group <- cut(total$p.adjust, breaks = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                    labels = c("<0.0001","0.0001-0.001", "0.001-0.01", "0.01-0.05", ">0.05"))

mycolor <- c("<0.0001"="#FF2121","0.0001-0.001"="#FF5968", 
             "0.001-0.01"="#FF8BA7", "0.01-0.05"="#FFBCE0")

#fig
library(ggplot2)

p <- ggplot(data=total,aes(x=compar,y=font_name))+
  geom_point(aes(size=RichFactor,color=p_group))+
  scale_radius(limits = c(0.05,0.8),breaks = c(0.05,0.25,0.45,0.65,0.75), 
               range = c(2, 8))+
  scale_color_manual(values = mycolor)+
  labs(x='',y='',size="Rich Factor",color="p.adjust")+
  theme_bw()+
  theme(panel.grid.major.x = element_line(linetype = "dashed"),
        axis.text.x = element_markdown(size = 8,hjust = 1,angle = 45),
        axis.text.y = element_markdown(size = 8,color = "black"),
        panel.border = element_rect(linewidth = 0.5,color = 'black'),
        axis.ticks.x = element_line(linewidth = 0.3), 
        axis.ticks.y = element_line(linewidth  = 0.3),
        legend.position = "top"
  )+coord_flip()
p
library(ggh4x)

p <- p + force_panelsizes(rows=unit(1.5, "in"),cols=unit(11.77, "in"))

ggsave("go_plot.pdf",plot = p,width = 20,height = 5)
#0.11*row

#figs for KEGG -----------------------
kegg_DCA_Mock <-read.csv("kegg_DCA_Mock001.csv",header = T,row.names = 1)
kegg_NE_Mock <-read.csv("kegg_NE_Mock001.csv",header = T,row.names = 1)
kegg_DCA_NE_Mock <-read.csv("kegg_DCA_NE_Mock001.csv",header = T,row.names = 1)

# add comparison info
kegg_DCA_Mock$compar <- "DCA vs Mock"
kegg_NE_Mock$compar <- "NE vs Mock"
kegg_DCA_NE_Mock$compar <- "DCA+NE vs Mock"

#combine 4 comparisons
total2 <- rbind(kegg_DCA_Mock,kegg_NE_Mock,kegg_DCA_NE_Mock)
total2$compar <- factor(total2$compar,levels = c("DCA vs Mock","NE vs Mock","DCA+NE vs Mock","DCA+NE vs NE"))

#make colorful name
library(ggtext)
library(glue)
library(dplyr)
library(paletteer)
#pick colors
d_palettes <- palettes_d_names
#my_color <- paletteer_d("ggsci::default_ucscgb",n=25)
my_color <- c("#B71C1C","#EF9A9A","#6A1B9A","#E1BEE7","#303F9F","#6b6882",
              "#0277BD","#4FC3F7","#00838F","#80DEEA","#2E7D32","#A5D6A7",
              "#9E9D24","#E6EE9C","#E65100","#FFB74D","#4E342E","#BCAAA4",
              "#66ff00","#9900CCFF","#00FFFFFF","#800000FF","#7986CB",
              "#000080","#bacac6","#392f41","#7c4b00","#4c221b")


my_color1 <- c("Amino acid metabolism"="#B71C1C","Carbohydrate metabolism"="#EF9A9A",
               "Cell growth and death"="#6A1B9A","Endocrine system"="#E1BEE7",
               "Global and overview maps"="#303F9F","Lipid metabolism"="#6b6882",
               "Membrane transport"="#0277BD","Metabolism of cofactors and vitamins"="#4FC3F7",
               "Metabolism of other amino acids"="#00838F","Metabolism of terpenoids and polyketides"="#80DEEA",
               "Nucleotide metabolism"="#2E7D32","Replication and repair"="#A5D6A7",
               "Signaling molecules and interaction"="#9E9D24","Transport and catabolism"="#E6EE9C",
               "Xenobiotics biodegradation and metabolism"="#bacac6"
               )




#add font color for KEGG 1st level 
for (i in 1:nrow(total2)) {
  if (is.na(total2$category[i])) {
    total2$color[i] = "#000000"  
  } else if (total2$category[i] == "Metabolism") {
    total2$color[i] = my_color[1]
  } else if (total2$category[i] == "Organismal Systems") {
    total2$color[i] = my_color[3]
  } else if (total2$category[i] == "Environmental Information Processing") {
    total2$color[i] = my_color[5]
  } else if (total2$category[i] == "Human Diseases") {
    total2$color[i] = my_color[7]
  } else if (total2$category[i] == "Cellular Processes") {
    total2$color[i] = my_color[9]
  } else if (total2$category[i] == "Genetic Information Processing") {
    total2$color[i] = my_color[11]
  } else {
    total2$color[i] = my_color[13]
  }
}

total2 <- mutate(total2,font_name=glue(("<span style='color:{color}'>{Description}")))

total2 <- total2[order(total2$category,total2$subcategory,decreasing = T),]

total2$font_name <- factor(total2$font_name,levels = unique(total2$font_name))

#modify p  
total2$p_group <- cut(total2$p.adjust, breaks = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                     labels = c("<0.0001","0.0001-0.001", "0.001-0.01", "0.01-0.05", ">0.05"))

mycolor2 <- c("<0.0001"="#FF2121","0.0001-0.001"="#FF5968", 
              "0.001-0.01"="#FF8BA7", "0.01-0.05"="#FFBCE0")



#fig
library(ggplot2)

p2 <- ggplot(data=total2,aes(x=compar,y=font_name))+
  geom_point(aes(size=RichFactor,color=p_group))+
  scale_x_discrete(limits=c(levels(total2$compar)))+
  scale_radius(limits = c(0.15,0.45),breaks = c(0.15,0.3,0.45), 
               range = c(1, 5))+
  scale_color_manual(values = mycolor2)+
  labs(x='',y='',size="Rich Factor",color="p.adjust")+
  theme_bw()+
  theme(panel.grid.major.x = element_line(linetype = "dashed"),
        axis.text.x = element_markdown(size = 7,hjust = 1,angle = 45,color = "black"),
        #axis.text.y = element_markdown(size = 7,color = "black"),
        axis.text.y = element_blank(),
        #axis.text.x = element_blank(),
        panel.border = element_rect(linewidth = 0.5,color = 'black'),
        #axis.ticks.x = element_line(linewidth = 0.3), 
        axis.ticks.y = element_line(linewidth  = 0.3),
        axis.ticks.x = element_blank(),
        legend.title = element_text(size = 8,color = "black")
  )
p2
#draw annotation
total2$x_name <- "name"

p3 <- ggplot(total2,aes(x=x_name,y=font_name))+
  geom_tile(aes(fill=subcategory))+
  scale_fill_manual(values = my_color1)+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(position="left",expand = c(0,0)) +
  labs(x="",y="",fill="KEGG Level2")+
  theme_minimal()+
  guides(fill = guide_legend(ncol = 1))+
  theme(axis.text.x =element_blank(),
        #axis.text.y =element_blank(),
        axis.text.y =element_markdown(size = 7),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.key.size = unit(0.6, "cm"),  #key size
        legend.key.height = unit(0.2, "cm"),
        legend.key.width = unit(0.3, "cm"),
        legend.text = element_text(size = 7, color = "black"), #legend label
        legend.spacing.x = unit(0.2, "cm"),#distance between key and label)
        legend.position = "right",
        legend.title = element_text(size = 8,color = "black"),
        plot.margin = margin(0,0,0,0)
  )
p3

library(aplot)
library(ggh4x)
p4 <- p2 %>% 
  insert_left(p3,width = 0.03)
p4
ggsave("kegg_plot001.pdf",plot = p4,width =7,height = 4.5)
# 7*4.5


#GSEA analysis----------------------------
library(clusterProfiler)
library(org.Gg.eg.db)  #genome annotation of chicken

#import data
DCA_Mock <- read.csv("DCAvsMock.csv",header = T)
NE_Mock <- read.csv("NEvsMock.csv",header = T)
DCA_NE_Mock <- read.csv("DCA+NEvsMock.csv",header = T)

#make gene list for GSEA
gene_DCA_Mock <- DCA_Mock[,3]
names(gene_DCA_Mock) <- as.character(DCA_Mock[,1])
gene_DCA_Mock <- sort(gene_DCA_Mock,decreasing = T)

gene_NE_Mock <- NE_Mock[,3]
names(gene_NE_Mock) <- as.character(NE_Mock[,1])
gene_NE_Mock <- sort(gene_NE_Mock,decreasing = T)

gene_DCA_NE_Mock <- DCA_NE_Mock[,3]
names(gene_DCA_NE_Mock) <- as.character(DCA_NE_Mock[,1])
gene_DCA_NE_Mock <- sort(gene_DCA_NE_Mock,decreasing = T)


#GSEA analysis for GO----------------------------------
set.seed(031424)
ego_DCA_Mock <- gseGO(geneList     = gene_DCA_Mock,
              OrgDb        = org.Gg.eg.db,
              ont          = "ALL",
              keyType      = "SYMBOL",
              minGSSize    = 10,
              maxGSSize    = 1000,
              pvalueCutoff = 1,
              pAdjustMethod = "fdr",
              seed = T,
              verbose      = T)

set.seed(031424)
ego_NE_Mock <- gseGO(geneList     = gene_NE_Mock,
                      OrgDb        = org.Gg.eg.db,
                      ont          = "ALL",
                      keyType      = "SYMBOL",
                      minGSSize    = 10,
                      maxGSSize    = 1000,
                      pvalueCutoff = 1,
                      pAdjustMethod = "fdr",
                      verbose      = T)

set.seed(031424)
ego_DCA_NE_Mock <- gseGO(geneList     = gene_DCA_NE_Mock,
                      OrgDb        = org.Gg.eg.db,
                      ont          = "ALL",
                      keyType      = "SYMBOL",
                      minGSSize    = 10,
                      maxGSSize    = 1000,
                      pvalueCutoff = 1,
                      pAdjustMethod = "fdr",
                      verbose      = T)


#extract and export results
gsea_go_DCA_Mock <- ego_DCA_Mock@result[which(abs(ego_DCA_Mock@result$NES)>1
                                            &ego_DCA_Mock@result$p.adjust<0.05),]
gsea_go_NE_Mock <- ego_NE_Mock@result[which(abs(ego_NE_Mock@result$NES)>1
                                            &ego_NE_Mock@result$p.adjust<0.05),]
gsea_go_DCA_NE_Mock <- ego_DCA_NE_Mock@result[which(abs(ego_DCA_NE_Mock@result$NES)>1
                                            &ego_DCA_NE_Mock@result$p.adjust<0.05),]

write.csv(gsea_go_DCA_Mock,file = "gsea_go_DCA_Mock.csv")
write.csv(gsea_go_NE_Mock,file = "gsea_go_NE_Mock.csv")
write.csv(gsea_go_DCA_NE_Mock,file = "gsea_go_DCA_NE_Mock.csv")


#GSEA analysis for KEGG----------------------------------
#gene ID transformation
keytypes(org.Gg.eg.db)  #check key type (gene ID type)

#make gene list for GSEA
library(dplyr)

name_DCA_Mock <- DCA_Mock[,1] %>% bitr(fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Gg.eg.db")
colnames(DCA_Mock) <- c("SYMBOL",colnames(DCA_Mock)[-1])
name_DCA_Mock <- left_join(name_DCA_Mock,DCA_Mock[,c(1,3)],by="SYMBOL")
gene_DCA_Mock_tr <- name_DCA_Mock[,3]
names(gene_DCA_Mock_tr) <- as.character(name_DCA_Mock[,2])
gene_DCA_Mock_tr <- sort(gene_DCA_Mock_tr,decreasing = T)

name_NE_Mock <- NE_Mock[,1] %>% bitr(fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Gg.eg.db")
colnames(NE_Mock) <- c("SYMBOL",colnames(NE_Mock)[-1])
name_NE_Mock <- left_join(name_NE_Mock,NE_Mock[,c(1,3)],by="SYMBOL")
gene_NE_Mock_tr <- name_NE_Mock[,3]
names(gene_NE_Mock_tr) <- as.character(name_NE_Mock[,2])
gene_NE_Mock_tr <- sort(gene_NE_Mock_tr,decreasing = T)

name_DCA_NE_Mock <- DCA_NE_Mock[,1] %>% bitr(fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Gg.eg.db")
colnames(DCA_NE_Mock) <- c("SYMBOL",colnames(DCA_NE_Mock)[-1])
name_DCA_NE_Mock <- left_join(name_DCA_NE_Mock,DCA_NE_Mock[,c(1,3)],by="SYMBOL")
gene_DCA_NE_Mock_tr <- name_DCA_NE_Mock[,3]
names(gene_DCA_NE_Mock_tr) <- as.character(name_DCA_NE_Mock[,2])
gene_DCA_NE_Mock_tr <- sort(gene_DCA_NE_Mock_tr,decreasing = T)


#GSEA analysis
set.seed(031524)
kk2_DCA_Mock <- gseKEGG(geneList  = gene_DCA_Mock_tr,
                        keyType = "ncbi-geneid",
                        organism     = 'gga',
                        minGSSize    = 10,
                        maxGSSize = 1000,
                        pvalueCutoff = 1,
                        pAdjustMethod = "fdr",
                        seed = T,
                        verbose      = FALSE)

set.seed(031524)
kk2_NE_Mock <- gseKEGG(geneList  = gene_NE_Mock_tr,
                        keyType = "ncbi-geneid",
                        organism     = 'gga',
                        minGSSize    = 10,
                        maxGSSize = 1000,
                        pvalueCutoff = 1,
                        pAdjustMethod = "fdr",
                        seed = T,
                        verbose      = FALSE)

set.seed(031524)
kk2_DCA_NE_Mock <- gseKEGG(geneList  = gene_DCA_NE_Mock_tr,
                        keyType = "ncbi-geneid",
                        organism     = 'gga',
                        minGSSize    = 10,
                        maxGSSize = 1000,
                        pvalueCutoff = 1,
                        pAdjustMethod = "fdr",
                        seed = T,
                        verbose      = FALSE)

               
#transform gene ID to human readable
kk2_DCA_Mock <- setReadable(kk2_DCA_Mock, OrgDb = org.Gg.eg.db, keyType="ENTREZID")
kk2_NE_Mock <- setReadable(kk2_NE_Mock, OrgDb = org.Gg.eg.db, keyType="ENTREZID")
kk2_DCA_NE_Mock <- setReadable(kk2_DCA_NE_Mock, OrgDb = org.Gg.eg.db, keyType="ENTREZID")

#extract results
gsea_kegg_DCA_Mock <- kk2_DCA_Mock@result[which(abs(kk2_DCA_Mock@result$NES)>1
                                          &kk2_DCA_Mock@result$p.adjust<0.05),]
gsea_kegg_NE_Mock <- kk2_NE_Mock@result[which(abs(kk2_NE_Mock@result$NES)>1
                                          &kk2_NE_Mock@result$p.adjust<0.05),]
gsea_kegg_DCA_NE_Mock <- kk2_DCA_NE_Mock@result[which(abs(kk2_DCA_NE_Mock@result$NES)>1
                                          &kk2_DCA_NE_Mock@result$p.adjust<0.05),]

#add level1 and level2 info
library(KEGGREST)
library(tidyr)

for (i in 1:nrow(gsea_kegg_DCA_Mock)){
  info <- keggGet(gsea_kegg_DCA_Mock$ID[i])
  if (length(info) > 0 && !is.null(info[[1]]$CLASS)) {
    gsea_kegg_DCA_Mock$info[i] <- info[[1]]$CLASS
  } else {
    gsea_kegg_DCA_Mock$info[i] <- NA
  }
}
gsea_kegg_DCA_Mock <- gsea_kegg_DCA_Mock %>% separate(info,sep = "; ",into=c("category","subcategory"))


for (i in 1:nrow(gsea_kegg_NE_Mock)){
  info = keggGet(gsea_kegg_NE_Mock$ID[i])
  if (length(info) > 0 && !is.null(info[[1]]$CLASS)) {
    gsea_kegg_NE_Mock$info[i] <- info[[1]]$CLASS
  } else {
    gsea_kegg_NE_Mock$info[i] <- NA
  }
}
gsea_kegg_NE_Mock <- gsea_kegg_NE_Mock %>% separate(info,sep = "; ",into=c("category","subcategory"))


for (i in 1:nrow(gsea_kegg_DCA_NE_Mock)){
  info = keggGet(gsea_kegg_DCA_NE_Mock$ID[i])
  if (length(info) > 0 && !is.null(info[[1]]$CLASS)) {
    gsea_kegg_DCA_NE_Mock$info[i] <- info[[1]]$CLASS
  } else {
    gsea_kegg_DCA_NE_Mock$info[i] <- NA
  }
}
gsea_kegg_DCA_NE_Mock <- gsea_kegg_DCA_NE_Mock %>% separate(info,sep = "; ",into=c("category","subcategory"))



#export
write.csv(gsea_kegg_DCA_Mock,file = "gsea_kegg_DCA_Mock.csv")
write.csv(gsea_kegg_NE_Mock,file = "gsea_kegg_NE_Mock.csv")
write.csv(gsea_kegg_DCA_NE_Mock,file = "gsea_kegg_DCA_NE_Mock.csv")


#figs for GSEA GO -----------------------
go_DCA_Mock <-read.csv("gsea_go_DCA_Mock.csv",header = T,row.names = 1)
go_NE_Mock <-read.csv("gsea_go_NE_Mock.csv",header = T,row.names = 1)
go_DCA_NE_Mock <-read.csv("gsea_go_DCA_NE_Mock.csv",header = T,row.names = 1)

# add comparison info
go_DCA_Mock$compar <- "DCA vs Mock"
go_NE_Mock$compar <- "NE vs Mock"
go_DCA_NE_Mock$compar <- "DCA+NE vs Mock"

#combine 4 comparisons
total <- rbind(go_DCA_Mock,go_NE_Mock,go_DCA_NE_Mock)
total$compar <- factor(total$compar,levels = c("DCA vs Mock","NE vs Mock","DCA+NE vs Mock"))

#make colorful name
library(ggtext)
library(glue)
library(dplyr)

for (i in 1:nrow(total)) {
  if (is.na(total$ONTOLOGY[i])) {
    total$color[i] = "#000000"  
  } else if (total$ONTOLOGY[i] == "BP") {
    total$color[i] = "#085f63"
  } else if (total$ONTOLOGY[i] == "MF") {
    total$color[i] = "#b25d25"
  } else {
    total$color[i] = "#ff5959"
  }
}

total <- mutate(total,font_name=glue(("<span style='color:{color}'>{Description}")))

total <- total[order(total$ONTOLOGY,total$font_name,decreasing = T),]

total$font_name <- factor(total$font_name,levels = unique(total$font_name))


#modify p  
total$p_group <- cut(total$p.adjust, breaks = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                     labels = c("<0.0001","0.0001-0.001", "0.001-0.01", "0.01-0.05", ">0.05"))

total$p_group2 <- cut(total$p.adjust, breaks = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                     labels = c("<0.0001","0.0001-0.001", "0.001-0.01", "0.01-0.05", ">0.05"))

mycolor1 <- c("<0.0001"="#FF2121","0.0001-0.001"="#FF5968", "0.001-0.01"="#FF8BA7", 
              "0.01-0.05"="#FFBCE0", ">0.05"="#E8E8E8")
mycolor2 <- c("<0.0001"="#009912","0.0001-0.001"="#00BC12", "0.001-0.01"="#89FF93", 
              "0.01-0.05"="#D0FAC7", ">0.05"="#E8E8E8")
#fig
library(ggplot2)

p <- ggplot(data = total,aes(x=compar,y=font_name))+
  geom_point(shape=21,stroke=0.01,
             aes(x=compar,y=font_name,size=abs(NES),fill=p_group2))+
  geom_point(data=subset(total, NES > 0),shape=19,stroke=0.1,
             aes(x=compar,y=font_name,size=abs(NES),color=p_group))+
  scale_fill_manual(values = mycolor2)+
  scale_color_manual(values = mycolor1)+
  scale_radius(limits = c(1.3,2.6),breaks = c(1.3,1.6,2.0,2.3,2.6), 
               range = c(2, 8))+
  labs(x='',y='')+
  guides(color = guide_legend(title = "Positive NES p.adjust"),
         fill = guide_legend(title = "Negative NES p.adjust"),
         size = guide_legend(title = "NES"))+
  theme_bw()+
  theme(panel.grid.major.x = element_line(linetype = "dashed"),
        axis.text.x = element_markdown(size = 8,hjust = 1,angle = 45),
        axis.text.y = element_markdown(size = 8,color = "black"),
        panel.border = element_rect(linewidth = 0.5,color = 'black'),
        axis.ticks.x = element_line(linewidth = 0.3), 
        axis.ticks.y = element_line(linewidth  = 0.3),
        legend.position = "top"
  )+coord_flip()
p
library(ggh4x)

p <- p + force_panelsizes(rows=unit(1.5, "in"),cols=unit(29.37, "in"))

ggsave("gsea_go_plot.pdf",plot = p,width = 35,height = 5)
#0.11*row

#specific gene set plot
library(enrichplot)
p1 <- gseaplot2(ego_, geneSetID ="GO:1990904",title = "DCA+NE vs NE")
p1
ggsave("1990904.pdf",plot = p1,width = 5,height = 4)



#figs for GSEA KEGG -----------------------
kegg_DCA_Mock <-read.csv("gsea_kegg_DCA_Mock.csv",header = T,row.names = 1)
kegg_NE_Mock <-read.csv("gsea_kegg_NE_Mock.csv",header = T,row.names = 1)
kegg_DCA_NE_Mock <-read.csv("gsea_kegg_DCA_NE_Mock.csv",header = T,row.names = 1)


# add comparison info
kegg_DCA_Mock$compar <- "DCA vs Mock"
kegg_NE_Mock$compar <- "NE vs Mock"
kegg_DCA_NE_Mock$compar <- "DCA+NE vs Mock"

#combine 4 comparisons
total2 <- rbind(kegg_DCA_Mock,kegg_NE_Mock,kegg_DCA_NE_Mock)
total2$compar <- factor(total2$compar,levels = c("DCA vs Mock","NE vs Mock","DCA+NE vs Mock"))

#make colorful name
library(ggtext)
library(glue)
library(dplyr)
library(paletteer)
#pick colors
length(unique((total2$subcategory)))
d_palettes <- palettes_d_names
#my_color <- paletteer_d("ggsci::default_ucscgb",n=25)
my_color <- c("#B71C1C","#EF9A9A","#6A1B9A","#E1BEE7","#303F9F","#6b6882",
              "#0277BD","#4FC3F7","#00838F","#80DEEA","#2E7D32","#A5D6A7",
              "#9E9D24","#E6EE9C","#E65100","#FFB74D","#4E342E","#BCAAA4",
              "#66ff00","#9900CCFF","#00FFFFFF","#800000FF","#7986CB",
              "#000080")

my_color1 <- c("Amino acid metabolism"="#B71C1C","Carbohydrate metabolism"="#EF9A9A",
               "Cell growth and death"="#6A1B9A","Endocrine system"="#E1BEE7",
               "Global and overview maps"="#303F9F","Lipid metabolism"="#6b6882",
               "Membrane transport"="#0277BD","Metabolism of cofactors and vitamins"="#4FC3F7",
               "Metabolism of other amino acids"="#00838F","Metabolism of terpenoids and polyketides"="#80DEEA",
               "Nucleotide metabolism"="#2E7D32","Replication and repair"="#A5D6A7",
               "Signaling molecules and interaction"="#9E9D24","Transport and catabolism"="#E6EE9C",
               "Xenobiotics biodegradation and metabolism"="#bacac6",
               "Cell motility"="#E65100","Energy metabolism"="#FFB74D",
               "Folding, sorting and degradation"="#4E342E",
               "Glycan biosynthesis and metabolism"="#BCAAA4",
               "Immune system"="#66ff00","Infectious disease: viral"="#9900CCFF",
               "Transcription"="#00FFFFFF","Translation"="#800000FF"
)


#add font color for KEGG 1st level 
for (i in 1:nrow(total2)) {
  if (is.na(total2$category[i])) {
    total2$color[i] = "#000000"  
  } else if (total2$category[i] == "Metabolism") {
    total2$color[i] = my_color[1]
  } else if (total2$category[i] == "Organismal Systems") {
    total2$color[i] = my_color[3]
  } else if (total2$category[i] == "Environmental Information Processing") {
    total2$color[i] = my_color[5]
  } else if (total2$category[i] == "Human Diseases") {
    total2$color[i] = my_color[7]
  } else if (total2$category[i] == "Cellular Processes") {
    total2$color[i] = my_color[9]
  } else if (total2$category[i] == "Genetic Information Processing") {
    total2$color[i] = my_color[11]
  } else {
    total2$color[i] = my_color[13]
  }
}


total2 <- mutate(total2,font_name=glue(("<span style='color:{color}'>{Description}")))

total2 <- total2[order(total2$category,total2$subcategory,decreasing = T),]

total2$font_name <- factor(total2$font_name,levels = unique(total2$font_name))

#modify p  
total2$p_group <- cut(total2$p.adjust, breaks = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                     labels = c("<0.0001","0.0001-0.001", "0.001-0.01", "0.01-0.05", ">0.05"))

total2$p_group2 <- cut(total2$p.adjust, breaks = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                      labels = c("<0.0001","0.0001-0.001", "0.001-0.01", "0.01-0.05", ">0.05"))

mycolor1 <- c("<0.0001"="#FF2121","0.0001-0.001"="#FF5968", "0.001-0.01"="#FF8BA7", 
              "0.01-0.05"="#FFBCE0", ">0.05"="#E8E8E8")
mycolor2 <- c("<0.0001"="#008C10","0.0001-0.001"="#45CC2D", "0.001-0.01"="#6CF093", 
              "0.01-0.05"="#C1F9C5", ">0.05"="#E8E8E8")

#fig
library(ggplot2)

data1 <- subset(total2, NES > 0)
data2 <- subset(total2, NES < 0)


p2 <- ggplot(total2,aes(x=compar,y=font_name))+
  geom_point(shape=21,stroke=0.01,
             aes(x=compar,y=font_name,size=abs(NES),fill=p_group2))+
  geom_point(data=subset(total2, NES > 0),shape=19,stroke=0.1,
             aes(x=compar,y=font_name,size=abs(NES),color=p_group))+
  scale_fill_manual(values = mycolor2)+
  scale_color_manual(values = mycolor1)+
  scale_radius(limits = c(1.4,2.6),breaks = c(1.4,1.8,2.2,2.6), 
               range = c(2, 8))+
  labs(x='',y='')+
  guides(color = guide_legend(title = "Positive NES p.adjust"),
         fill = guide_legend(title = "Negative NES p.adjust"),
         size = guide_legend(title = "NES"))+
  theme_bw()+
  theme(panel.grid.major.x = element_line(linetype = "dashed"),
        #axis.text.x = element_markdown(size = 8,hjust = 1,angle = 45),
        axis.text.y = element_markdown(size = 8,color = "black"),
        axis.text.x = element_blank(),
        panel.border = element_rect(linewidth = 0.5,color = 'black'),
        #axis.ticks.x = element_line(linewidth = 0.3), 
        axis.ticks.y = element_line(linewidth  = 0.3),
        axis.ticks.x = element_blank(),
        legend.position = "top"
  )+coord_flip()
p2
#draw annotation
total2$x_name <- "name"

p3 <- ggplot(total2,aes(x=x_name,y=font_name))+
  geom_tile(aes(fill=subcategory))+
  scale_fill_manual(values = my_color1)+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(position="left",expand = c(0,0)) +
  labs(x="",y="",fill="KEGG Level2")+
  theme_minimal()+
  guides(fill = guide_legend(ncol = 1))+
  theme(axis.text.x =element_markdown(size = 8,hjust=1,angle = 45),
        axis.text.y =element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.key.size = unit(0.6, "cm"),  #key size
        legend.key.height = unit(0.2, "cm"),
        legend.key.width = unit(0.3, "cm"),
        legend.text = element_text(size = 8, color = "black"), #legend label
        legend.spacing.x = unit(0.2, "cm"),#distance between key and label)
        legend.position = "right",
        legend.title = element_text(size = 9,color = "black"),
        plot.margin = margin(0,0,0,0)
  )+coord_flip()
p3
library(ggh4x)
library(cowplot)
p2 <- p2 + force_panelsizes(rows=unit(1.5, "in"),cols=unit(12.21, "in"))
p3 <- p3 + force_panelsizes(rows=unit(0.1, "in"),cols=unit(12.21, "in"))
p4 <- ggdraw()+
  draw_plot(p2,x=0,y=0.2)+
  draw_plot(p3,x=0.119,y=-0.198)

ggsave("gsea_kegg_plot.pdf",plot = p4,width =15,height = 5)
# 0.11*row

#specific gene set plot
library(enrichplot)
p5 <- gseaplot2(kk2_DCA_Mock, geneSetID =c(1:2),title = "DCA vs Mock")
ggsave("geneset_kegg.pdf",plot = p5,width = 5,height = 4)







