############ HCC Methylation ############
library(ggplot2)
library(tidyr)
library(reshape2) 
library(dplyr)
library(tibble)
library(ggpubr)
library(ggsci)

mycolour <- c('#af2337', '#ecc342', '#2967a0', '#96b437','#2f3c28', 
              '#da93ab','#e58932', '#80598f', '#7e331f', '#3b855a',
              '#c0b286', '#a9c9ed', '#ec977f', '#848482', '#604628',
              '#d26034', '#a64c6b', '#dbd245', '#eba83b', '#5d5092',
              '#222222', '#f2f3f4')
# basic statistics

## PCA

### 5mC

PCA <- read.table("path/WGBS_meth_PCA.txt", header = T, row.names=1)
PCA_2 <- t(PCA)
PCA_3 <- PCA_2[,which(apply(PCA_2,2,var)!=0)]
data.pca <- prcomp(PCA_3,scale = T)
pca.data = data.frame(data.pca$x)
pdf(file = "path/WGBS_PCA.pdf", width = 8, height = 3)
ggplot(pca.data, aes(x = PC1, y = PC2, color = c(rep("Neg_N",4),rep("Neg_T",4),
                                                   rep("Pos_N",4),rep("Pos_T",4)))) +
  geom_point(size = 1) +
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0) +
  stat_ellipse(aes(x = PC1, y = PC2), linetype = 2, size = 0.5, level = 0.95) + 
  theme_bw()+
  scale_color_manual(values = c("#af2337","#ecc342","#2967a0","#96b437"))+
  theme(axis.text.x = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(axis.text.y = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(axis.title.y = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(legend.text = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(legend.title = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5))
dev.off()

pca.data_2 <- pca.data[1:8,]
pdf(file = "path/WGBS_PCA_Neg.pdf", width = 6, height = 3)
ggplot(pca.data_2, aes(x = PC1, y = PC2, color = c(rep("Neg_N",4),rep("Neg_T",4)))) +
  geom_point(size = 1) +
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0) +
  stat_ellipse(aes(x = PC1, y = PC2), linetype = 2, size = 0.5, level = 0.95) + 
  theme_bw()+
  scale_color_manual(values = c("#af2337","#ecc342"))+
  theme(axis.text.x = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(axis.text.y = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(axis.title.y = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(legend.text = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(legend.title = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5))
dev.off()

### 5hmC

PCA_2 <- read.table("path/oxWGBS_meth_PCA.txt", header = T, row.names=1)
PCA_2_2 <- t(PCA_2)
PCA_2_3 <- PCA_2_2[,which(apply(PCA_2_2,2,var)!=0)]
data.pca_2 <- prcomp(PCA_2_3,scale = T)

pca.data_2 = data.frame(data.pca_2$x)
pdf(file = "path/oxWGBS_PCA.pdf", width = 8, height = 3)
ggplot(pca.data_2, aes(x = PC1, y = PC2, color = c(rep("Neg_N",4),rep("Neg_T",4),
                                                   rep("Pos_N",4),rep("Pos_T",4)))) +
  geom_point(size = 1) +
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0) +
  stat_ellipse(aes(x = PC1, y = PC2), linetype = 2, size = 0.5, level = 0.95) + 
  theme_bw()+
  scale_color_manual(values = c("#af2337","#ecc342","#2967a0","#96b437"))+
  theme(axis.text.x = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(axis.text.y = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(axis.title.y = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(legend.text = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(legend.title = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5))
dev.off()

pca.data_2_2 <- pca.data_2[1:8,]
pdf(file = "path/oxWGBS_PCA_Neg.pdf", width = 6, height = 3)
ggplot(pca.data_2_2, aes(x = PC1, y = PC2, color = c(rep("Neg_N",4),rep("Neg_T",4)))) +
  geom_point(size = 1) +
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0) +
  stat_ellipse(aes(x = PC1, y = PC2), linetype = 2, size = 0.5, level = 0.95) + 
  theme_bw()+
  scale_color_manual(values = c("#af2337","#ecc342"))+
  theme(axis.text.x = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(axis.text.y = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(axis.title.y = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(legend.text = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(legend.title = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5))
dev.off()

## Count whether methylation occurs

### 5mC
#### group level
WGBS_non <- read.table("path/WGBS_5mC_non_5mC.txt", header = T)
WGBS_non$sample <- factor(WGBS_non$sample,levels =c("Neg_N","Neg_T","Pos_N","Pos_T"))
pdf(file = "path/WGBS_5mC_non.pdf", width = 3, height = 3)
ggplot(WGBS_non, aes(x = sample, y = rate, fill = group))+
  geom_bar(stat="identity")+
  labs(x = "", y = "Proportion of modified cytosines (%)", title = "")+
  theme_bw()+
  theme()+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(values = c("#af2337","#2967a0"))+
  theme_bw()+theme(panel.grid =element_blank())+ theme(legend.text.align = 0)+
  theme(axis.text.x = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5, angle = 90))+
  theme(axis.text.y = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(axis.title.y = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(legend.text = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(legend.title = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5))
dev.off()

#### sample level
WGBS_non_sample <- read.table("path/WGBS_5mC_non_5mC_sample.txt", header = T)
WGBS_non_sample$sample <- factor(WGBS_non_sample$sample,levels =c("Neg_1_N", "Neg_2_N", "Neg_3_N", "Neg_4_N",
                                    "Neg_1_T", "Neg_2_T", "Neg_3_T", "Neg_4_T",
                                    "Pos_1_N", "Pos_2_N", "Pos_3_N", "Pos_4_N",
                                    "Pos_1_T", "Pos_2_T", "Pos_3_T", "Pos_4_T"))
pdf(file = "path/WGBS_5mC_non_sample.pdf", width = 4, height = 3)
ggplot(WGBS_non, aes(x = sample, y = rate, fill = group))+
  geom_bar(stat="identity")+
  labs(x = "", y = "Proportion of modified cytosines (%)", title = "")+
  theme_bw()+
  theme()+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(values = c("#af2337","#2967a0"))+
  theme_bw()+theme(panel.grid =element_blank())+ theme(legend.text.align = 0)+
  theme(axis.text.x = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5, angle = 90))+
  theme(axis.text.y = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(axis.title.y = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(legend.text = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(legend.title = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5))
dev.off()

### 5hmC

#### group level
oxWGBS_non <- read.table("path/oxWGBS_5mC_non_5mC.txt", header = T)
oxWGBS_non$sample <- factor(oxWGBS_non$sample,levels =c("Neg_N","Neg_T","Pos_N","Pos_T"))
pdf(file = "path/oxWGBS_5mC_non.pdf", width = 3, height = 3)
ggplot(oxWGBS_non, aes(x = sample, y = rate, fill = group))+
  geom_bar(stat="identity")+
  labs(x = "", y = "Proportion of modified cytosines (%)", title = "")+
  theme_bw()+
  theme()+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(values = c("#af2337","#2967a0"))+
  theme_bw()+theme(panel.grid =element_blank())+ theme(legend.text.align = 0)+
  theme(axis.text.x = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5, angle = 90))+
  theme(axis.text.y = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(axis.title.y = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(legend.text = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(legend.title = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5))
dev.off()

#### sample level
oxWGBS_non_sample <- read.table("path/oxWGBS_5mC_non_5mC_sample.txt", header = T)
oxWGBS_non_sample$sample <- factor(oxWGBS_non_sample$sample,levels =c("Neg_1_N", "Neg_2_N", "Neg_3_N", "Neg_4_N",
                                                    "Neg_1_T", "Neg_2_T", "Neg_3_T", "Neg_4_T",
                                                    "Pos_1_N", "Pos_2_N", "Pos_3_N", "Pos_4_N",
                                                    "Pos_1_T", "Pos_2_T", "Pos_3_T", "Pos_4_T"))
pdf(file = "path/oxWGBS_5mC_non_sample.pdf", width = 4, height = 3)
ggplot(oxWGBS_non_sample, aes(x = sample, y = rate, fill = group))+
  geom_bar(stat="identity")+
  labs(x = "", y = "Proportion of modified cytosines (%)", title = "")+
  theme_bw()+
  theme()+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(values = c("#af2337","#2967a0"))+
  theme_bw()+theme(panel.grid =element_blank())+ theme(legend.text.align = 0)+
  theme(axis.text.x = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5, angle = 90))+
  theme(axis.text.y = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(axis.title.y = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(legend.text = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(legend.title = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5))
dev.off()

## The distribution of methylation density on chromosomes

library("CMplot")

### 5mC

#### Neg_N

data <- read.table("path/meth_rate/WGBS_Neg_N.txt", header = T)
pdf(file = "path/WGBS_Neg_N_meth_site_distribution.pdf", height = 4, width = 6)
CMplot(
  data, plot.type="d",  bin.size=1e6, col=c("darkgreen", "yellow", "red"),
  file="jpg", dpi=300, file.output=F, verbose=TRUE
)
dev.off()
rm(data)

#### Neg_T

data <- read.table("path/meth_rate/WGBS_Neg_T.txt", header = T)
pdf(file = "path/WGBS_Neg_T_meth_site_distribution.pdf", height = 4, width = 6)
CMplot(
  data, plot.type="d",  bin.size=1e6, col=c("darkgreen", "yellow", "red"),
  file="jpg", dpi=300, file.output=F, verbose=TRUE
)
dev.off()
rm(data)

#### Neg_N_1

data <- read.table("path/meth_rate/WGBS_Neg_N_1.txt", header = T)
pdf(file = "path/WGBS_Neg_N_1_meth_site_distribution.pdf", height = 4, width = 6)
CMplot(
  data, plot.type="d",  bin.size=1e6, col=c("darkgreen", "yellow", "red"),
  file="jpg", dpi=300, file.output=F, verbose=TRUE
)
dev.off()
rm(data)

#### Neg_N_2

data <- read.table("path/meth_rate/WGBS_Neg_N_2.txt", header = T)
pdf(file = "path/WGBS_Neg_N_2_meth_site_distribution.pdf", height = 4, width = 6)
CMplot(
  data, plot.type="d",  bin.size=1e6, col=c("darkgreen", "yellow", "red"),
  file="jpg", dpi=300, file.output=F, verbose=TRUE
)
dev.off()
rm(data)

#### Neg_N_3

data <- read.table("path/meth_rate/WGBS_Neg_N_3.txt", header = T)
pdf(file = "path/WGBS_Neg_N_3_meth_site_distribution.pdf", height = 4, width = 6)
CMplot(
  data, plot.type="d",  bin.size=1e6, col=c("darkgreen", "yellow", "red"),
  file="jpg", dpi=300, file.output=F, verbose=TRUE
)
dev.off()
rm(data)

#### Neg_N_4

data <- read.table("path/meth_rate/WGBS_Neg_N_4.txt", header = T)
pdf(file = "path/WGBS_Neg_N_4_meth_site_distribution.pdf", height = 4, width = 6)
CMplot(
  data, plot.type="d",  bin.size=1e6, col=c("darkgreen", "yellow", "red"),
  file="jpg", dpi=300, file.output=F, verbose=TRUE
)
dev.off()
rm(data)

#### Neg_T_1

data <- read.table("path/meth_rate/WGBS_Neg_T_1.txt", header = T)
pdf(file = "path/WGBS_Neg_T_1_meth_site_distribution.pdf", height = 4, width = 6)
CMplot(
  data, plot.type="d",  bin.size=1e6, col=c("darkgreen", "yellow", "red"),
  file="jpg", dpi=300, file.output=F, verbose=TRUE
)
dev.off()
rm(data)

#### Neg_T_2

data <- read.table("path/meth_rate/WGBS_Neg_T_2.txt", header = T)
pdf(file = "path/WGBS_Neg_T_2_meth_site_distribution.pdf", height = 4, width = 6)
CMplot(
  data, plot.type="d",  bin.size=1e6, col=c("darkgreen", "yellow", "red"),
  file="jpg", dpi=300, file.output=F, verbose=TRUE
)
dev.off()
rm(data)

#### Neg_T_3

data <- read.table("path/meth_rate/WGBS_Neg_T_3.txt", header = T)
pdf(file = "path/WGBS_Neg_T_3_meth_site_distribution.pdf", height = 4, width = 6)
CMplot(
  data, plot.type="d",  bin.size=1e6, col=c("darkgreen", "yellow", "red"),
  file="jpg", dpi=300, file.output=F, verbose=TRUE
)
dev.off()
rm(data)

#### Neg_T_4

data <- read.table("path/meth_rate/WGBS_Neg_T_4.txt", header = T)
pdf(file = "path/WGBS_Neg_T_4_meth_site_distribution.pdf", height = 4, width = 6)
CMplot(
  data, plot.type="d",  bin.size=1e6, col=c("darkgreen", "yellow", "red"),
  file="jpg", dpi=300, file.output=F, verbose=TRUE
)
dev.off()
rm(data)

#### Pos_N

data <- read.table("path/meth_rate/WGBS_Pos_N.txt", header = T)
pdf(file = "path/WGBS_Pos_N_meth_site_distribution.pdf", height = 4, width = 6)
CMplot(
  data, plot.type="d",  bin.size=1e6, col=c("darkgreen", "yellow", "red"),
  file="jpg", dpi=300, file.output=F, verbose=TRUE
)
dev.off()
rm(data)

#### Pos_T

data <- read.table("path/meth_rate/WGBS_Pos_T.txt", header = T)
pdf(file = "path/WGBS_Pos_T_meth_site_distribution.pdf", height = 4, width = 6)
CMplot(
  data, plot.type="d",  bin.size=1e6, col=c("darkgreen", "yellow", "red"),
  file="jpg", dpi=300, file.output=F, verbose=TRUE
)
dev.off()
rm(data)

#### Pos_N_1

data <- read.table("path/meth_rate/WGBS_Pos_N_1.txt", header = T)
pdf(file = "path/WGBS_Pos_N_1_meth_site_distribution.pdf", height = 4, width = 6)
CMplot(
  data, plot.type="d",  bin.size=1e6, col=c("darkgreen", "yellow", "red"),
  file="jpg", dpi=300, file.output=F, verbose=TRUE
)
dev.off()
rm(data)

#### Pos_N_2

data <- read.table("path/meth_rate/WGBS_Pos_N_2.txt", header = T)
pdf(file = "path/WGBS_Pos_N_2_meth_site_distribution.pdf", height = 4, width = 6)
CMplot(
  data, plot.type="d",  bin.size=1e6, col=c("darkgreen", "yellow", "red"),
  file="jpg", dpi=300, file.output=F, verbose=TRUE
)
dev.off()
rm(data)

#### Pos_N_3

data <- read.table("path/meth_rate/WGBS_Pos_N_3.txt", header = T)
pdf(file = "path/WGBS_Pos_N_3_meth_site_distribution.pdf", height = 4, width = 6)
CMplot(
  data, plot.type="d",  bin.size=1e6, col=c("darkgreen", "yellow", "red"),
  file="jpg", dpi=300, file.output=F, verbose=TRUE
)
dev.off()
rm(data)

#### Pos_N_4

data <- read.table("path/meth_rate/WGBS_Pos_N_4.txt", header = T)
pdf(file = "path/WGBS_Pos_N_4_meth_site_distribution.pdf", height = 4, width = 6)
CMplot(
  data, plot.type="d",  bin.size=1e6, col=c("darkgreen", "yellow", "red"),
  file="jpg", dpi=300, file.output=F, verbose=TRUE
)
dev.off()
rm(data)

#### Pos_T_1

data <- read.table("path/meth_rate/WGBS_Pos_T_1.txt", header = T)
pdf(file = "path/WGBS_Pos_T_1_meth_site_distribution.pdf", height = 4, width = 6)
CMplot(
  data, plot.type="d",  bin.size=1e6, col=c("darkgreen", "yellow", "red"),
  file="jpg", dpi=300, file.output=F, verbose=TRUE
)
dev.off()
rm(data)

#### Pos_T_2

data <- read.table("path/meth_rate/WGBS_Pos_T_2.txt", header = T)
pdf(file = "path/WGBS_Pos_T_2_meth_site_distribution.pdf", height = 4, width = 6)
CMplot(
  data, plot.type="d",  bin.size=1e6, col=c("darkgreen", "yellow", "red"),
  file="jpg", dpi=300, file.output=F, verbose=TRUE
)
dev.off()
rm(data)

#### Pos_T_3

data <- read.table("path/meth_rate/WGBS_Pos_T_3.txt", header = T)
pdf(file = "path/WGBS_Pos_T_3_meth_site_distribution.pdf", height = 4, width = 6)
CMplot(
  data, plot.type="d",  bin.size=1e6, col=c("darkgreen", "yellow", "red"),
  file="jpg", dpi=300, file.output=F, verbose=TRUE
)
dev.off()
rm(data)

#### Pos_T_4

data <- read.table("path/meth_rate/WGBS_Pos_T_4.txt", header = T)
pdf(file = "path/WGBS_Pos_T_4_meth_site_distribution.pdf", height = 4, width = 6)
CMplot(
  data, plot.type="d",  bin.size=1e6, col=c("darkgreen", "yellow", "red"),
  file="jpg", dpi=300, file.output=F, verbose=TRUE
)
dev.off()
rm(data)

### 5hmC

#### Neg_N

data <- read.table("path/meth_rate/oxWGBS_Neg_N.txt", header = T)
pdf(file = "path/oxWGBS_Neg_N_meth_site_distribution.pdf", height = 4, width = 6)
CMplot(
  data, plot.type="d",  bin.size=1e6, col=c("darkgreen", "yellow", "red"),
  file="jpg", dpi=300, file.output=F, verbose=TRUE
)
dev.off()
rm(data)

#### Neg_T

data <- read.table("path/meth_rate/oxWGBS_Neg_T.txt", header = T)
pdf(file = "path/oxWGBS_Neg_T_meth_site_distribution.pdf", height = 4, width = 6)
CMplot(
  data, plot.type="d",  bin.size=1e6, col=c("darkgreen", "yellow", "red"),
  file="jpg", dpi=300, file.output=F, verbose=TRUE
)
dev.off()
rm(data)

#### Neg_N_1

data <- read.table("path/meth_rate/oxWGBS_Neg_N_1.txt", header = T)
pdf(file = "path/oxWGBS_Neg_N_1_meth_site_distribution.pdf", height = 4, width = 6)
CMplot(
  data, plot.type="d",  bin.size=1e6, col=c("darkgreen", "yellow", "red"),
  file="jpg", dpi=300, file.output=F, verbose=TRUE
)
dev.off()
rm(data)

#### Neg_N_2

data <- read.table("path/meth_rate/oxWGBS_Neg_N_2.txt", header = T)
pdf(file = "path/oxWGBS_Neg_N_2_meth_site_distribution.pdf", height = 4, width = 6)
CMplot(
  data, plot.type="d",  bin.size=1e6, col=c("darkgreen", "yellow", "red"),
  file="jpg", dpi=300, file.output=F, verbose=TRUE
)
dev.off()
rm(data)

#### Neg_N_3

data <- read.table("path/meth_rate/oxWGBS_Neg_N_3.txt", header = T)
pdf(file = "path/oxWGBS_Neg_N_3_meth_site_distribution.pdf", height = 4, width = 6)
CMplot(
  data, plot.type="d",  bin.size=1e6, col=c("darkgreen", "yellow", "red"),
  file="jpg", dpi=300, file.output=F, verbose=TRUE
)
dev.off()
rm(data)

#### Neg_N_4

data <- read.table("path/meth_rate/oxWGBS_Neg_N_4.txt", header = T)
pdf(file = "path/oxWGBS_Neg_N_4_meth_site_distribution.pdf", height = 4, width = 6)
CMplot(
  data, plot.type="d",  bin.size=1e6, col=c("darkgreen", "yellow", "red"),
  file="jpg", dpi=300, file.output=F, verbose=TRUE
)
dev.off()
rm(data)

#### Neg_T_1

data <- read.table("path/meth_rate/oxWGBS_Neg_T_1.txt", header = T)
pdf(file = "path/oxWGBS_Neg_T_1_meth_site_distribution.pdf", height = 4, width = 6)
CMplot(
  data, plot.type="d",  bin.size=1e6, col=c("darkgreen", "yellow", "red"),
  file="jpg", dpi=300, file.output=F, verbose=TRUE
)
dev.off()
rm(data)

#### Neg_T_2

data <- read.table("path/meth_rate/oxWGBS_Neg_T_2.txt", header = T)
pdf(file = "path/oxWGBS_Neg_T_2_meth_site_distribution.pdf", height = 4, width = 6)
CMplot(
  data, plot.type="d",  bin.size=1e6, col=c("darkgreen", "yellow", "red"),
  file="jpg", dpi=300, file.output=F, verbose=TRUE
)
dev.off()
rm(data)

#### Neg_T_3

data <- read.table("path/meth_rate/oxWGBS_Neg_T_3.txt", header = T)
pdf(file = "path/oxWGBS_Neg_T_3_meth_site_distribution.pdf", height = 4, width = 6)
CMplot(
  data, plot.type="d",  bin.size=1e6, col=c("darkgreen", "yellow", "red"),
  file="jpg", dpi=300, file.output=F, verbose=TRUE
)
dev.off()
rm(data)

#### Neg_T_4

data <- read.table("path/meth_rate/oxWGBS_Neg_T_4.txt", header = T)
pdf(file = "path/oxWGBS_Neg_T_4_meth_site_distribution.pdf", height = 4, width = 6)
CMplot(
  data, plot.type="d",  bin.size=1e6, col=c("darkgreen", "yellow", "red"),
  file="jpg", dpi=300, file.output=F, verbose=TRUE
)
dev.off()
rm(data)

### 5hmC

#### Pos_N

data <- read.table("path/meth_rate/oxWGBS_Pos_N.txt", header = T)
pdf(file = "path/oxWGBS_Pos_N_meth_site_distribution.pdf", height = 4, width = 6)
CMplot(
  data, plot.type="d",  bin.size=1e6, col=c("darkgreen", "yellow", "red"),
  file="jpg", dpi=300, file.output=F, verbose=TRUE
)
dev.off()
rm(data)

#### Pos_T

data <- read.table("path/meth_rate/oxWGBS_Pos_T.txt", header = T)
pdf(file = "path/oxWGBS_Pos_T_meth_site_distribution.pdf", height = 4, width = 6)
CMplot(
  data, plot.type="d",  bin.size=1e6, col=c("darkgreen", "yellow", "red"),
  file="jpg", dpi=300, file.output=F, verbose=TRUE
)
dev.off()
rm(data)

#### Pos_N_1

data <- read.table("path/meth_rate/oxWGBS_Pos_N_1.txt", header = T)
pdf(file = "path/oxWGBS_Pos_N_1_meth_site_distribution.pdf", height = 4, width = 6)
CMplot(
  data, plot.type="d",  bin.size=1e6, col=c("darkgreen", "yellow", "red"),
  file="jpg", dpi=300, file.output=F, verbose=TRUE
)
dev.off()
rm(data)

#### Pos_N_2

data <- read.table("path/meth_rate/oxWGBS_Pos_N_2.txt", header = T)
pdf(file = "path/oxWGBS_Pos_N_2_meth_site_distribution.pdf", height = 4, width = 6)
CMplot(
  data, plot.type="d",  bin.size=1e6, col=c("darkgreen", "yellow", "red"),
  file="jpg", dpi=300, file.output=F, verbose=TRUE
)
dev.off()
rm(data)

#### Pos_N_3

data <- read.table("path/meth_rate/oxWGBS_Pos_N_3.txt", header = T)
pdf(file = "path/oxWGBS_Pos_N_3_meth_site_distribution.pdf", height = 4, width = 6)
CMplot(
  data, plot.type="d",  bin.size=1e6, col=c("darkgreen", "yellow", "red"),
  file="jpg", dpi=300, file.output=F, verbose=TRUE
)
dev.off()
rm(data)

#### Pos_N_4

data <- read.table("path/meth_rate/oxWGBS_Pos_N_4.txt", header = T)
pdf(file = "path/oxWGBS_Pos_N_4_meth_site_distribution.pdf", height = 4, width = 6)
CMplot(
  data, plot.type="d",  bin.size=1e6, col=c("darkgreen", "yellow", "red"),
  file="jpg", dpi=300, file.output=F, verbose=TRUE
)
dev.off()
rm(data)

#### Pos_T_1

data <- read.table("path/meth_rate/oxWGBS_Pos_T_1.txt", header = T)
pdf(file = "path/oxWGBS_Pos_T_1_meth_site_distribution.pdf", height = 4, width = 6)
CMplot(
  data, plot.type="d",  bin.size=1e6, col=c("darkgreen", "yellow", "red"),
  file="jpg", dpi=300, file.output=F, verbose=TRUE
)
dev.off()
rm(data)

#### Pos_T_2

data <- read.table("path/meth_rate/oxWGBS_Pos_T_2.txt", header = T)
pdf(file = "path/oxWGBS_Pos_T_2_meth_site_distribution.pdf", height = 4, width = 6)
CMplot(
  data, plot.type="d",  bin.size=1e6, col=c("darkgreen", "yellow", "red"),
  file="jpg", dpi=300, file.output=F, verbose=TRUE
)
dev.off()
rm(data)

#### Pos_T_3

data <- read.table("path/meth_rate/oxWGBS_Pos_T_3.txt", header = T)
pdf(file = "path/oxWGBS_Pos_T_3_meth_site_distribution.pdf", height = 4, width = 6)
CMplot(
  data, plot.type="d",  bin.size=1e6, col=c("darkgreen", "yellow", "red"),
  file="jpg", dpi=300, file.output=F, verbose=TRUE
)
dev.off()
rm(data)

#### Pos_T_4

data <- read.table("path/meth_rate/oxWGBS_Pos_T_4.txt", header = T)
pdf(file = "path/oxWGBS_Pos_T_4_meth_site_distribution.pdf", height = 4, width = 6)
CMplot(
  data, plot.type="d",  bin.size=1e6, col=c("darkgreen", "yellow", "red"),
  file="jpg", dpi=300, file.output=F, verbose=TRUE
)
dev.off()
rm(data)


## DNA methylation levels of genes

library(ggplot2)
library(tidyverse)
library(reshape2)
library(ggsci)
library(Rmisc)
library(data.table)

### WGBS

#### All meth
scaler <- fread("path/meth_distribution/combine-scaleRegion-data-WGBS.gz",header = F)

# check
head(scaler[1:3,1:8])
        

clean_scaler <- scaler %>% select(c(-1:-3,-5,-6))
head(clean_scaler[1:3,1:8])
# wide to long
long_scaler <- melt(clean_scaler,id.vars = 'V4',value.name = 'signal') %>%
  select(-variable)

# add sample name
long_scaler$sample <- rep(c("Neg_N","Neg_T", "Pos_N", "Pos_T"),
                          each = nrow(scaler)*220)

# add x position
long_scaler$pos <- rep(c(1:220),each = nrow(scaler),times = 4)

# calculate means
filnal_scaler <- long_scaler %>%
  # remove na
  tidyr::drop_na() %>%
  dplyr::group_by(sample,pos) %>%
  # mean and 95% interval confidence
  dplyr::summarise(mean_signal = mean(signal),
                   sd = sd(signal),
                   upper = CI(signal,ci = 0.95)[1],
                   lower = CI(signal,ci = 0.95)[3])

# check
head(filnal_scaler)

pdf(file = "path/WGBS_methylation_levels_of_genes.pdf", height = 3, width = 4)
ggplot(filnal_scaler,aes(x = pos,y = mean_signal)) +
  geom_line(aes(color = sample),size = 1) +
  theme_classic(base_size = 14) +
  scale_color_manual(values = mycolour) +
  # x label
  scale_x_continuous(breaks = c(0,60,160,220),
                     labels = c('-3 kb','TSS','TES','+3 kb')) +
  xlab('') + ylab('Normalized signal') +
  theme(aspect.ratio = 0.8)
dev.off()

#### All meth
scaler <- fread("path/meth_distribution/combine-scaleRegion-data-WGBS-sample.gz",header = F)

# check
head(scaler[1:3,1:8])


clean_scaler <- scaler %>% select(c(-1:-3,-5,-6))
head(clean_scaler[1:3,1:8])
# wide to long
long_scaler <- melt(clean_scaler,id.vars = 'V4',value.name = 'signal') %>%
  select(-variable)

# add sample name
long_scaler$sample <- rep(c("Neg_N_1", "Neg_N_2", "Neg_N_3", "Neg_N_4", "Neg_T_1", "Neg_T_2", "Neg_T_3", "Neg_T_4", "Pos_N_1", "Pos_N_2", "Pos_N_3", "Pos_N_4", "Pos_T_1", "Pos_T_2", "Pos_T_3", "Pos_T_4"),
                          each = nrow(scaler)*220)

# add x position
long_scaler$pos <- rep(c(1:220),each = nrow(scaler),times = 16)

# calculate means
filnal_scaler <- long_scaler %>%
  # remove na
  tidyr::drop_na() %>%
  dplyr::group_by(sample,pos) %>%
  # mean and 95% interval confidence
  dplyr::summarise(mean_signal = mean(signal),
                   sd = sd(signal),
                   upper = CI(signal,ci = 0.95)[1],
                   lower = CI(signal,ci = 0.95)[3])

# check
head(filnal_scaler)

pdf(file = "path/WGBS-sample_methylation_levels_of_genes.pdf", height = 3, width = 5)
ggplot(filnal_scaler,aes(x = pos,y = mean_signal)) +
  geom_line(aes(color = sample),size = 1) +
  theme_classic(base_size = 14) +
  scale_color_manual(values = mycolour) +
  # x label
  scale_x_continuous(breaks = c(0,60,160,220),
                     labels = c('-3 kb','TSS','TES','+3 kb')) +
  xlab('') + ylab('Normalized signal') +
  theme(aspect.ratio = 0.8)
dev.off()
                              
### oxWGBS

#### All meth
scaler <- fread("path/meth_distribution/combine-scaleRegion-data-oxWGBS.gz",header = F)

# check
head(scaler[1:3,1:8])


clean_scaler <- scaler %>% select(c(-1:-3,-5,-6))
head(clean_scaler[1:3,1:8])
# wide to long
long_scaler <- melt(clean_scaler,id.vars = 'V4',value.name = 'signal') %>%
  select(-variable)

# add sample name
long_scaler$sample <- rep(c("Neg_N","Neg_T", "Pos_N", "Pos_T"),
                          each = nrow(scaler)*220)

# add x position
long_scaler$pos <- rep(c(1:220),each = nrow(scaler),times = 4)

# calculate means
filnal_scaler <- long_scaler %>%
  # remove na
  tidyr::drop_na() %>%
  dplyr::group_by(sample,pos) %>%
  # mean and 95% interval confidence
  dplyr::summarise(mean_signal = mean(signal),
                   sd = sd(signal),
                   upper = CI(signal,ci = 0.95)[1],
                   lower = CI(signal,ci = 0.95)[3])

# check
head(filnal_scaler)

pdf(file = "path/oxWGBS_methylation_levels_of_genes.pdf", height = 3, width = 4)
ggplot(filnal_scaler,aes(x = pos,y = mean_signal)) +
  geom_line(aes(color = sample),size = 1) +
  theme_classic(base_size = 14) +
  scale_color_manual(values = mycolour) +
  # x label
  scale_x_continuous(breaks = c(0,60,160,220),
                     labels = c('-3 kb','TSS','TES','+3 kb')) +
  xlab('') + ylab('Normalized signal') +
  theme(aspect.ratio = 0.8)
dev.off()

#### All meth
scaler <- fread("path/meth_distribution/combine-scaleRegion-data-oxWGBS-sample.gz",header = F)

# check
head(scaler[1:3,1:8])


clean_scaler <- scaler %>% select(c(-1:-3,-5,-6))
head(clean_scaler[1:3,1:8])
# wide to long
long_scaler <- melt(clean_scaler,id.vars = 'V4',value.name = 'signal') %>%
  select(-variable)

# add sample name
long_scaler$sample <- rep(c("Neg_N_1", "Neg_N_2", "Neg_N_3", "Neg_N_4", "Neg_T_1", "Neg_T_2", "Neg_T_3", "Neg_T_4", "Pos_N_1", "Pos_N_2", "Pos_N_3", "Pos_N_4", "Pos_T_1", "Pos_T_2", "Pos_T_3", "Pos_T_4"),
                          each = nrow(scaler)*220)

# add x position
long_scaler$pos <- rep(c(1:220),each = nrow(scaler),times = 16)

# calculate means
filnal_scaler <- long_scaler %>%
  # remove na
  tidyr::drop_na() %>%
  dplyr::group_by(sample,pos) %>%
  # mean and 95% interval confidence
  dplyr::summarise(mean_signal = mean(signal),
                   sd = sd(signal),
                   upper = CI(signal,ci = 0.95)[1],
                   lower = CI(signal,ci = 0.95)[3])

# check
head(filnal_scaler)

pdf(file = "path/oxWGBS-sample_methylation_levels_of_genes.pdf", height = 3, width = 5)
ggplot(filnal_scaler,aes(x = pos,y = mean_signal)) +
  geom_line(aes(color = sample),size = 1) +
  theme_classic(base_size = 14) +
  scale_color_manual(values = mycolour) +
  # x label
  scale_x_continuous(breaks = c(0,60,160,220),
                     labels = c('-3 kb','TSS','TES','+3 kb')) +
  xlab('') + ylab('Normalized signal') +
  theme(aspect.ratio = 0.8)
dev.off()

## DNA methylation levels of TE

library(ggplot2)
library(tidyverse)
library(reshape2)
library(ggsci)
library(Rmisc)
library(data.table)

### WGBS

#### All meth
scaler <- fread("path/meth_distribution/TE/TE_combine-scaleRegion-data-WGBS.gz",header = F)

# check
head(scaler[1:3,1:8])


clean_scaler <- scaler %>% select(c(-1:-3,-5,-6))
head(clean_scaler[1:3,1:8])
# wide to long
long_scaler <- melt(clean_scaler,id.vars = 'V4',value.name = 'signal') %>%
  select(-variable)

# add sample name
long_scaler$sample <- rep(c("Neg_N","Neg_T", "Pos_N", "Pos_T"),
                          each = nrow(scaler)*220)

# add x position
long_scaler$pos <- rep(c(1:220),each = nrow(scaler),times = 4)

# calculate means
filnal_scaler <- long_scaler %>%
  # remove na
  tidyr::drop_na() %>%
  dplyr::group_by(sample,pos) %>%
  # mean and 95% interval confidence
  dplyr::summarise(mean_signal = mean(signal),
                   sd = sd(signal),
                   upper = CI(signal,ci = 0.95)[1],
                   lower = CI(signal,ci = 0.95)[3])

# check
head(filnal_scaler)

pdf(file = "path/WGBS_methylation_levels_of_TE.pdf", height = 3, width = 4)
ggplot(filnal_scaler,aes(x = pos,y = mean_signal)) +
  geom_line(aes(color = sample),size = 1) +
  theme_classic(base_size = 14) +
  scale_color_manual(values = mycolour) +
  # x label
  scale_x_continuous(breaks = c(0,60,160,220),
                     labels = c('-3 kb','TSS','TES','+3 kb')) +
  xlab('') + ylab('Normalized signal') +
  theme(aspect.ratio = 0.8)
dev.off()

#### All meth
scaler <- fread("path/meth_distribution/TE/TE_combine-scaleRegion-data-WGBS-sample.gz",header = F)

# check
head(scaler[1:3,1:8])


clean_scaler <- scaler %>% select(c(-1:-3,-5,-6))
head(clean_scaler[1:3,1:8])
# wide to long
long_scaler <- melt(clean_scaler,id.vars = 'V4',value.name = 'signal') %>%
  select(-variable)

# add sample name
long_scaler$sample <- rep(c("Neg_N_1", "Neg_N_2", "Neg_N_3", "Neg_N_4", "Neg_T_1", "Neg_T_2", "Neg_T_3", "Neg_T_4", "Pos_N_1", "Pos_N_2", "Pos_N_3", "Pos_N_4", "Pos_T_1", "Pos_T_2", "Pos_T_3", "Pos_T_4"),
                          each = nrow(scaler)*220)

# add x position
long_scaler$pos <- rep(c(1:220),each = nrow(scaler),times = 16)

# calculate means
filnal_scaler <- long_scaler %>%
  # remove na
  tidyr::drop_na() %>%
  dplyr::group_by(sample,pos) %>%
  # mean and 95% interval confidence
  dplyr::summarise(mean_signal = mean(signal),
                   sd = sd(signal),
                   upper = CI(signal,ci = 0.95)[1],
                   lower = CI(signal,ci = 0.95)[3])

# check
head(filnal_scaler)

pdf(file = "path/WGBS-sample_methylation_levels_of_TE.pdf", height = 3, width = 5)
ggplot(filnal_scaler,aes(x = pos,y = mean_signal)) +
  geom_line(aes(color = sample),size = 1) +
  theme_classic(base_size = 14) +
  scale_color_manual(values = mycolour) +
  # x label
  scale_x_continuous(breaks = c(0,60,160,220),
                     labels = c('-3 kb','TSS','TES','+3 kb')) +
  xlab('') + ylab('Normalized signal') +
  theme(aspect.ratio = 0.8)
dev.off()

### oxWGBS

#### All meth
scaler <- fread("path/meth_distribution/TE/TE_combine-scaleRegion-data-oxWGBS.gz",header = F)

# check
head(scaler[1:3,1:8])


clean_scaler <- scaler %>% select(c(-1:-3,-5,-6))
head(clean_scaler[1:3,1:8])
# wide to long
long_scaler <- melt(clean_scaler,id.vars = 'V4',value.name = 'signal') %>%
  select(-variable)

# add sample name
long_scaler$sample <- rep(c("Neg_N","Neg_T", "Pos_N", "Pos_T"),
                          each = nrow(scaler)*220)

# add x position
long_scaler$pos <- rep(c(1:220),each = nrow(scaler),times = 4)

# calculate means
filnal_scaler <- long_scaler %>%
  # remove na
  tidyr::drop_na() %>%
  dplyr::group_by(sample,pos) %>%
  # mean and 95% interval confidence
  dplyr::summarise(mean_signal = mean(signal),
                   sd = sd(signal),
                   upper = CI(signal,ci = 0.95)[1],
                   lower = CI(signal,ci = 0.95)[3])

# check
head(filnal_scaler)

pdf(file = "path/oxWGBS_methylation_levels_of_TE.pdf", height = 3, width = 4)
ggplot(filnal_scaler,aes(x = pos,y = mean_signal)) +
  geom_line(aes(color = sample),size = 1) +
  theme_classic(base_size = 14) +
  scale_color_manual(values = mycolour) +
  # x label
  scale_x_continuous(breaks = c(0,60,160,220),
                     labels = c('-3 kb','TSS','TES','+3 kb')) +
  xlab('') + ylab('Normalized signal') +
  theme(aspect.ratio = 0.8)
dev.off()

#### All meth
scaler <- fread("path/meth_distribution/TE/TE_combine-scaleRegion-data-oxWGBS-sample.gz",header = F)

# check
head(scaler[1:3,1:8])


clean_scaler <- scaler %>% select(c(-1:-3,-5,-6))
head(clean_scaler[1:3,1:8])
# wide to long
long_scaler <- melt(clean_scaler,id.vars = 'V4',value.name = 'signal') %>%
  select(-variable)

# add sample name
long_scaler$sample <- rep(c("Neg_N_1", "Neg_N_2", "Neg_N_3", "Neg_N_4", "Neg_T_1", "Neg_T_2", "Neg_T_3", "Neg_T_4", "Pos_N_1", "Pos_N_2", "Pos_N_3", "Pos_N_4", "Pos_T_1", "Pos_T_2", "Pos_T_3", "Pos_T_4"),
                          each = nrow(scaler)*220)

# add x position
long_scaler$pos <- rep(c(1:220),each = nrow(scaler),times = 16)

# calculate means
filnal_scaler <- long_scaler %>%
  # remove na
  tidyr::drop_na() %>%
  dplyr::group_by(sample,pos) %>%
  # mean and 95% interval confidence
  dplyr::summarise(mean_signal = mean(signal),
                   sd = sd(signal),
                   upper = CI(signal,ci = 0.95)[1],
                   lower = CI(signal,ci = 0.95)[3])

# check
head(filnal_scaler)

pdf(file = "path/oxWGBS-sample_methylation_levels_of_TE.pdf", height = 3, width = 5)
ggplot(filnal_scaler,aes(x = pos,y = mean_signal)) +
  geom_line(aes(color = sample),size = 1) +
  theme_classic(base_size = 14) +
  scale_color_manual(values = mycolour) +
  # x label
  scale_x_continuous(breaks = c(0,60,160,220),
                     labels = c('-3 kb','TSS','TES','+3 kb')) +
  xlab('') + ylab('Normalized signal') +
  theme(aspect.ratio = 0.8)
dev.off()

## genebody distribution

### WGBS

#### group level

WGBS_species <- read.table("path/meth_genebody_distribution/WGBS_genebody.txt", header = T)
WGBS_species$sample <- factor(WGBS_species$sample,levels =c("Neg_N","Neg_T","Pos_N","Pos_T"))
WGBS_species$group <- factor(WGBS_species$group,levels =c("CDS_Exons","5'UTR_Exons","3'UTR_Exons","Introns", "Intergenic"))
pdf(file = "path/WGBS_genebody_distribution.pdf", width = 3, height = 3)
ggplot(WGBS_species, aes(x = sample, y = rate, fill = group))+
  geom_bar(stat="identity")+
  labs(x = "", y = "Percentage (%)", title = "")+
  theme_bw()+
  theme()+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(values = c('#af2337', '#ecc342', '#2967a0', '#96b437','#2f3c28'))+
  theme_bw()+theme(panel.grid =element_blank())+ theme(legend.text.align = 0)+
  theme(axis.text.x = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5, angle = 90))+
  theme(axis.text.y = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(axis.title.y = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(legend.text = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(legend.title = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5))
dev.off()

oxWGBS_species <- read.table("path/meth_genebody_distribution/oxWGBS_genebody.txt", header = T)
oxWGBS_species$sample <- factor(oxWGBS_species$sample,levels =c("Neg_N","Neg_T","Pos_N","Pos_T"))
oxWGBS_species$group <- factor(oxWGBS_species$group,levels =c("CDS_Exons","5'UTR_Exons","3'UTR_Exons","Introns", "Intergenic"))
pdf(file = "path/oxWGBS_genebody_distribution.pdf", width = 3, height = 3)
ggplot(oxWGBS_species, aes(x = sample, y = rate, fill = group))+
  geom_bar(stat="identity")+
  labs(x = "", y = "Percentage (%)", title = "")+
  theme_bw()+
  theme()+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(values = c('#af2337', '#ecc342', '#2967a0', '#96b437','#2f3c28'))+
  theme_bw()+theme(panel.grid =element_blank())+ theme(legend.text.align = 0)+
  theme(axis.text.x = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5, angle = 90))+
  theme(axis.text.y = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(axis.title.y = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(legend.text = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(legend.title = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5))
dev.off()

#### Sample level

WGBS_species <- read.table("path/meth_genebody_distribution/WGBS_genebody_sample.txt", header = T)
WGBS_species$sample <- factor(WGBS_species$sample,levels =c("Neg_1_N", "Neg_2_N", "Neg_3_N", "Neg_4_N",
                                                            "Neg_1_T", "Neg_2_T", "Neg_3_T", "Neg_4_T",
                                                            "Pos_1_N", "Pos_2_N", "Pos_3_N", "Pos_4_N",
                                                            "Pos_1_T", "Pos_2_T", "Pos_3_T", "Pos_4_T"))
WGBS_species$group <- factor(WGBS_species$group,levels =c("CDS_Exons","5'UTR_Exons","3'UTR_Exons","Introns", "Intergenic"))
pdf(file = "path/WGBS_genebody_sample_distribution.pdf", width = 5, height = 3)
ggplot(WGBS_species, aes(x = sample, y = rate, fill = group))+
  geom_bar(stat="identity")+
  labs(x = "", y = "Percentage (%)", title = "")+
  theme_bw()+
  theme()+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(values = c('#af2337', '#ecc342', '#2967a0', '#96b437','#2f3c28'))+
  theme_bw()+theme(panel.grid =element_blank())+ theme(legend.text.align = 0)+
  theme(axis.text.x = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5, angle = 90))+
  theme(axis.text.y = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(axis.title.y = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(legend.text = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(legend.title = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5))
dev.off()

oxWGBS_species <- read.table("path/meth_genebody_distribution/oxWGBS_genebody_sample.txt", header = T)
oxWGBS_species$sample <- factor(oxWGBS_species$sample,levels =c("Neg_1_N", "Neg_2_N", "Neg_3_N", "Neg_4_N",
                                                                "Neg_1_T", "Neg_2_T", "Neg_3_T", "Neg_4_T",
                                                                "Pos_1_N", "Pos_2_N", "Pos_3_N", "Pos_4_N",
                                                                "Pos_1_T", "Pos_2_T", "Pos_3_T", "Pos_4_T"))
oxWGBS_species$group <- factor(oxWGBS_species$group,levels =c("CDS_Exons","5'UTR_Exons","3'UTR_Exons","Introns", "Intergenic"))
pdf(file = "path/oxWGBS_genebody_sample_distribution.pdf", width = 5, height = 3)
ggplot(oxWGBS_species, aes(x = sample, y = rate, fill = group))+
  geom_bar(stat="identity")+
  labs(x = "", y = "Percentage (%)", title = "")+
  theme_bw()+
  theme()+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(values = c('#af2337', '#ecc342', '#2967a0', '#96b437','#2f3c28'))+
  theme_bw()+theme(panel.grid =element_blank())+ theme(legend.text.align = 0)+
  theme(axis.text.x = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5, angle = 90))+
  theme(axis.text.y = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(axis.title.y = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(legend.text = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(legend.title = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5))
dev.off()

# DMR and DhMR

## DMR and DhMR distribution

### PosN_vs_PosT

WGBS_species <- read.table("path/PosN_vs_PosT_DMR_genebody_distribution.txt", header = T)
WGBS_species$sample <- factor(WGBS_species$sample,levels =c("DMR", "DhMR"))
WGBS_species$group <- factor(WGBS_species$group,levels =c("CDS_Exons","5'UTR_Exons","3'UTR_Exons","Introns", "Intergenic"))
pdf(file = "path/Part2/PosN_vs_PosT_DMR_genebody_distribution.pdf", width = 3, height = 3)
ggplot(WGBS_species, aes(x = sample, y = rate, fill = group))+
  geom_bar(stat="identity")+
  labs(x = "", y = "Percentage (%)", title = "")+
  theme_bw()+
  theme()+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(values = c('#af2337', '#ecc342', '#2967a0', '#96b437','#2f3c28'))+
  theme_bw()+theme(panel.grid =element_blank())+ theme(legend.text.align = 0)+
  theme(axis.text.x = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5, angle = 90))+
  theme(axis.text.y = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(axis.title.y = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(legend.text = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(legend.title = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5))
dev.off()

### NegN_vs_NegT

WGBS_species <- read.table("path/NegN_vs_NegT_DMR_genebody_distribution.txt", header = T)
WGBS_species$sample <- factor(WGBS_species$sample,levels =c("DMR", "DhMR"))
WGBS_species$group <- factor(WGBS_species$group,levels =c("CDS_Exons","5'UTR_Exons","3'UTR_Exons","Introns", "Intergenic"))
pdf(file = "path/Part2/NegN_vs_NegT_DMR_genebody_distribution.pdf", width = 3, height = 3)
ggplot(WGBS_species, aes(x = sample, y = rate, fill = group))+
  geom_bar(stat="identity")+
  labs(x = "", y = "Percentage (%)", title = "")+
  theme_bw()+
  theme()+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(values = c('#af2337', '#ecc342', '#2967a0', '#96b437','#2f3c28'))+
  theme_bw()+theme(panel.grid =element_blank())+ theme(legend.text.align = 0)+
  theme(axis.text.x = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5, angle = 90))+
  theme(axis.text.y = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(axis.title.y = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(legend.text = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(legend.title = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5))
dev.off()

### NegT_vs_PosT

WGBS_species <- read.table("path/NegT_vs_PosT_DMR_genebody_distribution.txt", header = T)
WGBS_species$sample <- factor(WGBS_species$sample,levels =c("DMR", "DhMR"))
WGBS_species$group <- factor(WGBS_species$group,levels =c("CDS_Exons","5'UTR_Exons","3'UTR_Exons","Introns", "Intergenic"))
pdf(file = "path/Part2/NegT_vs_PosT_DMR_genebody_distribution.pdf", width = 3, height = 3)
ggplot(WGBS_species, aes(x = sample, y = rate, fill = group))+
  geom_bar(stat="identity")+
  labs(x = "", y = "Percentage (%)", title = "")+
  theme_bw()+
  theme()+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(values = c('#af2337', '#ecc342', '#2967a0', '#96b437','#2f3c28'))+
  theme_bw()+theme(panel.grid =element_blank())+ theme(legend.text.align = 0)+
  theme(axis.text.x = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5, angle = 90))+
  theme(axis.text.y = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(axis.title.y = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(legend.text = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(legend.title = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5))
dev.off()

### NegN_vs_PosN

WGBS_species <- read.table("path/NegN_vs_PosN_DMR_genebody_distribution.txt", header = T)
WGBS_species$sample <- factor(WGBS_species$sample,levels =c("DMR", "DhMR"))
WGBS_species$group <- factor(WGBS_species$group,levels =c("CDS_Exons","5'UTR_Exons","3'UTR_Exons","Introns", "Intergenic"))
pdf(file = "path/Part2/NegN_vs_PosN_DMR_genebody_distribution.pdf", width = 3, height = 3)
ggplot(WGBS_species, aes(x = sample, y = rate, fill = group))+
  geom_bar(stat="identity")+
  labs(x = "", y = "Percentage (%)", title = "")+
  theme_bw()+
  theme()+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(values = c('#af2337', '#ecc342', '#2967a0', '#96b437','#2f3c28'))+
  theme_bw()+theme(panel.grid =element_blank())+ theme(legend.text.align = 0)+
  theme(axis.text.x = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5, angle = 90))+
  theme(axis.text.y = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(axis.title.y = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(legend.text = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(legend.title = element_text(size = 10, face = "bold", vjust = 0.5, hjust = 0.5))
dev.off()


## DMR and DhMR annotation

### GO and KEGG


library("ChIPseeker")
library("GenomicFeatures")
ara_TxDb <- makeTxDbFromGFF("/home/yangwl/data/reference/hg38/Homo_sapiens.GRCh38.109.chr_patch_hapl_scaff.gtf")
library("clusterProfiler")
library("org.Hs.eg.db")

#### DMR

##### PosN_vs_PosT
files_DMR = list(DMR = "/home/yangwl/data/HCC_bulk/result/WGBS/Pos/circos/WGBS_PosN_vs_PosT_DMR_peaks.narrowPeak")
peakAnno_DMR <- lapply(files_DMR, annotatePeak, 
                      TxDb=ara_TxDb,tssRegion=c(-3000, 3000))

gene = lapply(peakAnno_DMR, function(i) as.data.frame(i)$geneId)

gene_2 <- unlist(gene, use.names = FALSE )
write.table(gene_2, file = "/home/yangwl/data/HCC_bulk/result/WGBS/Pos/circos/WGBS_PosN_vs_PosT_DMR_associated_genes.txt", col.names=FALSE, row.names=FALSE, quote=FALSE)

##### NegN_vs_NegT
files_DMR = list(DMR = "/home/yangwl/data/HCC_bulk/result/WGBS/Neg/circos/WGBS_NegN_vs_NegT_DMR_peaks.narrowPeak")
peakAnno_DMR <- lapply(files_DMR, annotatePeak, 
                       TxDb=ara_TxDb,tssRegion=c(-3000, 3000))

gene = lapply(peakAnno_DMR, function(i) as.data.frame(i)$geneId)

gene_2 <- unlist(gene, use.names = FALSE )
write.table(gene_2, file = "/home/yangwl/data/HCC_bulk/result/WGBS/Neg/circos/WGBS_NegN_vs_NegT_DMR_associated_genes.txt", col.names=FALSE, row.names=FALSE, quote=FALSE)

##### NegT_vs_PosT
files_DMR = list(DMR = "/home/yangwl/data/HCC_bulk/result/WGBS/Neg_T_vs_Pos_T/cricos/WGBS_NegT_vs_PosT_DMR_peaks.narrowPeak")
peakAnno_DMR <- lapply(files_DMR, annotatePeak, 
                       TxDb=ara_TxDb,tssRegion=c(-3000, 3000))

gene = lapply(peakAnno_DMR, function(i) as.data.frame(i)$geneId)

gene_2 <- unlist(gene, use.names = FALSE )
write.table(gene_2, file = "/home/yangwl/data/HCC_bulk/result/WGBS/Neg_T_vs_Pos_T/cricos/WGBS_NegT_vs_PosT_DMR_associated_genes.txt", col.names=FALSE, row.names=FALSE, quote=FALSE)

##### NegN_vs_PosN
files_DMR = list(DMR = "/home/yangwl/data/HCC_bulk/result/WGBS/Neg_N_vs_Pos_N/circos/WGBS_NegN_vs_PosN_DMR_peaks.narrowPeak")
peakAnno_DMR <- lapply(files_DMR, annotatePeak, 
                       TxDb=ara_TxDb,tssRegion=c(-3000, 3000))

gene = lapply(peakAnno_DMR, function(i) as.data.frame(i)$geneId)

gene_2 <- unlist(gene, use.names = FALSE )
write.table(gene_2, file = "/home/yangwl/data/HCC_bulk/result/WGBS/Neg_N_vs_Pos_N/circos/WGBS_NegN_vs_PosN_DMR_associated_genes.txt", col.names=FALSE, row.names=FALSE, quote=FALSE)


#### DhMR
##### NegN_vs_PosT
files_DMR = list(DMR = "/home/yangwl/data/HCC_bulk/result/oxWGBS/Pos/circos/oxWGBS_PosN_vs_PosT_DMR_peaks.narrowPeak")
peakAnno_DMR <- lapply(files_DMR, annotatePeak, 
                       TxDb=ara_TxDb,tssRegion=c(-3000, 3000))

gene = lapply(peakAnno_DMR, function(i) as.data.frame(i)$geneId)

gene_2 <- unlist(gene, use.names = FALSE )
write.table(gene_2, file = "/home/yangwl/data/HCC_bulk/result/oxWGBS/Pos/circos/oxWGBS_PosN_vs_PosT_DMR_associated_genes.txt", col.names=FALSE, row.names=FALSE, quote=FALSE)

##### NegN_vs_NegT
files_DMR = list(DMR = "/home/yangwl/data/HCC_bulk/result/oxWGBS/Neg/circos/oxWGBS_NegN_vs_NegT_DMR_peaks.narrowPeak")
peakAnno_DMR <- lapply(files_DMR, annotatePeak, 
                       TxDb=ara_TxDb,tssRegion=c(-3000, 3000))

gene = lapply(peakAnno_DMR, function(i) as.data.frame(i)$geneId)

gene_2 <- unlist(gene, use.names = FALSE )
write.table(gene_2, file = "/home/yangwl/data/HCC_bulk/result/oxWGBS/Neg/circos/oxWGBS_NegN_vs_NegT_DMR_associated_genes.txt", col.names=FALSE, row.names=FALSE, quote=FALSE)

##### NegT_vs_PosT
files_DMR = list(DMR = "/home/yangwl/data/HCC_bulk/result/oxWGBS/Neg_T_vs_Pos_T/cricos/oxWGBS_NegT_vs_PosT_DMR_peaks.narrowPeak")
peakAnno_DMR <- lapply(files_DMR, annotatePeak, 
                       TxDb=ara_TxDb,tssRegion=c(-3000, 3000))

gene = lapply(peakAnno_DMR, function(i) as.data.frame(i)$geneId)

gene_2 <- unlist(gene, use.names = FALSE )
write.table(gene_2, file = "/home/yangwl/data/HCC_bulk/result/oxWGBS/Neg_T_vs_Pos_T/cricos/oxWGBS_NegT_vs_PosT_DMR_associated_genes.txt", col.names=FALSE, row.names=FALSE, quote=FALSE)

##### NegN_vs_PosN
files_DMR = list(DMR = "/home/yangwl/data/HCC_bulk/result/oxWGBS/Neg_N_vs_Pos_N/circos/oxWGBS_NegN_vs_PosN_DMR_peaks.narrowPeak")
peakAnno_DMR <- lapply(files_DMR, annotatePeak, 
                       TxDb=ara_TxDb,tssRegion=c(-3000, 3000))

gene = lapply(peakAnno_DMR, function(i) as.data.frame(i)$geneId)

gene_2 <- unlist(gene, use.names = FALSE )
write.table(gene_2, file = "/home/yangwl/data/HCC_bulk/result/oxWGBS/Neg_N_vs_Pos_N/circos/oxWGBS_NegN_vs_PosN_DMR_associated_genes.txt", col.names=FALSE, row.names=FALSE, quote=FALSE)


# TCGA The expression of DMR or DhMR associated genes (meth)

library(TCGAplot)
gene <- c("ARL6IP4",
          "CCDC200",
          "DHDH",
          "DUX4L50",
          "ESPNP",
          "FCGBP",
          "KCNQ3",
          "KCNS1",
          "LRATD2",
          "MEX3B",
          "MIR1257",
          "MIR9-1HG",
          "MTND3P22",
          "OGFOD3",
          "POU4F1",
          "PSMA8",
          "PSME2P6",
          "RARA",
          "RNA5-8SN1",
          "RNA5-8SN2",
          "RNA5-8SN3",
          "SEL1L3",
          "SLC22A23",
          "SPATA18",
          "TEC",
          "TNNT3",
          "WDR97",
          "ZNF285")
for (i in gene) { gene_methylation_scatter("LIHC", i)
  methy_kmplot("LIHC", i)
}
pdf(file = "path/Part3/OS_ARL6IP4.pdf", height = 6, width = 6)
tcga_kmplot("LIHC", "ARL6IP4")
dev.off()
pdf(file = "path/Part3/OS_CCDC200.pdf", height = 6, width = 6)
tcga_kmplot("LIHC", "CCDC200")
dev.off()
pdf(file = "path/Part3/OS_DHDH.pdf", height = 6, width = 6)
tcga_kmplot("LIHC", "DHDH")
dev.off()
pdf(file = "path/Part3/OS_DUX4L50.pdf", height = 6, width = 6)
tcga_kmplot("LIHC", "DUX4L50")
dev.off()
pdf(file = "path/Part3/OS_ESPNP.pdf", height = 6, width = 6)
tcga_kmplot("LIHC", "ESPNP")
dev.off()
pdf(file = "path/Part3/OS_FCGBP.pdf", height = 6, width = 6)
tcga_kmplot("LIHC", "FCGBP")
dev.off()
pdf(file = "path/Part3/OS_KCNQ3.pdf", height = 6, width = 6)
tcga_kmplot("LIHC", "KCNQ3")
dev.off()
pdf(file = "path/Part3/OS_KCNS1.pdf", height = 6, width = 6)
tcga_kmplot("LIHC", "KCNS1")
dev.off()
pdf(file = "path/Part3/OS_LRATD2.pdf", height = 6, width = 6)
tcga_kmplot("LIHC", "LRATD2")
dev.off()
pdf(file = "path/Part3/OS_MEX3B.pdf", height = 6, width = 6)
tcga_kmplot("LIHC", "MEX3B")
dev.off()
pdf(file = "path/Part3/OS_MIR1257.pdf", height = 6, width = 6)
tcga_kmplot("LIHC", "MIR1257")
dev.off()
pdf(file = "path/Part3/OS_MIR9-1HG.pdf", height = 6, width = 6)
tcga_kmplot("LIHC", "MIR9-1HG")
dev.off()
pdf(file = "path/Part3/OS_MTND3P22.pdf", height = 6, width = 6)
tcga_kmplot("LIHC", "MTND3P22")
dev.off()
pdf(file = "path/Part3/OS_OGFOD3.pdf", height = 6, width = 6)
tcga_kmplot("LIHC", "OGFOD3")
dev.off()
pdf(file = "path/Part3/OS_POU4F1.pdf", height = 6, width = 6)
tcga_kmplot("LIHC", "POU4F1")
dev.off()
pdf(file = "path/Part3/OS_PSMA8.pdf", height = 6, width = 6)
tcga_kmplot("LIHC", "PSMA8")
dev.off()
pdf(file = "path/Part3/OS_PSME2P6.pdf", height = 6, width = 6)
tcga_kmplot("LIHC", "PSME2P6")
dev.off()
pdf(file = "path/Part3/OS_RARA.pdf", height = 6, width = 6)
tcga_kmplot("LIHC", "RARA")
dev.off()
pdf(file = "path/Part3/OS_RNA5-8SN1.pdf", height = 6, width = 6)
tcga_kmplot("LIHC", "RNA5-8SN1")
dev.off()
pdf(file = "path/Part3/OS_RNA5-8SN2.pdf", height = 6, width = 6)
tcga_kmplot("LIHC", "RNA5-8SN2")
dev.off()
pdf(file = "path/Part3/OS_RNA5-8SN3.pdf", height = 6, width = 6)
tcga_kmplot("LIHC", "RNA5-8SN3")
dev.off()
pdf(file = "path/Part3/OS_SEL1L3.pdf", height = 6, width = 6)
tcga_kmplot("LIHC", "SEL1L3")
dev.off()
pdf(file = "path/Part3/OS_SLC22A23.pdf", height = 6, width = 6)
tcga_kmplot("LIHC", "SLC22A23")
dev.off()
pdf(file = "path/Part3/OS_SPATA18.pdf", height = 6, width = 6)
tcga_kmplot("LIHC", "SPATA18")
dev.off()
pdf(file = "path/Part3/OS_TEC.pdf", height = 6, width = 6)
tcga_kmplot("LIHC", "TEC")
dev.off()
pdf(file = "path/Part3/OS_TNNT3.pdf", height = 6, width = 6)
tcga_kmplot("LIHC", "TNNT3")
dev.off()
pdf(file = "path/Part3/OS_WDR97.pdf", height = 6, width = 6)
tcga_kmplot("LIHC", "WDR97")
dev.off()
pdf(file = "path/Part3/OS_ZNF285.pdf", height = 6, width = 6)
tcga_kmplot("LIHC", "ZNF285")
dev.off()

# Integrate meth and RNA

## The expression of gene involveed with meth

RNA <- read.csv("/home/yangwl/data/HCC_bulk/result/RNA_last/HCC_Pos_Neg_expression_TPM_gene.csv", header=TRUE, row.names=1)
NegT_PosT <- RNA[ ,c(2,4,6,8,10,12,14,16,18,20,22,23,25,27,29,31,33,35,37,39)]
meth_gene <- c("DNMT1", "DNMT3A", "DNMT3B", "DNMT3L", "TET1", "TET2", "TET3")
RNA_meth <- t(NegT_PosT[meth_gene,])
res <- data.frame(RNA_meth)%>%
  mutate(group = c(rep('HCC_Neg_T',10),rep('HCC_Pos_T',10)))%>%
  rownames_to_column("sample")%>%
  pivot_longer(cols = colnames(.)[2:8],
               names_to = "genes",
               values_to = 'value')
res$group <- factor(res$group,levels =c("HCC_Pos_T","HCC_Neg_T"))

pdf(file = "path/Part4/meth_genes_expression.pdf", height = 4, width = 6)
ggplot(res,aes(genes,value,fill = group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "", y = "TPM") +
  theme(legend.position = "top") + 
  theme(
        panel.grid=element_blank(),axis.text.x = element_text(angle=90,vjust = 0.5,size = 14,colour = 'black'),
        axis.text.y = element_text(size = 14,colour = 'black'))+
  scale_fill_nejm()+
  stat_compare_means(aes(group = group),label = "p.format",size=3,method = "kruskal.test")
dev.off()

## The expression of DMR or DhMR associated genes (RNA)

meth_associated_gene <- c("ARL6IP4","CCDC200","DHDH","DUX4L50","ESPNP","FCGBP","KCNQ3","KCNS1","LRATD2",
                          "MEX3B","MIR1257","MIR9-1HG","MTND3P22","OGFOD3","POU4F1","PSMA8","PSME2P6",
                          "RARA","RNA5-8SN1","RNA5-8SN2","RNA5-8SN3","SEL1L3","SLC22A23","SPATA18",
                          "TEC","TNNT3","WDR97","ZNF285")

RNA_meth_associated <- t(NegT_PosT[meth_associated_gene,])
res <- data.frame(RNA_meth_associated)%>%
  mutate(group = c(rep('HCC_Neg_T',10),rep('HCC_Pos_T',10)))%>%
  rownames_to_column("sample")%>%
  pivot_longer(cols = colnames(.)[2:29],
               names_to = "genes",
               values_to = 'value')
res$group <- factor(res$group,levels =c("HCC_Pos_T","HCC_Neg_T"))

pdf(file = "path/Part4/meth_associated_genes_expression.pdf", height = 5, width = 17)
ggplot(res,aes(genes,value,fill = group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "", y = "TPM") +
  theme(legend.position = "top") + 
  theme(
    panel.grid=element_blank(),axis.text.x = element_text(angle=90,vjust = 0.5,size = 14,colour = 'black'),
    axis.text.y = element_text(size = 14,colour = 'black'))+
  scale_fill_nejm()+
  stat_compare_means(aes(group = group),label = "p.format",size=3,method = "kruskal.test")
dev.off()