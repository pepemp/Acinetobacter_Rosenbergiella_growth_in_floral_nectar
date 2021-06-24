library(readxl) 
library("beanplot")
library(RColorBrewer)

## DATA LOADING
# load optical density z-scores (growth scores) table
Zscores <- read_excel("./data/Abt_ANs_data.xlsx", sheet = "Z-scores")
Zscores <- as.data.frame(Zscores) 
# set isolate column as index
rownames(Zscores) <- Zscores$Isolate 
Zscores$Isolate <- NULL 
# the dataframe gets 24 columns (12 artificials nectars  x 2 reads (3 y 7 días))
nb_col = ncol(Zscores) 
View(Zscores)



## PART 4 - Data distribution

### 4.1. Does the data fit a normal distribution? 
  
# Shapiro-Wilk normality test
shapiro.p.values = matrix(0L, nrow = 1, ncol = nb_col)
for(i in 1:nb_col) {
  shapiro.p.values[i] = shapiro.test(Zscores[ , i])$p.value
}
# Adjust p-values according to Holm's method
shapiro.p.values.corregidos = p.adjust(shapiro.p.values, method = "holm")
shapiro.p.values.corregidos = matrix(
  shapiro.p.values.corregidos, 
  nrow = 1, 
  ncol = nb_col, 
  dimnames = list(NULL,colnames(Zscores))
)
print("p-values of Shapiro test corrected by Holm's method:") 
print(shapiro.p.values.corregidos) 
# 14 variables follow a normal distribution
  
### 4.2. Generation of bean plots, violin plots or box-and-whiskers plots

colorBlindJacksonLab <- list("#009292", "#e69f00", "#6db6ff", "#b4eb34")
# loading metadata to extract Species and Origin
Metadata <- read_excel("./data/Abt_ANs_data.xlsx", sheet = "Abt_metadata")
Metadata <- as.data.frame(Metadata[,1:3]) 
# set isolate column as index
rownames(Metadata) <- Metadata$Isolate
# in the aim to keep same row order in Z-scores and Metadata, we order both dataframe by row names:
Metadata <- Metadata[ order(row.names(Metadata)), ] 
Zscores <- Zscores[ order(row.names(Zscores)), ] 


# beanplots according to species

par(mfrow=c(3, 2)) # culture media 1 to 6 at 3 days
bp_AN01_3d<-Zscores$AN01_3d~Metadata$Species 
beanplot(bp_AN01_3d,col=colorBlindJacksonLab[1:2], bw="nrd0", xlab="Species", ylab="AN01 (3 days)")
bp_AN02_3d<-Zscores$AN02_3d~Metadata$Species
beanplot(bp_AN02_3d,col=colorBlindJacksonLab[1:2], bw="nrd0", xlab="Species", ylab="AN02 (3 days)")
bp_AN03_3d<-Zscores$AN03_3d~Metadata$Species
beanplot(bp_AN03_3d,col=colorBlindJacksonLab[1:2], bw="nrd0", xlab="Species", ylab="AN03 (3 days)")
bp_AN04_3d<-Zscores$AN04_3d~Metadata$Species
beanplot(bp_AN04_3d,col=colorBlindJacksonLab[1:2], bw="nrd0", xlab="Species", ylab="AN04 (3 days)")
bp_AN05_3d<-Zscores$AN05_3d~Metadata$Species
beanplot(bp_AN05_3d,col=colorBlindJacksonLab[1:2], bw="nrd0", xlab="Species", ylab="AN05 (3 days)")
bp_AN06_3d<-Zscores$AN06_3d~Metadata$Species
beanplot(bp_AN06_3d,col=colorBlindJacksonLab[1:2], bw="nrd0", xlab="Species", ylab="AN06 (3 days)")

par(mfrow=c(3, 2)) # culture media 7 to 12 at 3 days
bp_AN07_3d<-Zscores$AN07_3d~Metadata$Species
beanplot(bp_AN07_3d,col=colorBlindJacksonLab[1:2], bw="nrd0", xlab="Species", ylab="AN07 (3 days)")
bp_AN08_3d<-Zscores$AN08_3d~Metadata$Species
beanplot(bp_AN08_3d,col=colorBlindJacksonLab[1:2], bw="nrd0", xlab="Species", ylab="AN08 (3 days)")
bp_AN09_3d<-Zscores$AN09_3d~Metadata$Species
beanplot(bp_AN09_3d,col=colorBlindJacksonLab[1:2], bw="nrd0", xlab="Species", ylab="AN09 (3 days)")
bp_AN10_3d<-Zscores$AN10_3d~Metadata$Species
beanplot(bp_AN10_3d,col=colorBlindJacksonLab[1:2], bw="nrd0", xlab="Species", ylab="AN10 (3 days)")
bp_AN11_3d<-Zscores$AN11_3d~Metadata$Species
beanplot(bp_AN11_3d,col=colorBlindJacksonLab[1:2], bw="nrd0", xlab="Species", ylab="AN11 (3 days)")
bp_AN12_3d<-Zscores$AN12_3d~Metadata$Species
beanplot(bp_AN12_3d,col=colorBlindJacksonLab[1:2], bw="nrd0", xlab="Species", ylab="AN12 (3 days)")

par(mfrow=c(3, 2)) # culture media 1 to 6 at 7 days
bp_AN01_7d<-Zscores$AN01_7d~Metadata$Species 
beanplot(bp_AN01_7d,col=colorBlindJacksonLab[1:2], bw="nrd0", xlab="Species", ylab="AN01 (7 days)")
bp_AN02_7d<-Zscores$AN02_7d~Metadata$Species
beanplot(bp_AN02_7d,col=colorBlindJacksonLab[1:2], bw="nrd0", xlab="Species", ylab="AN02 (7 days)")
bp_AN03_7d<-Zscores$AN03_7d~Metadata$Species
beanplot(bp_AN03_7d,col=colorBlindJacksonLab[1:2], bw="nrd0", xlab="Species", ylab="AN03 (7 days)")
bp_AN04_7d<-Zscores$AN04_7d~Metadata$Species
beanplot(bp_AN04_7d,col=colorBlindJacksonLab[1:2], bw="nrd0", xlab="Species", ylab="AN04 (7 days)")
bp_AN05_7d<-Zscores$AN05_7d~Metadata$Species
beanplot(bp_AN05_7d,col=colorBlindJacksonLab[1:2], bw="nrd0", xlab="Species", ylab="AN05 (7 days)")
bp_AN06_7d<-Zscores$AN06_7d~Metadata$Species
beanplot(bp_AN06_7d,col=colorBlindJacksonLab[1:2], bw="nrd0", xlab="Species", ylab="AN06 (7 days)")

par(mfrow=c(3, 2)) # culture media 7 to 12 at 7 days
bp_AN07_7d<-Zscores$AN07_7d~Metadata$Species
beanplot(bp_AN07_7d,col=colorBlindJacksonLab[1:2], bw="nrd0", xlab="Species", ylab="AN07 (7 days)")
bp_AN08_7d<-Zscores$AN08_7d~Metadata$Species
beanplot(bp_AN08_7d,col=colorBlindJacksonLab[1:2], bw="nrd0", xlab="Species", ylab="AN08 (7 days)")
bp_AN09_7d<-Zscores$AN09_7d~Metadata$Species
beanplot(bp_AN09_7d,col=colorBlindJacksonLab[1:2], bw="nrd0", xlab="Species", ylab="AN09 (7 days)")
bp_AN10_7d<-Zscores$AN10_7d~Metadata$Species
beanplot(bp_AN10_7d,col=colorBlindJacksonLab[1:2], bw="nrd0", xlab="Species", ylab="AN10 (7 days)")
bp_AN11_7d<-Zscores$AN11_7d~Metadata$Species
beanplot(bp_AN11_7d,col=colorBlindJacksonLab[1:2], bw="nrd0", xlab="Species", ylab="AN11 (7 days)")
bp_AN12_7d<-Zscores$AN12_7d~Metadata$Species
beanplot(bp_AN12_7d,col=colorBlindJacksonLab[1:2], bw="nrd0", xlab="Species", ylab="AN12 (7 days)")


## beanplots según el hábitat

par(mfrow=c(3, 2)) # culture media 1 to 6 at 3 days
bp_AN01_3d<-Zscores$AN01_3d~Metadata$Origin
beanplot(bp_AN01_3d,col=colorBlindJacksonLab[3:4], bw="nrd0", xlab="Habitat", ylab="AN01 (3 days)")
bp_AN02_3d<-Zscores$AN02_3d~Metadata$Origin
beanplot(bp_AN02_3d,col=colorBlindJacksonLab[3:4], bw="nrd0", xlab="Habitat", ylab="AN02 (3 days)")
bp_AN03_3d<-Zscores$AN03_3d~Metadata$Origin
beanplot(bp_AN03_3d,col=colorBlindJacksonLab[3:4], bw="nrd0", xlab="Habitat", ylab="AN03 (3 days)")
bp_AN04_3d<-Zscores$AN04_3d~Metadata$Origin
beanplot(bp_AN04_3d,col=colorBlindJacksonLab[3:4], bw="nrd0", xlab="Habitat", ylab="AN04 (3 days)")
bp_AN05_3d<-Zscores$AN05_3d~Metadata$Origin
beanplot(bp_AN05_3d,col=colorBlindJacksonLab[3:4], bw="nrd0", xlab="Habitat", ylab="AN05 (3 days)")
bp_AN06_3d<-Zscores$AN06_3d~Metadata$Origin
beanplot(bp_AN06_3d,col=colorBlindJacksonLab[3:4], bw="nrd0", xlab="Habitat", ylab="AN06 (3 days)")

par(mfrow=c(3, 2)) # culture media 7 to 12 at 3 days
bp_AN07_3d<-Zscores$AN07_3d~Metadata$Origin 
beanplot(bp_AN07_3d,col=colorBlindJacksonLab[3:4], bw="nrd0", xlab="Habitat", ylab="AN07 (3 days)")
bp_AN08_3d<-Zscores$AN08_3d~Metadata$Origin
beanplot(bp_AN08_3d,col=colorBlindJacksonLab[3:4], bw="nrd0", xlab="Habitat", ylab="AN08 (3 days)")
bp_AN09_3d<-Zscores$AN09_3d~Metadata$Origin
beanplot(bp_AN09_3d,col=colorBlindJacksonLab[3:4], bw="nrd0", xlab="Habitat", ylab="AN09 (3 days)")
bp_AN10_3d<-Zscores$AN10_3d~Metadata$Origin
beanplot(bp_AN10_3d,col=colorBlindJacksonLab[3:4], bw="nrd0", xlab="Habitat", ylab="AN10 (3 days)")
bp_AN11_3d<-Zscores$AN11_3d~Metadata$Origin
beanplot(bp_AN11_3d,col=colorBlindJacksonLab[3:4], bw="nrd0", xlab="Habitat", ylab="AN11 (3 days)")
bp_AN12_3d<-Zscores$AN12_3d~Metadata$Origin
beanplot(bp_AN12_3d,col=colorBlindJacksonLab[3:4], bw="nrd0", xlab="Habitat", ylab="AN12 (3 days)")

par(mfrow=c(3, 2)) # culture media 1 to 6 at 7 days
bp_AN01_7d<-Zscores$AN01_7d~Metadata$Origin 
beanplot(bp_AN01_7d,col=colorBlindJacksonLab[3:4], bw="nrd0", xlab="Habitat", ylab="AN01 (7 days)")
bp_AN02_7d<-Zscores$AN02_7d~Metadata$Origin
beanplot(bp_AN02_7d,col=colorBlindJacksonLab[3:4], bw="nrd0", xlab="Habitat", ylab="AN02 (7 days)")
bp_AN03_7d<-Zscores$AN03_7d~Metadata$Origin
beanplot(bp_AN03_7d,col=colorBlindJacksonLab[3:4], bw="nrd0", xlab="Habitat", ylab="AN03 (7 days)")
bp_AN04_7d<-Zscores$AN04_7d~Metadata$Origin
beanplot(bp_AN04_7d,col=colorBlindJacksonLab[3:4], bw="nrd0", xlab="Habitat", ylab="AN04 (7 days)")
bp_AN05_7d<-Zscores$AN05_7d~Metadata$Origin
beanplot(bp_AN05_7d,col=colorBlindJacksonLab[3:4], bw="nrd0", xlab="Habitat", ylab="AN05 (7 days)")
bp_AN06_7d<-Zscores$AN06_7d~Metadata$Origin
beanplot(bp_AN06_7d,col=colorBlindJacksonLab[3:4], bw="nrd0", xlab="Habitat", ylab="AN06 (7 days)")

par(mfrow=c(3, 2)) # culture media 7 to 12 at 7 days
bp_AN07_7d<-Zscores$AN07_7d~Metadata$Origin
beanplot(bp_AN07_7d,col=colorBlindJacksonLab[3:4], bw="nrd0", xlab="Habitat", ylab="AN07 (7 days)")
bp_AN08_7d<-Zscores$AN08_7d~Metadata$Origin
beanplot(bp_AN08_7d,col=colorBlindJacksonLab[3:4], bw="nrd0", xlab="Habitat", ylab="AN08 (7 days)")
bp_AN09_7d<-Zscores$AN09_7d~Metadata$Origin
beanplot(bp_AN09_7d,col=colorBlindJacksonLab[3:4], bw="nrd0", xlab="Habitat", ylab="AN09 (7 days)")
bp_AN10_7d<-Zscores$AN10_7d~Metadata$Origin
beanplot(bp_AN10_7d,col=colorBlindJacksonLab[3:4], bw="nrd0", xlab="Habitat", ylab="AN10 (7 days)")
bp_AN11_7d<-Zscores$AN11_7d~Metadata$Origin
beanplot(bp_AN11_7d,col=colorBlindJacksonLab[3:4], bw="nrd0", xlab="Habitat", ylab="AN11 (7 days)")
bp_AN12_7d<-Zscores$AN12_7d~Metadata$Origin
beanplot(bp_AN12_7d,col=colorBlindJacksonLab[3:4], bw="nrd0", xlab="Habitat", ylab="AN12 (7 days)")


### 4.3. Analysis of the pairwise correlations between traits and generation of correlograms

library("Hmisc")
library(WriteXLS)

# Spearman correlation
CM <- rcorr(as.matrix(Zscores), type="spearman")
r <- as.data.frame(CM$r) # Spearman correlation matrix
P <- as.data.frame(CM$P) # p-values associated to the correlation matrix
P[upper.tri(P)] <- NA

# Adjust p-values according to Holm's method
Padj <- p.adjust(as.matrix(P), method="holm")
Padj <- matrix(Padj,ncol=24,nrow=24)
Padj <- as.data.frame(Padj)

# Save in an excel file
library(xlsx)
write.xlsx(r, file="./results/Abt_zscores_correlation.xlsx", sheetName="Spearman corr coef", row.names=TRUE)
write.xlsx(P, file="./results/Abt_zscores_correlation.xlsx", sheetName="p-values", append=TRUE, row.names=TRUE)
write.xlsx(Padj, file="./results/Abt_zscores_correlation.xlsx", sheetName="adjusted p-values", append=TRUE, row.names=TRUE)


# Correlograms
library(corrplot)
par(mfrow=c(1, 1))
corrplot(
  CM$r, method= "color", type = "full", order = "original", 
  p.mat = as.matrix(Padj), sig.level = 0.05, insig="blank", 
  addgrid.col = "grey", col = brewer.pal("Spectral", n=7), tl.col = "black"
) # the lower diagonal shows significant tests


## PART 5 - Principal Component Analysis (PCA) of the growth data.

# Create PCA
pca <- princomp(Zscores, cor = FALSE) # "cor = FALSE" => covariance matrix
summary(pca)
# The 2 first principal components explain 55.09% of the variance. 

### 5.1. Scree plots
variance <- (pca / sum(pca)) * 100
v <- round(((pca$sdev^2 / sum(pca$sdev^2)) * 100), 2)
x <- barplot(v, ylim=c(0,60), ylab="% Variance", xlab="Component", col="limegreen", xaxt="n")
text(cex=1, x=x, y=-4.5, names(v), xpd=TRUE, srt=45)

# common theme for plots
theme.plot = theme(
  axis.text=element_text(size=12), 
  axis.title=element_text(size=12), 
  plot.title = element_text(size=12, face="bold"), 
  axis.title.x = element_text(size=12, face="bold"), 
  axis.title.y = element_text(size=12, face="bold"), 
  panel.border = element_rect(fill=NA, colour = "black", linetype="solid"), 
  panel.background = element_rect(fill = "white", colour = "black", linetype="solid"),
  axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
)

### 5.2. Loading plots
library(ggrepel) 
loadings<-cbind(pca$loadings[,"Comp.1"], pca$loadings[,"Comp.2"])
colnames(loadings)<-c("Component1", "Component2")
loadings<-as.data.frame(loadings)
loadings$label <- rownames(loadings)
loadings

ggplot(loadings, aes(x=Component1, y=Component2, label = label)) + 
  geom_point(size=2, color="#920000") + labs(title="", x="Component 1", y = "Component 2") +
  theme.plot + 
  geom_text_repel(point.padding = 0.1, segment.color = 'grey50', color = 'black')

ggplot(data = loadings, aes(x = Component1, y = Component2, label = label)) + 
  coord_equal() + 
  geom_segment(
    data = loadings, aes(x = 0, y = 0, xend = Component1, yend = Component2), 
    arrow = arrow(length = unit(1/2, 'picas')), color = 'firebrick') + 
  labs(title="", x="Component 1", y = "Component 2") + 
  theme.plot + 
  geom_text(
    data = loadings, 
    aes(
      label = label, 
      x = Component1, 
      y = Component2, 
      angle = with(loadings, (180/pi) * atan(Component2 /  Component1)), 
      hjust = with(loadings, (1 - 1.5 * sign(Component1)) / 2)), 
    color = 'firebrick') + 
  xlim(c(-0.7,0.6)) + ylim(c(-0.5,0.6)) #Loadings plot with arrows


# Bar plots for the loadings
c.pc1 <- ifelse(loadings[,1] > 0, yes="#009E73", no="#E69F00")
c.pc2 <- ifelse(loadings[,2] > 0, yes="#009E73", no="#E69F00")

bp1 <- ggplot(loadings, aes(x = label, y=Component1)) + 
  geom_bar(stat="identity", fill=c.pc1, colour = "black") + 
  labs(title="Conventional PCA", x="", y = "Loadings PC1") + 
  theme.plot 

bp2 <- ggplot(
  loadings, aes(x = label, y=Component2)) + 
  geom_bar(stat="identity", fill=c.pc2, colour = "black") + 
  labs(title="Conventional PCA", x="", y = "Loadings PC2") + 
  theme.plot 

library(gridExtra)
grid.arrange(bp1, bp2, nrow = 2) 


### 5.3. Scores plots
scores<-cbind(pca$scores[,"Comp.1"], pca$scores[,"Comp.2"]) # to retain the loadings for the first two PCs
colnames(scores)<-c("Component1", "Component2") # to rename the columns
scores<-as.data.frame(scores) # to transform the matrix into a data.frame
scores$label<-rownames(scores) # to create a new column for the labels of points in the plots
scores$species<-Metadata$Species # to create a new column for the species
scores$habitat<-Metadata$Origin # to create a new column for the habitat
head(scores)
tail(scores)

ggplot(scores, aes(x=Component1, y=Component2, label = label, color = species)) + 
  geom_point(size=2) + labs(title="", x="Component 1", y = "Component 2") + 
  theme.plot +
  scale_color_manual(name="",values=unlist(colorBlindJacksonLab[1:2])) +
  stat_ellipse(type = "norm", level = 0.68) # By species
# the ellipses corresponding to the 2 species are cleary separated


ggplot(scores, aes(x=Component1, y=Component2, label = label, color = habitat)) + 
  geom_point(size=2) + labs(title="", x="Component 1", y = "Component 2") + 
  theme.plot +
  scale_color_manual(name="",values=unlist(colorBlindJacksonLab[3:4])) + 
  stat_ellipse(type = "norm", level = 0.68) # By habitat
# the ellipses corresponding to the 2 habitats are overlapped