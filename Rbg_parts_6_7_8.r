## Import phylogenetic tree
library(phytools)

tree<-read.tree("./data/Rbg_NJ_tree_concatamer.nwk")
outgroup.name<-"JN808190_JF745803_JF745804_Phaseolibacter_flectens"
outgroup<-match(outgroup.name,tree$tip.label) # establish the outgroup
rtree<-root(tree,outgroup,resolve.root=TRUE) # rooting the tree
tree2<-drop.tip(rtree, outgroup.name)

is.ultrametric(tree2) # Check if the tree is ultrametric
tree3<-force.ultrametric(tree2)
is.ultrametric(tree3)
# Add a slight quantity in branch length to avoid zeros:
tree3$edge.length<-tree3$edge.length + 10^-10 
tree3$edge.length[tree3$edge.length==0]

# Replace tree leaves names by isolate names
library(readxl) 
setwd("C:/Users/moral/Documents/UPM/Stage")
accession.numbers <- read_excel("Rbg_Genbank_accession_numbers.xlsx")
accession.numbers <- as.data.frame(accession.numbers) #Para transformar en data.frame
accession.numbers$Species <- lapply(accession.numbers$Species, gsub, pattern = " ", replacement = "_", fixed = TRUE) # reconstruct tree leaves names
rownames(accession.numbers) <- paste(accession.numbers$atpD, accession.numbers$gyrB, accession.numbers$rpoB, accession.numbers$Species, sep="_") #set it as index
for (i in 1:length(tree3$tip.label)) { # replace tree leaves names by the corresponding isolate name
  tree3$tip.label[i] <- accession.numbers[tree3$tip.label[i],]$Isolate
}

## Import phenotypic data
traits <- read_excel("./data/Rbg_ANs_data.xlsx", sheet = "Z-scores")
traits <- as.data.frame(traits)
# set isolate as index:
rownames(traits) <- traits$Isolate
traits$Isolate <- NULL
traits <- traits[ order(row.names(traits)), ] # sort row by names


## Import metadata (origin of isolates)
# Loading metadata to extract Species and Origin
Metadata <- read_excel("./data/Rbg_ANs_data.xlsx", sheet = "Abt_metadata")
Metadata <- as.data.frame(Metadata[,1:3]) 
# set isolate as index:
rownames(Metadata) <- Metadata$Isolate 
# in the aim to keep same row order in Zscores and Metadata, we order both dataframe by row names
Metadata <- Metadata[ order(row.names(Metadata)), ] # sort row by names
groups<-Metadata[,"Origin"] # Create an object "groups" to store origins (insect vs nectar)
groups<-as.factor(groups)
names(groups)<-row.names(Metadata)


## PART 6 - Analysis of the pairwise correlations between traits using PICs and generation of new correlograms

### 6.1. Compute phylogenetic independent contrasts (PICs)
PICs<-list() # Generate empty list to store PICs
for (i in 1:ncol(traits)) { 
     t<-traits[,i]
     PICs[[i]]<-pic(t,tree3)
}

df <- data.frame(matrix(unlist(PICs), ncol=ncol(traits), byrow=F))

### 6.2. CORRELATIONS: Analysis of the pairwise correlations between traits and generation of correlograms
library(RColorBrewer)
library("Hmisc")
CM <- rcorr(as.matrix(df_PICs), type="spearman")
r <- as.data.frame(CM$r) # Spearman correlation matrix
P <- as.data.frame(CM$P) # p-values associated to the correlation matrix
P[upper.tri(P)] <- NA

# Adjust p-values according to Holm's method
Padj <- p.adjust(as.matrix(P), method="holm")
Padj <- matrix(Padj,ncol=24,nrow=24)
Padj <- as.data.frame(Padj)

# CORRELOGRAMS
library(corrplot)
par(mfrow=c(1, 1))
corrplot(
  CM$r, method= "color", type = "full", order = "original", 
  p.mat = as.matrix(Padj), sig.level = 0.05, insig="blank", 
  addgrid.col = "grey", col = brewer.pal("Spectral", n=7), tl.col = "black"
) ## We can see that the correlation between same conditions at 3d and 7d are higher than the average

# Save in an excel file
library(xlsx)
write.xlsx(r, file="./results/Rbg_PICs_correlation.xlsx", sheetName="Spearman corr coef", row.names=TRUE)
write.xlsx(P, file="./results/Rbg_PICs_correlation.xlsx", sheetName="p-values", append=TRUE, row.names=TRUE)
write.xlsx(Padj, file="./results/Rbg_PICs_correlation.xlsx", sheetName="adjusted p-values", append=TRUE, row.names=TRUE)



## PART 7 - PHYLOGENETIC ANOVA TESTS
aov<-list() # Generate an empty list to store ANOVA phylogenetic results
for (i in 1:ncol(traits)) {
  a<-as.vector(traits[,i])
  names(a)<-row.names(traits)
  aov[[i]] <-phylANOVA(tree3,groups,a)
}
aov

# store results in a table:
df_aov <- matrix(NA,ncol(traits),3,dimnames=list(NULL,c("test","F value","Pr(>F)"))) 
for (i in 1:ncol(traits)) {
  df_aov[i,]<-c(colnames(traits)[i],aov[[i]]$F,aov[[i]]$Pf)
}
df_aov<-as.data.frame(df_aov) 
# Holm's correction of p-values:
pval<-p.adjust(as.numeric(as.character(df_aov[,3])),method="holm") 
pval

df_aov$adjP<-pval # Add corrected p-values to the table of results
df_aov
# Here we can see that all corrected p-values are equal to one. 
# There is hence no significant difference between groups of isolates 
# (e.g., insect vs. nectar isolates) by phylogenetic ANOVA.



## PART 8 - PHYLOGENETIC PRINCIPAL COMPONENT ANALYSIS OF THE GROWTH DATA 

phylo_pca<-phyl.pca(tree3, traits, method="BM", mode="cov") # Phylogenetic PCA  
summary(phylo_pca) # detailed results

### 8.1. Scree plot 
variance<-((phylo_pca$Eval/sum(phylo_pca$Eval))*100)
v<-round(diag(variance),2)
for (i in 1:length(names(v))) {
  names(v)[i] <- paste0("Comp.", i)
}
x<-barplot(v, ylim=c(0,60), ylab="% Variance", xlab="Component", col="limegreen", xaxt="n")
text(cex=1, x=x, y=-4.5, names(v), xpd=TRUE, srt=45)
# The 2 first principal components explain 45.08% of the variance. 


### 8.2. Loadings plot
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

library(ggrepel) 
loadings<-cbind(phylo_pca$L[,"PC1"], phylo_pca$L[,"PC2"])
colnames(loadings)<-c("Component1", "Component2")
loadings<-as.data.frame(loadings)
loadings$label <- rownames(loadings)
loadings

ggplot(loadings, aes(x=Component1, y=Component2, label = label)) + 
  geom_point(size=2, color="#0072B2") + labs(title="", x="Component 1", y = "Component 2") + 
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
  xlim(c(-1.6, 1)) + ylim(c(-1.4, 1.3)) #Loadings plot with arrows


# Bar plots for the loadings
c.pc1 <- ifelse(loadings[,1] > 0, yes="#009E73", no="#E69F00")
c.pc2 <- ifelse(loadings[,2] > 0, yes="#009E73", no="#E69F00")

bp1 <- ggplot(loadings, aes(x = label, y=Component1)) + 
  geom_bar(stat="identity", fill=c.pc1, colour = "black") + 
  labs(title="Phylogenetic PCA", x="", y = "Loadings PC1") +
  theme.plot 

bp2 <- ggplot(
  loadings, aes(x = label, y=Component2)) + 
  geom_bar(stat="identity", fill=c.pc2, colour = "black") + 
  labs(title="Phylogenetic PCA", x="", y = "Loadings PC2") +
  theme.plot 

library(gridExtra)
grid.arrange(bp1, bp2, nrow = 2) 


### 8.3. PCA scores plots by species 
scores<-cbind(phylo_pca$S[,"PC1"], phylo_pca$S[,"PC2"]) # to retain the loadings for the first two PCs
colnames(scores)<-c("Component1", "Component2") # to rename the columns
scores<-as.data.frame(scores) # to transform the matrix into a data.frame
scores$label<-rownames(scores) # to create a new column for the labels of points in the plots
scores$species<-Metadata$Species # to create a new column for the species
scores$habitat<-Metadata$Origin # to create a new column for the habitat
head(scores)
tail(scores)

colorBlindGrey8 <- list("#999999", "#E69F00", "#56B4E9", "#009E73", 
                        "#F0E442", "#0072B2", "#CC79A7", "#D55E00")

ggplot(scores, aes(x=Component1, y=Component2, label = label, color = species)) + 
  geom_point(size=2) + labs(title="", x="Component 1", y = "Component 2") + 
  theme.plot +
  scale_color_manual(name="",values=unlist(colorBlindGrey8[1:5])) +
  stat_ellipse(type = "norm", level = 0.68) # By species
# ellipses are overlapped, there is no significant differences between species


ggplot(scores, aes(x=Component1, y=Component2, label = label, color = habitat)) + 
  geom_point(size=2) + labs(title="", x="Component 1", y = "Component 2") + 
  theme.plot +
  scale_color_manual(name="",values=unlist(colorBlindGrey8[6:7])) + 
  stat_ellipse(type = "norm", level = 0.68) # By habitat
# ellipses are overlapped, there is no significant differences between habitats