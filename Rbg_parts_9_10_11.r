library(readxl) 

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
library(stringr)
accession.numbers <- read_excel(".data/Rbg_ANs_data.xlsx", sheet = "Abt_metadata")
accession.numbers <- as.data.frame(accession.numbers) 
rownames(accession.numbers) <- accession.numbers$rpoB # set it as index
tree3$tip.label = str_split_fixed(tree3$tip.label, ".1Acin", 2)[, 1]
for (i in 1:length(tree3$tip.label)) { # replace tree leaves names by the corresponding isolate name
  tree3$tip.label[i] <- accession.numbers[tree3$tip.label[i],]$Isolate
}



## PART 9 - PHYLOGENETIC SIGNAL METRICS

library(picante)

# For calculation of K
x<-multiPhylosignal(traits,tree3,reps=1000)
x$Padj<-p.adjust(x$PIC.variance.P,method="holm") # Add adjusted p-values
x ## 12 tests with significant corrected p-values

# For calculation of lambda
# empty matrix to store lambda results:
lambda<-matrix(NA,ncol(traits),5,dimnames=list(NULL,c("test","lambda","logL","logL0","P")))
# Fill lambda dataframe
for (i in 1:24){
  a<-phylosig(tree3,traits[,i],method="lambda",test=T,nsim=1000)
  lambda[i,]<-unlist(c(colnames(traits)[i],a$lambda,a$logL,a$logL0,a$P)) #Esta lÃ???nea es para rellenar la matriz de resultados creada anteriormente
}
# Add corrected p-values
lambda<-as.data.frame(lambda) #Convertimos el objeto 'lambda' en data.frame
lambda$Padj<-p.adjust(as.numeric(as.character(lambda$P)),method="holm") #Para aÃ±adir la columna de valores P ajustados
lambda ## ¡Muy bien! A mí me salen 14 test con P corregida significativa. Y nota que tienes bastantes valores de lambda próximos a 1 (15 son mayores de 0.8).



## PART 10 - PHYLOGENETIC HEATMAP

library(viridis) 
phylo.heatmap(tree3, traits, fsize=c(1,0.75,1),colors=plasma(200)) 


## PART 11 - CONTMAP
par(mfrow=c(2,3)) #Super figure of 3 x 2 subfigures

# contmap for isolates grown in culture media from 1 to 6 during 3 days
for (i in 1:6){ 
  a<-traits[,i] # i-th column of phenotypic table 
  names(a)<- row.names(traits)
  obj<-contMap(tree3,a,fsize=c(0.8,1),lwd=2,plot=TRUE) # Application of function ConMap to columna 'a'
  obj<-setMap(obj,colors=viridis(100)) 
  plot(obj,direction="rightwards",leg.txt="Growth score",sig=2,legend=0.05,mar=c(0.1,0.1,1.1,0.1))
  title(main=colnames(as.data.frame(traits))[i],font.main=3) # Set title
}

# contmap for isolates grown in culture media from 7 to 12 during 3 days
for (i in 7:12){
  a<-traits[,i] 
  names(a)<- row.names(traits)
  obj<-contMap(tree3,a,fsize=c(0.3,0.5),lwd=2,plot=FALSE) 
  obj<-setMap(obj,colors=viridis(100)) 
  plot(obj,direction="rightwards",leg.txt="Growth score",sig=2,legend=0.05,mar=c(0.1,0.1,1.1,0.1))
  title(main=colnames(as.data.frame(traits))[i],font.main=3) 
}

# contmap for isolates grown in culture media from 1 to 6 during 7 days
for (i in 13:18){
  a<-traits[,i] 
  names(a)<- row.names(traits)
  obj<-contMap(tree3,a,fsize=c(0.3,0.5),lwd=2,plot=FALSE) 
  obj<-setMap(obj,colors=viridis(100)) 
  plot(obj,direction="rightwards",leg.txt="Growth score",sig=2,legend=0.05,mar=c(0.1,0.1,1.1,0.1))
  title(main=colnames(as.data.frame(traits))[i],font.main=3) 
}

# contmap for isolates grown in culture media from 7 to 12 during 7 days
for (i in 19:24){
  a<-traits[,i] 
  names(a)<- row.names(traits)
  obj<-contMap(tree3,a,fsize=c(0.3,0.5),lwd=2,plot=FALSE) 
  obj<-setMap(obj,colors=viridis(100)) 
  plot(obj,direction="rightwards",leg.txt="Growth score",sig=2,legend=0.05,mar=c(0.1,0.1,1.1,0.1))
  title(main=colnames(as.data.frame(traits))[i],font.main=3)
}
# In some graphics, the phylogenetic signal of the analyzed character stands out.
