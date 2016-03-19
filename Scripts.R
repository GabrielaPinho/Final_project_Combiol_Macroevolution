#################                  Arvore gerada com BEAST
library(geiger)
library(phytools)
library(ape)
library(ggtree)
library(phylolm)
library(nlme)

setwd ("/home/gabriela/Dropbox/Documents/Doutorado/Courses_conferences/ComBiol_Macroevolution/Lab/Projeto_Mike/Supplementary_material")

Tree_3genes <- read.nexus("3genes_final.tree")
Tree_cytb <- read.nexus("Cytb_final.tree")
Tree_3genes <- ladderize(Tree_3genes)
Tree_cytb <- ladderize(Tree_cytb)

par(mfrow= c(1,2))
plot(Tree_3genes)
plot(Tree_cytb)

#topologies are not equal, so I am using the tree with more information (3 genes)

#####                                                  Test body size evolution model
tree <- Tree_3genes

Ge2014 <- read.csv("Ge_2014_initial.csv")
head(Ge2014)
length(Ge2014$Species)
rownames(Ge2014) <- Ge2014$Species #changing the rownames() in the traits data.frame to the tip labels
Ge2014 <- Ge2014[match(tree$tip.label, rownames(Ge2014)),]
cbind(Ge2014,tree$tip.label)

Body_length <- Ge2014$Body_length
names(Body_length) <- tree$tip.label

Body_length.bm <- fitContinuous(tree, Body_length, model = "BM")
Body_length.ou <- fitContinuous(tree, Body_length, model = "OU")
Body_length.eb <- fitContinuous(tree, Body_length, model = "EB")
Body_length.wn <- fitContinuous(tree, Body_length, model = "white")
Body_lengthAIC <- c(Body_length.bm$opt$aic , Body_length.ou$opt$aic, Body_length.eb$opt$aic, Body_length.wn$opt$aic)

names(Body_lengthAIC) <- c("Brownian Motion", "Ornstein–Uhlenbeck", "Early Burst", "White noise")
Body_lengthAIC # you want the model with the lowest AIC

#######################################################################################################
#### BAMM body length

#install.packages("BAMMtools")
library(BAMMtools)
library(coda)

Tree_3genes <- read.nexus("3genes_final.tree")
Tree_3genes <- ladderize(Tree_3genes)
tree <- Tree_3genes
length(tree$tip.label)
plot(tree)

is.binary.tree(tree) #check if it is binary (no politomies)
#tree<- multi2di(tree) #to make binary

is.ultrametric(tree) #to check if ultrametric

#min(tree$edge.length)# to check for 0's
#sum(tree$edge.length == 0)

Ge2014 <- read.csv("Ge_2014_initial.csv")
head(Ge2014)
length(Ge2014$Species)
rownames(Ge2014) <- Ge2014$Species #changing the rownames() in the traits data.frame to the tip labels
Ge2014 <- Ge2014[match(tree$tip.label, rownames(Ge2014)),]
cbind(Ge2014,tree$tip.label)

Body_length <- log(Ge2014[6])
head (Body_length)
write.csv(Body_length, file = "Body_length_BAMM_log.txt", row.names= TRUE)
write.tree(tree, "Tree_BAMM.tre")

setBAMMpriors(phy = tree, traits = "Body_length_BAMM_log.txt")
# correr o bam do terminal (bamm -c controlfile)

##### After BAMM results
setwd("/home/gabriela/Dropbox/Documents/Doutorado/Courses_conferences/ComBiol_Macroevolution/Lab/Projeto_Mike/Working_directory/Bamm_run1_rawvalues")
edata <- getEventData(tree, eventdata = "event_data.txt", burnin=0.1, type= "trait")
summary(edata)

mcmcout <- read.csv("mcmc_out.txt", header = T)
plot(mcmcout$logLik ~ mcmcout$generation)

burnstart <- floor(0.1 * nrow(mcmcout)) # Discard the first 10% of samples as burnin
postburn <- mcmcout[burnstart:nrow(mcmcout), ]
effectiveSize(postburn$N_shifts) 

effectiveSize(postburn$logLik) # 8322

Bayes_Factor <- computeBayesFactors("mcmc_out.txt", expectedNumberOfShifts = 1, burnin = 0.1)
plotPrior("mcmc_out.txt", expectedNumberOfShifts = 1)
#write.table(Bayes_Factor, file = "Table_ShiftBayesfactor.txt", col.names = TRUE, row.names = TRUE)

par(mfrow = c(1,2))
s <- plot.bammdata(edata, labels = T, font = 6, cex = 0.7, lwd = 3)
addBAMMlegend(s, location = "left", nTicks = 1)
axisPhylo()
mtext("time before present", side =1, adj = 0, padj= 4)
par(font = 1)
marg_probs <- marginalShiftProbsTree(edata)
plot.phylo(marg_probs, font= 6, cex = 0.7)
add.scale.bar(ask = TRUE, font = 0.9)
# add.scale.bar(ask=TRUE)

########################## test if sociality is correlated with body size
setwd ("/home/gabriela/Dropbox/Documents/Doutorado/Courses_conferences/ComBiol_Macroevolution/Lab/Projeto_Mike/Supplementary_material")

Tree_3genes <- read.nexus("3genes_final.tree")
Tree_3genes <- ladderize(Tree_3genes)
tree <- Tree_3genes
tree$tip.label
tree <- drop.tip(tree, tip= c("Ictidomys_mexicanus", "Spermophilus_alashanicus", "Spermophilus_dauricus", "Spermophilus_erythrogenys", "Spermophilus_relictus",   "Xerospermophilus_perotensis", "Notocitellus_adocetus", "Notocitellus_annulatus", "Otospermophilus_atricapillus", "Spermophilus_major", "Spermophilus_musicus", "Urocitellus_brunneus", "Urocitellus_undulatus", "Urocitellus_washingtoni", "Xerospermophilus_spilosoma", "Spermophilus_pallidicauda"))
length(tree$tip.label)

Ge2014_hib <- read.csv("Ge_2014_Bodysize_sociality.csv")
rownames(Ge2014_hib) <- Ge2014_hib$Species #changing the rownames() in the traits data.frame to the tip labels
length(Ge2014_hib$Species)
head (Ge2014_hib)
Ge2014_hib <- Ge2014_hib[match(tree$tip.label, rownames(Ge2014_hib)),]
cbind(Ge2014_hib, tree$tip.label)

Body_length <- Ge2014_hib$Body_length
names(Body_length) <- tree$tip.label
length(Body_length)

Socialidade <- read.csv("Sociality.txt", sep= ",", header = FALSE)
rownames(Socialidade) <- Socialidade$V1
Socialidade <- Socialidade[match(tree$tip.label, rownames(Socialidade)),]
cbind (Socialidade, tree$tip.label)
State <- Socialidade$V2 #V2 'e a classificacao ampla
State <- as.factor(State)
names(State) <- tree$tip.label

##### Phylogenetic Generalized Least Squares (PGLS)

data <- as.data.frame (cbind (Body_length, State))
pglsModel <- gls(Body_length ~ State, correlation = corBrownian(phy = tree), data=data, method = "ML")
summary(pglsModel)
#anova(pglsModel)
#coef(pglsModel)

par(mfrow = c(1,1))
plot(Body_length ~ State, xlab = "Sociality", ylab = "Body length") 
abline(a = coef(pglsModel)[1], b = coef(pglsModel)[2], lty = "solid", col = "red")

# following a OU model
pglsModelLambda <- gls(Body_length ~ State, correlation = corPagel(1, phy =tree, fixed = FALSE), data = data, method = "ML")
summary(pglsModelLambda)
plot(Body_length ~ State, xlab = "Sociality", ylab = "Body length")
abline(a = coef(pglsModelLambda)[1], b = coef(pglsModelLambda)[2])

##### Pra testar mais de uma coisa junta

tree <- Tree_3genes
tree <- drop.tip(tree, tip= c("Callospermophilus_madrensis", "Cynomys_gunnisoni", "Ictidomys_mexicanus", "Spermophilus_alashanicus", "Spermophilus_dauricus", "Spermophilus_erythrogenys", "Spermophilus_pallidicauda", "Spermophilus_relictus", "Urocitellus_elegans", "Urocitellus_townsendi", "Xerospermophilus_perotensis", "Notocitellus_adocetus", "Notocitellus_annulatus", "Otospermophilus_atricapillus", "Spermophilus_major", "Spermophilus_musicus", "Urocitellus_brunneus", "Urocitellus_undulatus", "Urocitellus_washingtoni", "Xerospermophilus_spilosoma", "Xerospermophilus_tereticaudus"))
length(tree$tip.label)

Ge2014_hib <- read.csv("Ge_2014_Hibernation_sociality.csv")
rownames(Ge2014_hib) <- Ge2014_hib$Species #changing the rownames() in the traits data.frame to the tip labels
Ge2014_hib <- Ge2014_hib[match(tree$tip.label, rownames(Ge2014_hib)),]
length(Ge2014_hib$Species)
cbind (tree$tip.label, Ge2014_hib)

Body_length <- Ge2014_hib$Body_length
names(Body_length) <- tree$tip.label
length(Body_length)

hibernation <- Ge2014_hib$Hibernation
names(hibernation) <- tree$tip.label
length(hibernation)

Socialidade <- read.csv("Sociality_hib_soc.csv", sep= ",", header = FALSE)
length(Socialidade$V1)
rownames(Socialidade) <- Socialidade$V1
Socialidade <- Socialidade[match(tree$tip.label, rownames(Socialidade)),]
cbind (tree$tip.label, Socialidade)
State <- Socialidade$V2 #V2 'e a classificacao ampla
State <- as.factor(State)
names(State) <- tree$tip.label

##### Phylogenetic Generalized Least Squares (PGLS)

library(ape)
library(geiger)
library(nlme)
library(phytools)

## Body vs sociality
data <- as.data.frame (cbind (Body_length, State))
pglsModel <- gls(Body_length ~ State, correlation = corBrownian(phy = tree), data=data, method = "ML")
summary(pglsModel)

#anova(pglsModel)
#coef(pglsModel)

par(mfrow = c(1,1))
plot(Body_length ~ State, , xlab = "Sociality", ylab = "Contrasts in ln(hibernation)", main = "Analysis using PICs") #mudar isso
abline(a = coef(pglsModel)[1], b = coef(pglsModel)[2])

# following a OU model
pglsModelLambda <- gls(Body_length ~ State, correlation = corPagel(1, phy =tree, fixed = FALSE), data = data, method = "ML")
summary(pglsModelLambda)
plot(Body_length ~ State, , xlab = "Sociality", ylab = "Contrasts in ln(hibernation)", main = "Analysis using PICs") #mudar isso
abline(a = coef(pglsModelLambda)[1], b = coef(pglsModelLambda)[2])


## Sociality_hibernation

cbind(tree$tip.label, hibernation)
data <- as.data.frame (cbind (hibernation, State))
pglsModel <- gls(hibernation ~ State, correlation = corBrownian(phy = tree), data=data, method = "ML")
summary(pglsModel)
#anova(pglsModel)
#coef(pglsModel)

par(mfrow = c(1,1))
plot(hibernation ~ State, , xlab = "Sociality", ylab = "Contrasts in ln(hibernation)", main = "Analysis using PICs") #mudar isso
abline(a = coef(pglsModel)[1], b = coef(pglsModel)[2])

# following a OU model
pglsModelLambda <- gls(hibernation ~ State, correlation = corPagel(1, phy =tree, fixed = FALSE), data = data, method = "ML")
summary(pglsModelLambda)
plot(hibernation ~ State, , xlab = "Sociality", ylab = "Contrasts in ln(hibernation)", main = "Analysis using PICs") #mudar isso
abline(a = coef(pglsModelLambda)[1], b = coef(pglsModelLambda)[2])


## Body size vs hibernation

data <- as.data.frame (cbind (hibernation, Body_length))
pglsModel <- gls(hibernation ~ Body_length, correlation = corBrownian(phy = tree), data=data, method = "ML")
summary(pglsModel)
#anova(pglsModel)
#coef(pglsModel)

par(mfrow = c(1,1))
plot(hibernation ~ Body_length, , xlab = "Body length", ylab = "Hibernation length") #mudar isso
abline(a = coef(pglsModel)[1], b = coef(pglsModel)[2], lty = "solid", col = "red")
abline(a = coef(pglsModelLambda)[1], b = coef(pglsModelLambda)[2], lty = "solid", col = "blue")

# following a OU model
pglsModelLambda <- gls(hibernation ~ Body_length, correlation = corPagel(1, phy =tree, fixed = FALSE), data = data, method = "ML")
summary(pglsModelLambda)
plot(hibernation ~ Body_length, , xlab = "Sociality", ylab = "Contrasts in ln(hibernation)", main = "Analysis using PICs") #mudar isso
abline(a = coef(pglsModelLambda)[1], b = coef(pglsModelLambda)[2])

#### todo mundo junto
data <- as.data.frame (cbind (hibernation, Body_length, State))
pglsModel3 <- gls(Body_length ~ State * hibernation, correlation = corBrownian(phy = tree), data = data, method = "ML")
anova(pglsModel3)
summary (pglsModel3)

data <- as.data.frame (cbind (hibernation, Body_length, State))
pglsModel4 <- gls(Body_length ~ State * hibernation, correlation = corPagel(1, phy =tree, fixed = FALSE), data = data, method = "ML")
anova(pglsModel4)
summary (pglsModel4)


### Body size and hibernation are not related, but if I remove the species that store food?


tree <- Tree_3genes
tree <- drop.tip(tree, tip= c("Callospermophilus_madrensis", "Cynomys_gunnisoni", "Ictidomys_mexicanus", "Spermophilus_alashanicus", "Spermophilus_dauricus", "Spermophilus_erythrogenys", "Spermophilus_pallidicauda", "Spermophilus_relictus", "Urocitellus_elegans", "Urocitellus_townsendi", "Xerospermophilus_perotensis", "Notocitellus_adocetus", "Notocitellus_annulatus", "Otospermophilus_atricapillus", "Spermophilus_major", "Spermophilus_musicus", "Urocitellus_brunneus", "Urocitellus_undulatus", "Urocitellus_washingtoni", "Xerospermophilus_spilosoma", "Xerospermophilus_tereticaudus", "Ammospermophilus_harrisii", "Cynomys_ludovicianus", "Cynomys_mexicanus", "Otospermophilus_beecheyi", "Otospermophilus_variegatus", "Poliocitellus_franklinii", "Spermophilus_fulvus", "Spermophilus_suslicus", "Spermophilus_xanthoprymnus"))
length(tree$tip.label)

Ge2014_hib <- read.csv("Ge_2014_Hibernation_body_store.csv")
rownames(Ge2014_hib) <- Ge2014_hib$Species #changing the rownames() in the traits data.frame to the tip labels
Ge2014_hib <- Ge2014_hib[match(tree$tip.label, rownames(Ge2014_hib)),]
length(Ge2014_hib$Species)
cbind (tree$tip.label, Ge2014_hib)

Body_length <- Ge2014_hib$Body_length
names(Body_length) <- tree$tip.label
length(Body_length)

hibernation <- Ge2014_hib$Hibernation
names(hibernation) <- tree$tip.label
length(hibernation)

data <- as.data.frame (cbind (hibernation, Body_length))
pglsModel <- gls(hibernation ~ Body_length, correlation = corBrownian(phy = tree), data=data, method = "ML")
summary(pglsModel)
#anova(pglsModel)
#coef(pglsModel)

par(mfrow = c(1,1))
plot(hibernation ~ Body_length, , xlab = "Sociality", ylab = "Contrasts in ln(hibernation)", main = "Analysis using PICs") #mudar isso
abline(a = coef(pglsModel)[1], b = coef(pglsModel)[2])

# following a OU model
pglsModelLambda <- gls(hibernation ~ Body_length, correlation = corPagel(1, phy =tree, fixed = FALSE), data = data, method = "ML")
summary(pglsModelLambda)
plot(hibernation ~ Body_length, , xlab = "Sociality", ylab = "Contrasts in ln(hibernation)", main = "Analysis using PICs") #mudar isso
abline(a = coef(pglsModelLambda)[1], b = coef(pglsModelLambda)[2])
