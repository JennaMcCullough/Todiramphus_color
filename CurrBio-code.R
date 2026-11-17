# Entire code for McCullough et al. 2026 
#
# Written by Jenna McCullough 
# 17 November 2026 

# There are five main parts to this code: 
# 1. Calculating diversification with CLaDS and ES metrics. 
# 
# 2. Creation of excel files from spectra that have averaged X,Y,Z coordinates 
# that represent the average location of each specimen of a taxon in 
# colorspace. Each value is averaged across multiple specimens. This is 
# calculated from the specs4_metaclean4-final-extra-cleaning.rda file that is
# available on Dryad. 
#
# 3. Estimating color evolutionary rates with RRphylo. 
#
# 4. Preparing the adaptive explanatory variables: A) calculating dichromatism 
# as well as testing whether colordist vs. Euclidian distances give us 
# the same general answers. B) calculating island size metrics.  
# 
# 5. Testing drivers of rapid color evolution with linear and phylogenetic 
# linear models. We will also plot these (representing many of the plots in
# figures S7 <- S14. 
#
# 6. Sensitivity tests for heterozygosity. We specifically will look at 
# different ways to account for the potential bias that short branches may
# have when trying to understand the rate of trait changes. We address three 
# questions: 1. Are heterozygosity estimates are uncorrelated with taxon age
# (i.e. terminal branch length in the tree)? 2. Do relationships persist when
# using colour rates derived from lambda-scaled trees? 3. Do relationships 
# persist when including taxon age (i.e. terminal branch length in the tree) as
# an additional predictor in the multipredictor models to account for the fact
# that short branches are invariably inferred to have faster rates? 

# note: Sometimes I don't produce the best, nicest code. Sometimes I output a 
# CSV file from R and then manually add it into an excel file. Bear with me! 
# I also will comment out analyses that take a long time just in case someone 
# might accidentally run it. I'll tell you when that happens. 
# Also for all generation of figures, I always have a hashtag before the PDF. 
# This is because I don't want a figure to accidentally get overwritten. 
# Uncomment and run at your leisure 


# lots of libraries to load for all the analyses 
library(pavo)
library(RRphylo)
library(multcompView)
library(phytools)
library(ape)
library(ggplot2)
library(viridis)
library(terra)
library(dplyr)
library(grDevices)
library(gridExtra)
library(geomorph)
library(treeplyr)
library(phylolm)
library(cowplot)
library(mvMORPH)
library(RPANDA)
library(coda)
library(phangorn)
library(RColorBrewer)
library(fields)
library(maps)
library(gdm)
library(tidyr)
library(pbapply)


#################################################################
# PART 1: Calculating diversification with CLaDS and ES metrics #
#################################################################

phy <- read.nexus("Todi-timetree.trans.nex") 

# For both with and without extinction, we use an inverse gamma prior with shape
#parameter 1 and rate parameter 0.1 for σ, and a flat prior for all other 
# parameters.

#This took several days. Uncomment and run at your leisure 
#sampler <- fit_ClaDS(phy, sample_fraction=1, iterations=200000, thin = 100, 
#                    file_name = "mcmc_phy-200K-thin100.rda", model_id="ClaDS2", nCPU = 6)

# I included the results for this in diversification-output/ 

load("diversification-output/mcmc_phy-200K-thin100.rda")
# Get mean posterior diversification rates
maps <- getMAPS_ClaDS(mcmcSampler, thin = 1)
cat(paste0("sigma = ", maps[1],
           "\n alpha = ", maps[2],
           "\n epsilon = ", maps[3],
           "\n l_0 = ", maps[4], "\n"))

# Output figure. I added finishing touches in Adobe Illustrator. 

#pdf("cladsphylo-200Kgen-V1.pdf", width=8, height=8)
par(mar=c(5,0,0,0))
plot_ClaDS_phylo(phy, maps[-(1:4)], show.tip=T, cex=.5, type="fan")
plot_ClaDS_phylo(phy, maps[-(1:4)], show.tip=T, cex=.5)
#dev.off()

# Make vector with diversification rates for use later
specCLADS <- maps[-(1:4)][match(1:Ntip(phy), phy$edge[,2])]
names(specCLADS) <- phy$tip
head(specCLADS)

# Also let's calculate equal splits (ES) speciation metric
specES <- compute_es(phy)

ggplot(data=NULL, aes(x=specCLADS, y=specES)) +
  geom_point() + 
  # ggrepel::geom_text_repel(aes(label=names(specES)), cex=1) + 
  geom_abline(slope=1, intercept=0, lty=2)
# Positive relationship, but lots of scatter..

#make csv files
#write.csv(specES, file="todi96-divrates-ES.csv")
#write.csv(specCLADS, file="todi96-divrates-clads.csv")

# I then added these manually to the file that I will plot and do lms with, 
# color-data/sex-difference-colrate-lumrate.csv 
# This is also available in supplementary table S7! 

# We are going to be using this file a lot throughout the code.  
sex_diffs <- read.csv("color-data/sex-difference-colrate-lumrate.csv")

# comparing the two, this creates figure S16 
#pdf("es-clads-comparison-unlabelled.pdf", width=6.5, height=4)
ggplot(data=sex_diffs, aes(x=clads, y=es)) +
  geom_point(size=3, stroke = 0.5, color = sex_diffs$sp_col) +
  geom_abline(slope=1, intercept=0, lty=2) 
#dev.off() 

# ok we have our diversification results, now lets do color. 

################################################
# PART 2: Calculating color evolutionary rates #
################################################

# The raw color data was spread across 12 different batches and two additional
# published datasets (Eliason et al. 2019, 2023). Additionally, I personally 
# looked at all 53,000 spectra and had to do some quality control on those 
# spectra. Sometimes a patch would have been mis-ordered or we needed to drop 
# some wonky averaged spectra. It was 2,000 individual lines of code that I will
# not include here. If this is something that you would like, I would be happy
# to give it to you. But I don't think its nessary for replication. the dataset
# think its necessary to replicate the analyses In my 
# opinion, it is superfluous for anyone who'd want to replicate these analyses. 
# Therefore, I am not including the raw dataset here. 

# Instead, I am providing the cleaned dataset with averaged spectra. This is 
# the dataset that all the color analyses were performed on. We will use this R
# data file, which includes the spectra and the associated metadata. It is 
# available on Dryad 

# Make chromatic and achromatic data sets for males and females from the
# averaged spectra. This final dataset went through quality control and some
# of the averaged scans with issues were dropped.  

# you can download the r data file on dryad. Then include it into your working
# directory. 

load("specs4_metaclean4-final-extra-cleaning.rda")

# you can see all these data like this! 
# Plot spectra according to their human-perceived colors
#plot(specs4, col=spec2rgb(specs4))
#str(specs4)

# Plot spectra in avian tetrahedral colorspace 3D coordinates (x,y,z)
# one of these replicates figure 1A
plot(x~y, colordist, col=spec2rgb(specs4), pch=16)
plot(x~z, colordist, col=spec2rgb(specs4), pch=16)
plot(y~z, colordist, col=spec2rgb(specs4), pch=16)

#make averaged dataset for males and females for macroevolutionary analyses.
vm.uv <- vismodel(specs4, visual="avg.uv", relative=FALSE, achro="bt.dc")
X <- colspace(vm.uv)

X$spp <- meta_clean4$Taxon[match(rownames(X), meta_clean4$specID)]
X$sex <- meta_clean4$Sex[match(rownames(X), meta_clean4$specID)]
X$patch <- meta_clean4$V2[match(rownames(X), meta_clean4$specID)]

# average XYZ coordinates across specimens to get an average X, Y, and Z value
# for each patch and each species. There's issue with missing data though. This
# was an issue for three taxa: T. lazuli, T. tristrami bennetti, and T. chloris
# davisoni. 

# First We are extracting the average XYZ coordinates (chromatic dataset) for 
# all the patches here
chrom_m <- X %>% group_by(spp,patch) %>% filter(sex=="m") %>%
  select(x,y,z) %>%
  summarize_if(is.numeric, mean) %>%
  gather(coord,value,-spp,-patch) %>%
  unite(newsex,coord,patch) %>%
  spread(newsex, value) %>% na.omit

chrom_f <- X %>% group_by(spp,patch) %>% filter(sex=="f") %>%
  select(x,y,z) %>%
  summarize_if(is.numeric, mean) %>%
  gather(coord,value,-spp,-patch) %>%
  unite(newsex,coord,patch) %>%
  spread(newsex, value) %>% na.omit

#write.csv(chrom_f, file="tax-f-avg-specs3.xyz.csv", row.names = FALSE)
#write.csv(chrom_m, file="tax-m-avg-specs3.xyz.csv", row.names = FALSE)

# however, because two taxa do not have any individuals for the lore patch, this
# function drops those two taxa from the dataset! Lets make another dummy 
# dataset so we can add their 14 patches manually them in manually and then put
# "NA" in their lore patch.
X_noLores <- X
X_noLores = X_noLores[!grepl("1_lores", X_noLores$patch),]

chrom_m_Nolores <- X_noLores %>% group_by(spp,patch) %>% filter(sex=="m") %>%
  select(x,y,z) %>%
  summarize_if(is.numeric, mean) %>%
  gather(coord,value,-spp,-patch) %>%
  unite(newsex,coord,patch) %>%
  spread(newsex, value) %>% na.omit

chrom_f_Nolores <- X_noLores %>% group_by(spp,patch) %>% filter(sex=="f") %>%
  select(x,y,z) %>%
  summarize_if(is.numeric, mean) %>%
  gather(coord,value,-spp,-patch) %>%
  unite(newsex,coord,patch) %>%
  spread(newsex, value) %>% na.omit

#write.csv(chrom_f_Nolores, file="output/tax-f-avg-specs3-nolores.xyz.csv", row.names = TRUE)
#write.csv(chrom_m_Nolores, file="output/tax-m-avg-specs3-nolores.xyz.csv", row.names = TRUE)

# in excel, I added those two taxa manually. I suggest just using the file I 
# included on dryad (and github). I'm specifically referring to the files
# with names like "tax95-f-achrom-avg-specs3-nolores.lum.csv"  

# I want to also extract the achromatic variable "luminance" to include for
# RRphylo. Going to do both the no lores and with lores as above. 
achrom_m_noLores <- X_noLores%>% group_by(spp,patch) %>% filter(sex=="m") %>%
  select(lum) %>%
  summarize_if(is.numeric, mean) %>%
  gather(coord,value,-spp,-patch) %>%
  unite(newsex,coord,patch) %>%
  spread(newsex, value) %>% na.omit

achrom_f_noLores <- X_noLores %>% group_by(spp,patch) %>% filter(sex=="f") %>%
  select(lum) %>%
  summarize_if(is.numeric, mean) %>%
  gather(coord,value,-spp,-patch) %>%
  unite(newsex,coord,patch) %>%
  spread(newsex, value) %>% na.omit
#write.csv(achrom_f_noLores, file="output/tax-f-avg-specs3-nolores.lum.csv", row.names = TRUE)
#write.csv(achrom_m_noLores, file="output/tax-m-avg-specs3-nolores.lum.csv", row.names = TRUE)

achrom_m_allpatch <- X%>% group_by(spp,patch) %>% filter(sex=="m") %>%
  select(lum) %>%
  summarize_if(is.numeric, mean) %>%
  gather(coord,value,-spp,-patch) %>%
  unite(newsex,coord,patch) %>%
  spread(newsex, value) %>% na.omit

achrom_f_allpatch <- X %>% group_by(spp,patch) %>% filter(sex=="f") %>%
  select(lum) %>%
  summarize_if(is.numeric, mean) %>%
  gather(coord,value,-spp,-patch) %>%
  unite(newsex,coord,patch) %>%
  spread(newsex, value) %>% na.omit
#write.csv(achrom_f_allpatch, file="output/tax-f-avg-specs3-allpatch.lum.csv", row.names = TRUE)
#write.csv(achrom_m_allpatch, file="output/tax-m-avg-specs3-allpatch.lum.csv", row.names = TRUE)

############################################################
# PART 3: Estimating color evolutionary rates with RRphylo #
############################################################ 

tre <- read.nexus("Todi-timetree.trans.nex") 
plot(tre, cex=0.5)

# estimating color evolutionary rates using RRphylo. As discussed in the 
#methods, we did not analyze the white headed morphotypes of two saurophagus 
#subspecies, anachoreta and admiralitatis. The males and females datasets also
# differ because T. sacer aneityum does not have any female specimens. We 
#have to use different numbered phylogenies because of this. 

# the way these are named overwrites the last analysis. Be aware. 

# Chromatic colors (xyz coordinates)
# females first 
fem95_14_xyz <- read.csv("color-data/tax95-f-chrom-avg-specs3.xyz.nolores.csv")
#i've already removed the lores and the white-headed saurophagus morphotypes
rownames(fem95_14_xyz) <- fem95_14_xyz$spp
fem95_14_xyz = subset(fem95_14_xyz, select = -c(spp) )
# we don't have data for female Todiramphus_sacer_aneityum
tre_f95 <- drop.tip(tre,"Todiramphus_sacer_aneityum", rooted = is.rooted(tre))

RR <- RRphylo(tree=tre_f95, y=fem95_14_xyz)
pRR <- plotRR(RR,y=fem95_14_xyz)
id <- c((Ntip(tre_f95)+1):(Ntip(tre_f95)+Nnode(tre_f95)), 1:Ntip(tre_f95))
erates <- RR$rates[,1][match(tre_f95$edge[,2], id)]
res <- search.shift(RR, status.type="clade")
res

#producing the figure that are eventually included in Figure S4 
#pdf(file="rrphylo.14patch.xyzcoord.rates.pdf",width=10, height=14)
pRR$plotRRrates(tree.args=list(edge.width=4),color.pal=viridis,
                colorbar.args = list(x="bottomleft",direction="h",labs.adj=0.7,xpd=TRUE))
#dev.off()

#If you wanted the phenotype, which I'm not really sure what that means here, 
# you could run this. This isn't included in the study though. 
#pdf(file="rrphylo.14patch.xyzcoord.phen.pdf",width=10, height=14)
pRR$plotRRphen(variable=1,tree.args=list(edge.width=4),color.pal=viridis,
               colorbar.args = list(x="bottomleft",labs.adj=0.7,xpd=TRUE))
#dev.off()

# I downloaded the rates like this 
#write.csv(erates, "14patch.f.xyzcoord.erates.csv")

# the appropriate data was combined in a subsequent file by hand in excel. It's
# in our sex_diffs$colrate.female column 

#male chromatic color. 
mal96_14_xyz <- read.csv("color-data/tax96-m-chrom-avg-specs3.xyz.nolores.csv")
rownames(mal96_14_xyz) <- mal96_14_xyz$spp
mal96_14_xyz = subset(mal96_14_xyz, select = -c(spp) )

#we don't drop aneityum here so we use the full tree 
RR <- RRphylo(tree=tre, y=mal96_14_xyz)
pRR <- plotRR(RR,y=mal96_14_xyz)
id <- c((Ntip(tre)+1):(Ntip(tre)+Nnode(tre)), 1:Ntip(tre))
erates <- RR$rates[,1][match(tre$edge[,2], id)]
res <- search.shift(RR, status.type="clade")
res
#another phenotype version of plotting if you're interested. 
#pdf(file="rrphylo.14patch.m.xyzcoord.phen.pdf",width=10, height=14)
pRR$plotRRphen(variable=1,tree.args=list(edge.width=4),color.pal=viridis,
               colorbar.args = list(x="bottomleft", labs.adj=0.7,xpd=TRUE))
#dev.off()

#This produces a part of figure S3. 
#pdf(file="rrphylo.14patch.m.xyzcoord.rates.pdf",width=10, height=14)
pRR$plotRRrates(tree.args=list(edge.width=4),color.pal=viridis,
                colorbar.args = list(x="bottomleft",direction="h",labs.adj=0.7,xpd=TRUE))
#dev.off()

#This produces the radial phylogeny from figure 3. I added all the kingfishers. 
#pdf(file="rrphylo.14patch.m.xyzcoord.rates.FAN.pdf",width=10, height=10)
pRR$plotRRrates(tree.args=list(edge.width=4, type="fan"),color.pal=viridis,
                colorbar.args = list(x="bottomleft",direction="h",labs.adj=0.7,xpd=TRUE))
#dev.off()

#saving those rates then I added them to our trusty sex_diffs$colrate.male 
#write.csv(erates, "14patch.m.xyzcoord.erates.csv")

# ACHROMATIC colors 
#females first. 
fem95_14_lum <- read.csv("color-data/tax95-f-achrom-avg-specs3-nolores.lum.csv")
rownames(fem95_14_lum) <- fem95_14_lum$spp
fem95_14_lum = subset(fem95_14_lum, select = -c(spp) )

RR <- RRphylo(tree=tre_f95, y=fem95_14_lum)
pRR <- plotRR(RR,y=fem95_14_lum)
id <- c((Ntip(tre_f95)+1):(Ntip(tre_f95)+Nnode(tre_f95)), 1:Ntip(tre_f95))
erates <- RR$rates[,1][match(tre_f95$edge[,2], id)]
res <- search.shift(RR, status.type="clade")
res
#producing the other half of figure S4. 
#pdf(file="rrphylo.14patch.f.lum.rates.pdf",width=10, height=14)
pRR$plotRRrates(tree.args=list(edge.width=4),color.pal=viridis,
                colorbar.args = list(x="bottomleft",direction="h",labs.adj=0.7,xpd=TRUE))
#dev.off()
# adding the rates to sex_diffs$lumrate.female by hand. 
#write.csv(erates, "14patch.f.luminance.erates.csv")

#male achromatic colors
mal96_14_lum <- read.csv("color-data/tax96-m-achrom-avg-specs3-nolores.lum.csv")
rownames(mal96_14_lum) <- mal96_14_lum$spp
mal96_14_lum = subset(mal96_14_lum, select = -c(spp) )

RR <- RRphylo(tree=tre, y=mal96_14_lum)
pRR <- plotRR(RR,y=mal96_14_lum)
id <- c((Ntip(tre)+1):(Ntip(tre)+Nnode(tre)), 1:Ntip(tre))
erates <- RR$rates[,1][match(tre$edge[,2], id)]
res <- search.shift(RR, status.type="clade")
res
# this produces a part of figure S3. 
#pdf(file="rrphylo.14patch.m.lum.rates.pdf",width=10, height=14)
pRR$plotRRrates(tree.args=list(edge.width=4),color.pal=viridis,
                colorbar.args = list(x="bottomleft",direction="h",labs.adj=0.7,xpd=TRUE))
#dev.off()
# adding the rates to sex_diffs$lumrate.male by hand. 
#write.csv(erates, "14patch.m.luminance.erates.csv")

#######################################################################
# lms and plotting the relationship between color and diversification # 
#######################################################################

#we are again using the final dataset here for simplicity.

#rate of color evolution against rate of diversification. 
# The sex_diffs$sp_col column gives a color to each of the species level 
# taxonomy. I didn't like the color so I changed them in Illustrator. You get
# the gist in this color scheme though. 

#if you want to know the identity of any point, you can add this line to the
# ggplot box of text. I don't do it for each of them for simplicity though. 
# add this with no hashtag: 
# + geom_text(aes(label = branch), size = 1)

# Male chromatic and achromatic rates by diversification metric. We also do 
# the linear models (which are included in Table S1). 

#Figure S5 
S5A_lm <- lm(colrate.male~clads, data=sex_diffs)
summary(S5A_lm)
coef(S5A_lm)
S5A <- ggplot(data=sex_diffs, aes(y = colrate.male, x = clads)) + geom_point() + 
  geom_point(color = sex_diffs$sp_col,  size=3)+
  geom_smooth(method="lm")  +
  ggtitle("chromatic rate vs. clads diversification rate male ")

S5B_lm <- lm(lumrate.male~clads, data=sex_diffs)
summary(S5B_lm)
coef(S5B_lm)
S5B <- ggplot(data=sex_diffs, aes(y = lumrate.male, x = clads)) + geom_point() + 
  geom_point(color = sex_diffs$sp_col,  size=3)+
  geom_smooth(method="lm")  +
  ggtitle("achromatic rate vs. clads diversification rate male ")

S5C_lm <- lm(colrate.male~es, data=sex_diffs)
summary(S5C_lm)
coef(S5C_lm)
S5C <- ggplot(data=sex_diffs, aes(y = colrate.male, x = es)) + geom_point() + 
  geom_point(color = sex_diffs$sp_col,  size=3)+
  geom_smooth(method="lm")  +
  ggtitle("chromatic rate vs. es diversification rate male ")

S5D_lm <- lm(lumrate.male~es, data=sex_diffs)
summary(S5D_lm)
coef(S5D_lm)
S5D <- ggplot(data=sex_diffs, aes(y = lumrate.male, x = es)) + geom_point() + 
  geom_point(color = sex_diffs$sp_col,  size=3)+
  geom_smooth(method="lm")  +
  ggtitle("achromatic rate vs. es diversification rate male ")

#pdf(file="S5_clads_v_ES_male.pdf", width=7, height=9)
grid.arrange(S5A, S5B, S5C, S5D, ncol = 2)
#dev.off()

#Figure S6, females!
S6A_lm <- lm(colrate.female~clads, data=sex_diffs)
summary(S6A_lm)
coef(S6A_lm)
S6A <- ggplot(data=sex_diffs, aes(y = colrate.female, x = clads)) + geom_point() + 
  geom_point(color = sex_diffs$sp_col,  size=3)+
  geom_smooth(method="lm")  +
  ggtitle("chromatic rate vs. clads diversification rate female ")

S6B_lm <- lm(lumrate.female~clads, data=sex_diffs)
summary(S6B_lm)
coef(S6B_lm)
S6B <- ggplot(data=sex_diffs, aes(y = lumrate.female, x = clads)) + geom_point() + 
  geom_point(color = sex_diffs$sp_col,  size=3)+
  geom_smooth(method="lm")  +
  ggtitle("achromatic rate vs. clads diversification rate female ")

S6C_lm <- lm(colrate.female~es, data=sex_diffs)
summary(S6C_lm)
coef(S6C_lm)
S6C <- ggplot(data=sex_diffs, aes(y = colrate.female, x = es)) + geom_point() + 
  geom_point(color = sex_diffs$sp_col,  size=3)+
  geom_smooth(method="lm")  +
  ggtitle("chromatic rate vs. es diversification rate female ")

S6D_lm <- lm(lumrate.female~es, data=sex_diffs)
summary(S6D_lm)
coef(S6D_lm)
S6D <- ggplot(data=sex_diffs, aes(y = lumrate.female, x = es)) + geom_point() + 
  geom_point(color = sex_diffs$sp_col,  size=3)+
  geom_smooth(method="lm")  +
  ggtitle("achromatic rate vs. es diversification rate female ")

#pdf(file="S6_clads_v_ES_female.pdf", width=7, height=9)
grid.arrange(S6A, S6B, S6C, S6D, ncol = 2)
#dev.off()

#########################################################
# PART 4:  Preparing the adaptive explanatory variables #
#########################################################

# We have three general explanations for color rates. 
#1. Sexual selection, as measured by sexual dichromatism 
#2. Competition, as measured by sympatry 
#3. Population size, measured by island size and heterozygosity 

#1. CALCULATING CHROMATIC AND ACHROMATIC DICHROMATISM FROM SPECTRAL DATA

#these are R functions written by Chad Eliason 
source("R/mod_pavo_functions.R")
source("R/coldist_FAST.r")

#in case this isn't loaded already 
load("specs4_metaclean4-final-extra-cleaning.rda")

#the names of our patches 
patches <- c('1_lores','2_crown','3_eyebrow','4_mustache','5_nape','6_mantle','7_back','8_rump','9_tail','10_coverts','11_wing','12_throat','13_breast','14_belly','15_vent')

#create vector of taxon names. 
Taxon_vec <- c(t(unique(meta_clean4$Taxon)))
## i want to remove NA for the wl column
Taxon_vec <- Taxon_vec[-c(1)]
is.vector(Taxon_vec)

id_vec <- meta_clean4[,1]
sex_vec <- meta_clean4[,3]
patch_vec <- meta_clean4[,2]
taxon_vec <- meta_clean4[,4]

# First, we are calculating dichromatism in a receptor-noise-based 
# receptor-noise model of Vorobyev et al. (1998) to calculate color distances 
# (measured as 'just noticeable distances, JNDs). This will allow us to get 
# dichromatism both in chromatic and achromatic colors.

# make a visual model from the specs. 
vm.uv <- vismodel(specs3, visual="avg.uv", relative=FALSE, achro="bt.dc")

# what are the distances between points? 
# this takes about 8 minutes to run:
coldist_patches <- pblapply(patches, function(x) {
  # x=patches[1]
  cd <- pavo::coldist(subset(vm.uv, meta_clean4$V2==x), achromatic=TRUE)
  # sanity checks-
  # cd[grep("NA", cd$patch1),]
  # setdiff(cd$patch1, newmeta$specID)
  # setdiff(cd$patch2, newmeta$specID)
  cd$sex1 <- meta_clean4$Sex[match(cd$patch1, meta_clean4$specID)]
  cd$sex2 <- meta_clean4$Sex[match(cd$patch2, meta_clean4$specID)]
  cd$taxon1 <- meta_clean4$Taxon[match(cd$patch1, meta_clean4$specID)]
  cd$taxon2 <- meta_clean4$Taxon[match(cd$patch2, meta_clean4$specID)]
  cd <- separate(cd, patch1, into=c('id1','patch1'), sep="\\.")
  cd <- separate(cd, patch2, into=c('id2','patch2'), sep="\\.")
  cd_sub <- subset(cd, (sex1=="m"&sex2=="f" | sex1=="f"&sex2=="m") & taxon1 == taxon2)
  cd_sub %>% group_by(taxon1) %>% summarize(dichrom_color = mean(dS), dichrom_lum = mean(dL))
})
#save(coldist_patches,file="coldist.rda")

# getting average of all the patches across the body: 
# I'm not always the fanciest coder. Some the patches are missing species so I
# just decided to do it by hand. 

# we have 15 tibbles of patches 
patch1 <- coldist_patches[[1]]
patch2 <- coldist_patches[[2]]
patch3 <- coldist_patches[[3]]
patch4 <- coldist_patches[[4]]
patch5 <- coldist_patches[[5]]
patch6 <- coldist_patches[[6]]
patch7 <- coldist_patches[[7]]
patch8 <- coldist_patches[[8]]
patch9 <- coldist_patches[[9]]
patch10 <- coldist_patches[[10]]
patch11 <- coldist_patches[[11]]
patch12 <- coldist_patches[[12]]
patch13 <- coldist_patches[[13]]
patch14 <- coldist_patches[[14]]
patch15 <- coldist_patches[[15]]
# the numbers of taxa in these subsets differ 
# patch 1 is because two subspecies don't have lores 
# but i don't know which ones at this point 
nrow(patch1)#95 taxa
nrow(patch2)#97 taxa
nrow(patch3)#97 taxa
nrow(patch4)#97 taxa
nrow(patch5)#97 taxa
nrow(patch6)# 96 taxa
nrow(patch7)# 96 taxa
nrow(patch8)#97 taxa
nrow(patch9)#97 taxa
nrow(patch10)# 96 taxa
nrow(patch11)#97 taxa
nrow(patch12)#97 taxa
nrow(patch13)#97 taxa
nrow(patch14)# 96 taxa
nrow(patch15)# 96 taxa
# missing: 
# bennetti and lazuli in patch 1 lores only
# Kalbaensis in patches 6,7,13,14
# stictolaemus in patch 10 coverts only 

# uncomment when you want to run 
#write.csv(coldist_patches[[1]],  file="color-data/coldist_patches_receptor-noise/patch1.csv")
#write.csv(coldist_patches[[2]],  file="color-data/coldist_patches_receptor-noise/patch2.csv")
#write.csv(coldist_patches[[3]],  file="color-data/coldist_patches_receptor-noise/patch3.csv")
#write.csv(coldist_patches[[4]],  file="color-data/coldist_patches_receptor-noise/patch4.csv")
#write.csv(coldist_patches[[5]],  file="color-data/coldist_patches_receptor-noise/patch5.csv")
#write.csv(coldist_patches[[6]],  file="color-data/coldist_patches_receptor-noise/patch6.csv")
#write.csv(coldist_patches[[7]],  file="color-data/coldist_patches_receptor-noise/patch7.csv")
#write.csv(coldist_patches[[8]],  file="color-data/coldist_patches_receptor-noise/patch8.csv")
#write.csv(coldist_patches[[9]],  file="color-data/coldist_patches_receptor-noise/patch9.csv")
#write.csv(coldist_patches[[10]],  file="color-data/coldist_patches_receptor-noise/patch10.csv")
#write.csv(coldist_patches[[11]],  file="color-data/coldist_patches_receptor-noise/patch11.csv")
#write.csv(coldist_patches[[12]],  file="color-data/coldist_patches_receptor-noise/patch12.csv")
#write.csv(coldist_patches[[13]],  file="color-data/coldist_patches_receptor-noise/patch13.csv")
#write.csv(coldist_patches[[14]],  file="color-data/coldist_patches_receptor-noise/patch14.csv")
#write.csv(coldist_patches[[15]], file="color-data/coldist_patches_receptor-noise/patch15.csv")

# all of these are available in color-data/coldist_patches_receptor-noise/ 
# I combined them by hand into a single excel file... sometimes I just like 
# making more work for myself. I have the final excel file 15-patches-combined.
# csv also in that folder structure. In that final file, I averaged values to
# make a single value of chromatic ('color dichromatism') and achromatic 
# ('luminance dichromatism') dichromatism. These are the columns 
# "AvgColorDichromatism" and "AvgLuminanceDichromatism"from Table S7. I 
# shortened these names for the sex_diffs file: sex_diffs$AvgDicCol and 
# sex_diffs$AvgDicLum. 

# COMPARISON WITH EUCLIDIAN DISTANCES ESTIMATED FROM TETRAHEDRAL COLORSPACE 
# But is the way I calculated dichromatism the *best* way? Is calculating 
# dichromatism in a receptor-noise model of Vorobyev et al. (1998) an equitable
# comparison if I used XYZ coordinates for the RRphylo analyses? We could 
# compare another way to estimate dichromatism by using the Euclidian distances
# between points in tetrahedral colorspace. However, this method can't capture
# achromatic differences. 

# question: are the two ways of calculating chromatic dichromatism correlated?
# we calculate a different visual model, this time with relative as TRUE. 
vm.uvREL <- vismodel(specs4, visual="avg.uv", relative=TRUE)
tcs.REL <- colspace(vm.uvREL, space = "tcs")

#are these the same length? 
all(rownames(tcs.REL) == meta_clean4$specID)
head(rownames(tcs.REL))
head(meta_clean4$specID)
#no. Lets drop the wl column 
meta_REL <- subset(meta_clean4, specID != "wl")
tcs.REL_clean <- tcs.REL[rownames(tcs.REL) %in% meta_REL$specID,]
meta_REL_clean <- meta_REL[match(rownames(tcs.REL_clean), meta_REL$specID),]
all(rownames(tcs.REL_clean) == meta_REL_clean$specID)
sum(is.na(meta_REL_clean$specID))

# ok now we can calculate differences in the sexes. 
coldist_patchesREL <- pblapply(patches, function(x) {
  tcs_sub <- subset(tcs.REL_clean, meta_REL_clean$V2 == x)
  pairs <- t(combn(seq_len(nrow(tcs_sub)), 2))
  distREL <- data.frame(   patch1 = rownames(tcs_sub)[pairs[,1]],   patch2 = rownames(tcs_sub)[pairs[,2]] )
  distREL$dS <- sqrt((tcs_sub$x[pairs[,1]] - tcs_sub$x[pairs[,2]])^2 +    (tcs_sub$y[pairs[,1]] - tcs_sub$y[pairs[,2]])^2 +    (tcs_sub$z[pairs[,1]] - tcs_sub$z[pairs[,2]])^2)
  distREL$sex1 <- meta_REL_clean$Sex[   match(distREL$patch1, meta_REL_clean$specID) ]
  distREL$sex2 <- meta_REL_clean$Sex[   match(distREL$patch2, meta_REL_clean$specID) ]
  distREL$taxon1 <- meta_REL_clean$Taxon[    match(distREL$patch1, meta_REL_clean$specID)]
  distREL$taxon2 <- meta_REL_clean$Taxon[match(distREL$patch2, meta_REL_clean$specID)]
  distREL <- separate(  distREL,  patch1,  into = c("id1","patch1"),  sep = "\\.")
  distREL <- separate( distREL, patch2, into = c("id2","patch2"),sep = "\\.")
  distREL_sub <- subset(  distREL,  (sex1=="m" & sex2=="f" |     sex1=="f" & sex2=="m") &    taxon1 == taxon2)
  distREL_sub %>% group_by(taxon1) %>% summarize(   dichrom_colorspan = mean(dS) )
})
#save(coldist_patchesREL,file="coldist_patchesREL.rda")

# again, I'm not always a fancy coder, so I again saved these individually 
# we have 15 tibbles of patches 
patch1_REL <- coldist_patchesREL[[1]]
patch2_REL <- coldist_patchesREL[[2]]
patch3_REL <- coldist_patchesREL[[3]]
patch4_REL <- coldist_patchesREL[[4]]
patch5_REL <- coldist_patchesREL[[5]]
patch6_REL <- coldist_patchesREL[[6]]
patch7_REL <- coldist_patchesREL[[7]]
patch8_REL <- coldist_patchesREL[[8]]
patch9_REL <- coldist_patchesREL[[9]]
patch10_REL <- coldist_patchesREL[[10]]
patch11_REL <- coldist_patchesREL[[11]]
patch12_REL <- coldist_patchesREL[[12]]
patch13_REL <- coldist_patchesREL[[13]]
patch14_REL <- coldist_patchesREL[[14]]
patch15_REL <- coldist_patchesREL[[15]]
# the numbers of taxa in these subsets differ 
nrow(patch1_REL)#94 taxa
nrow(patch2_REL)#97 taxa
nrow(patch3_REL)#97 taxa
nrow(patch4_REL)#97 taxa
nrow(patch5_REL)#97 taxa
nrow(patch6_REL)# 97 taxa  
nrow(patch7_REL)# 96 taxa
nrow(patch8_REL)#97 taxa
nrow(patch9_REL)#97 taxa
nrow(patch10_REL)# 96 taxa  
nrow(patch11_REL)#97 taxa
nrow(patch12_REL)#97 taxa
nrow(patch13_REL)#97 taxa
nrow(patch14_REL)# 96 taxa
nrow(patch15_REL)# 96 taxa
# patch 1 is because two subspecies don't have lores 
# bennetti and lazuli in patch 1 lores only
# Kalbaensis in patches 6,7,13,14
# stictolaemus in patch 10 coverts only 
# i deleted these by hand in the excel dataset. 
# then saved each of these files. 
# uncomment when you want to run 
#write.csv(coldist_patchesREL[[1]],  file="color-data/coldist_patches_distance/coldist_patchesREL/patch1.csv")
#write.csv(coldist_patchesREL[[2]],  file="color-data/coldist_patches_distance/coldist_patchesREL/patch2.csv")
#write.csv(coldist_patchesREL[[3]],  file="color-data/coldist_patches_distance/coldist_patchesREL/patch3.csv")
#write.csv(coldist_patchesREL[[4]],  file="color-data/coldist_patches_distance/coldist_patchesREL/patch4.csv")
#write.csv(coldist_patchesREL[[5]],  file="color-data/coldist_patches_distance/coldist_patchesREL/patch5.csv")
#write.csv(coldist_patchesREL[[6]],  file="color-data/coldist_patches_distance/coldist_patchesREL/patch6.csv")
#write.csv(coldist_patchesREL[[7]],  file="color-data/coldist_patches_distance/coldist_patchesREL/patch7.csv")
#write.csv(coldist_patchesREL[[8]],  file="color-data/coldist_patches_distance/coldist_patchesREL/patch8.csv")
#write.csv(coldist_patchesREL[[9]],  file="color-data/coldist_patches_distance/coldist_patchesREL/patch9.csv")
#write.csv(coldist_patchesREL[[10]],  file="color-data/coldist_patches_distance/coldist_patchesREL/patch10.csv")
#write.csv(coldist_patchesREL[[11]],  file="color-data/coldist_patches_distance/coldist_patchesREL/patch11.csv")
#write.csv(coldist_patchesREL[[12]],  file="color-data/coldist_patches_distance/coldist_patchesREL/patch12.csv")
#write.csv(coldist_patchesREL[[13]],  file="color-data/coldist_patches_distance/coldist_patchesREL/patch13.csv")
#write.csv(coldist_patchesREL[[14]],  file="color-data/coldist_patches_distance/coldist_patchesREL/patch14.csv")
#write.csv(coldist_patchesREL[[15]],  file="color-data/coldist_patches_distance/coldist_patchesREL/patch15.csv")

# then I combined them by hand and made an average of the 15 patches in that 
#excel file. In this dataset for this comparison i've kept kalbaensis. 

# to make it easier for us now, I included that column in our handy dandy data
# file sex_diffs as sex_diffs$AvgDicCol_distance 

# are they correlated? 
cor.test(sex_diffs$AvgDicCol,sex_diffs$AvgDicCol_distance,method = "spearman",exact = FALSE)
# 16 November 2026: I realized the day this study was published online as I 
# write this that I accidentally calculated the correlation of these two metrics
# with a slightly different dataset of specs (44 fewer measurements than the 
# "final" dataset that I have you working from on this version). Ironically, 
# this "final" version has a stronger correlation than the one reported in
# the supplementary figure S15 (a spearman's P of 0.797 versus this which is 
# a spearman's P of 0.803). This is why the numbers aren't exact. Sorry for
# this small mistake! Whoops. We all make mistakes. Let's carry on. 

# ok, yes, they are correlated. Let's plot this. 
figS15 <- ggplot(sex_diffs,aes(x = AvgDicCol, y = AvgDicCol_distance)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "Sexual dichromatism (JND)",y = "Sexual dichromatism (Euclidean dist)")+
  theme_classic()

#pdf(file="comparison_JND_EUC_coldists.pdf", width=5, height=3)
figS15
#dev.off()

# now that we know that it doesn't particularly matter if we use the model 
# of receptor-noise model vs. Euclidian distances, I will continue on with 
# receptor-noise because this model allows me to get a metric for achromatic 
# dichromatism in addition to chromatic. 

# 2. PREPARING THE SYMPATRY DATA
# the sympatry data is estimated by a per-island number of potentially 
# co-occuring Todiramphus taxa. This was a simple count from Figure S1 for
# breeding ranges, which is based on our knowledge of breeding ranges of all
# taxa from species distribution maps, taxonomic records, and data from regional
# bird guides. I compared these distribution maps for non-breeding ranges, 
# for which 40 taxa experience potential differences in sympatry in the 
# non-breeding season. 56 experience no change. Because there is such a 
# difference between taxa that are always isolated from another todiramphus 
# taxon and those that experience some level of sympatry, we transformed 
# our sympatry data. The raw data is already in sex_diffs. 

row.names(sex_diffs) <- sex_diffs$branch
sex_diffs$symp_breed_4thTran<- sex_diffs$symp_breed^(1/4) 
sex_diffs$symp_all_4thTran<- sex_diffs$symp_all^(1/4) 


# 3. CALCULATING ISLAND SIZES
library(sf)
library(sp)
library(terra)
library(PleistoDist)

#load shapefile
# to download the shape file needed for this code, you need to download the 
# Todiramphus_specific.zip file from dryad, unzip it and put it in your
# working directory. In case you have never worked with these types of files
# before, DO NOT delete all the other seemingly superfluous files that do not 
# end with .shp. It will not work if you do that. 

islands <- st_read("Todiramphus_specific.shp")

#define Equi7 Oceania projection
prj_oceania <- st_crs("+proj=aeqd +lat_0=-19.5 +lon_0=131.5 +x_0=6988408.536 +y_0=7654884.537 +datum=WGS84 +units=m +no_defs +type=crs")

#project islands to Equal Earth Asia-Pacific
islands_projected <- st_transform(islands,crs="EPSG:8859")

islands_projected

#load taxa names as list
taxa <- c("alb_alb","alb_orii","alb_owston","albonotatu","aus_aus","aus_dam","aus_inter","aus_odites","aus_tring","chl_azelus","chl_chlori","chl_chloro","chl_collar","chl_davis","chl_occip","chl_terao","cinnamo","colonus","diops","enigma","farquhari","funebris","gertrudae","godeffroyi","lazuli","leucopy","pelewensis","recurviros","reichen","ruficollar","sac_amoe","sac_aneity","sac_brachy","sac_erroma","sac_eximi","sac_juliae","sac_mala","sac_manuae","sac_marin","sac_melano","sac_ornatu","sac_pavuvu","sac_pealei","sac_regina","sac_sacer","sac_santo","sac_solomo","sac_soror","sac_tannen","sac_torres","sac_tucopi","sac_utupua","sac_vicina","sac_viti","sanc_cana","sanc_macmi","sanc_norf","sanc_vagan","sauro_adm","sauro_ana","tris_alber","tris_benne","tris_matt","tris_novae","tris_nusae","tris_stres","tris_tris","tutus_aitu","tutus_mauk","tutus_tutu","ven_venera","ven_young","winch_alf","winch_mind","winch_nesy","winch_nigo","winch_winc")

#create data frame
output <- tibble(.rows = 77)
output$taxon <- taxa
ta <- c()
ma <- c()
cha <- c()
centroids <- data.frame()
ch <- data.frame()

#run for loop that iterates through taxa names, subsets islands, and calculates area
for (i in taxa) {
  #subset relevant features
  islandssubset <- subset(islands_projected,islands_projected[[i]] == "1")
  
  #calculate centroid
  centroids <- rbind(centroids,st_centroid(islandssubset))
  
  #calculate total area
  total_area <- sum(st_area(islandssubset))
  ta <- c(ta,total_area)
  
  #calculate mean area
  mean_area <- mean(st_area(islandssubset))
  ma <- c(ma,mean_area)
  
  #calculate minimum convex hull area
  conv_hull <- st_as_sf(st_convex_hull(st_union(islandssubset)))
  conv_hull$name = i
  ch <- rbind(ch,conv_hull)
  convex_hull_area <- st_area(conv_hull)
  cha <- c(cha,convex_hull_area)
}

output$totalarea <- ta
output$meanarea <- ma
output$convexhull <- cha

#write output file with taxon areas
#write.csv(output,"Todiramphus_areas.csv")

# We have total area and mean area (in meters squared). I log transformed these
# in excel before adding them into the sex_diffs dataset as 
# sex_diffs$log_totalarea and sex_diffs$log_meanarea 

# The heterozygosity data comes directly from McCullough et al. 2025. It is
# in Table S6. It is in our dataset here as sex_diffs$het 

# we will compare whether or not the inclusion of historical toepad samples 
# impacts our results for the heterozygosity dataset. To make this easier, 
# we will drop any samples that are from toepads. Because we are going to plot
# these with a red border later, column sex_diffs$tp is coded like so: 
# toepad sample = "red" 
# tissue sample = "black"

sex_diffs_noTP = sex_diffs
sex_diffs_noTP = sex_diffs_noTP[!sex_diffs_noTP$tp %in% c("red"),]
View(sex_diffs_noTP)

######################################################################
# PART 5: testing adaptive and neutral hypotheses of color evolution #
######################################################################
# We have three general explanations for color rates. 
#1. Competition, as measured by sympatry 
#   a. breeding sympatry (fourth root transformed)
#   b. year round sympatry (fourth root transformed)
#2. Sexual selection, as measured by sexual dichromatism 
#   a. average chromatic dichromatism
#   b. average achromatic dichromatism
#3. Population size, measured by island size and heterozygosity 
#    a. island size total area (log)
#    b. island size mean area (log)
#    c. heterozygosity for all samples (including those from historical toepads)
#    c. heterozygosity for tissue samples only 

# We are going to do linear and phylogenetic linear models for each of these
# comparisons because there are two different opinions about using 
# phylogenetic linear models when dealing with a trait rate (like color rate). 
# Some view a rate as inherently taking phylogenetic history into account, 
# therefore doing a phylogenetic linear model would be taking evolutionary 
# relatedness into account twice. This is the reason why we use traditional 
# linear models. Others view this as superflous and that we still need to take
# the evolutionary relationship into account explicitly with a phylogenetic 
# linear model. This is why we did both linear and phylogenetic linear models. 

# for each of the analyses below, I will have the linear model (corresponding
# to output in table S3), the phylogenetic linear model (table S2), and then
# plot that relationship and name it for which figure it appears in. 

# for the plots, if you ever want to get the points labelled, you can do so with
# the addition of this into the ggplot code block: 
# + geom_text(aes(label = branch), size = 1)

#1. COMPETITION  
# Male chromatic rates and breeding sympatry
m_1a_color_phylolm <- phylolm(colrate.male~symp_breed_4thTran, data=sex_diffs, phy = tre ,model = "lambda")
summary(m_1a_color_phylolm)
m_1a_color_lm <- lm(colrate.male~symp_breed_4thTran, data=sex_diffs)
summary(m_1a_color_lm)

FigS8A <- ggplot(data=sex_diffs, aes(x = symp_breed_4thTran, y = colrate.male)) + geom_point() + 
  geom_point(aes(fill=sp_col), size=4, pch=21, stroke = 0.5)+theme(legend.position = "none")+
  ggtitle("breeding sympatry vs. male col rate")+
  scale_y_continuous(breaks = seq(0, 0.4, by=0.1), limits=c(0,0.46)) + 
  scale_x_continuous(breaks = seq(0, 2.5, by=0.5), limits=c(0,2.7)) 

# Male chromatic rates and non-breeding sympatry
m_1b_color_phylolm <- phylolm(colrate.male~symp_all_4thTran, data=sex_diffs, phy = tre ,model = "lambda")
summary(m_1b_color_phylolm)
m_1b_color_lm <- lm(colrate.male~symp_all_4thTran, data=sex_diffs)
summary(m_1b_color_lm)

FigS8B <- ggplot(data=sex_diffs, aes(x = symp_all_4thTran, y = colrate.male)) + geom_point() + 
  geom_point(aes(fill=sp_col), size=4, pch=21, stroke = 0.53)+theme(legend.position = "none")+
  ggtitle("year round sympatry vs. male col rate")+ 
  scale_y_continuous(breaks = seq(0, 0.4, by=0.1), limits=c(0,0.46)) + 
  scale_x_continuous(breaks = seq(0, 2.5, by=0.5), limits=c(0,2.7)) 

# Male achromatic rates and breeding sympatry
m_1a_lum_phylolm <- phylolm(lumrate.male~symp_breed_4thTran, data=sex_diffs, phy = tre ,model = "lambda")
summary(m_1a_lum_phylolm)
m_1a_lum_lm <- lm(lumrate.male~symp_breed_4thTran, data=sex_diffs)
summary(m_1a_lum_lm)

FigS9A <- ggplot(data=sex_diffs, aes(x = symp_breed_4thTran, y = lumrate.male)) + geom_point() + 
  geom_point(aes(fill=sp_col), size=4, pch=21, stroke = 0.5)+theme(legend.position = "none")+
  ggtitle("breeding sympatry vs. male lum rate")+ geom_smooth(method = "lm", se = FALSE) +
  scale_y_continuous(breaks = seq(0, 0.5, by=0.1), limits=c(0,0.53)) + 
  scale_x_continuous(breaks = seq(0, 2.5, by=0.5), limits=c(0,2.7)) 

# Male achromatic rates and non-breeding sympatry
m_1b_lum_phylolm <- phylolm(lumrate.male~symp_all_4thTran, data=sex_diffs, phy = tre ,model = "lambda")
summary(m_1b_lum_phylolm)
m_1b_lum_lm <- lm(lumrate.male~symp_all_4thTran, data=sex_diffs)
summary(m_1b_lum_lm)

FigS9B <- ggplot(data=sex_diffs, aes(x = symp_all_4thTran, y = lumrate.male)) + geom_point() + 
  geom_point(aes(fill=sp_col), size=4, pch=21, stroke = 0.53)+theme(legend.position = "none")+
  ggtitle("year round sympatry vs. male lum rate")+  geom_smooth(method = "lm", se = FALSE) +
  scale_y_continuous(breaks = seq(0, 0.5, by=0.1), limits=c(0,0.53)) + 
  scale_x_continuous(breaks = seq(0, 2.5, by=0.5), limits=c(0,2.7)) 

# Female chromatic rates and breeding sympatry
f_1a_color_phylolm <- phylolm(colrate.female~symp_breed_4thTran, data=sex_diffs, phy = tre ,model = "lambda")
summary(f_1a_color_phylolm)
f_1a_color_lm <- lm(colrate.female~symp_breed_4thTran, data=sex_diffs)
summary(f_1a_color_lm)

FigS8C <- ggplot(data=sex_diffs, aes(x = symp_breed_4thTran, y = colrate.female)) + geom_point() + 
  geom_point(aes(fill=sp_col), size=4, pch=21, stroke = 0.5)+theme(legend.position = "none")+
  ggtitle("breeding sympatry vs. female col rate")+ 
  scale_y_continuous(breaks = seq(0, 0.4, by=0.1), limits=c(0,0.46)) + 
  scale_x_continuous(breaks = seq(0, 2.5, by=0.5), limits=c(0,2.7)) 

# Female chromatic rates and non-breeding sympatry
f_1b_color_phylolm <- phylolm(colrate.female~symp_all_4thTran, data=sex_diffs, phy = tre ,model = "lambda")
summary(f_1b_color_phylolm)
f_1b_color_lm <- lm(colrate.female~symp_all_4thTran, data=sex_diffs)
summary(f_1b_color_lm)

FigS8D <- ggplot(data=sex_diffs, aes(x = symp_all_4thTran, y = colrate.female)) + geom_point() + 
  geom_point(aes(fill=sp_col), size=4, pch=21, stroke = 0.53)+theme(legend.position = "none")+
  ggtitle("year round sympatry vs. female col rate")+ 
  scale_y_continuous(breaks = seq(0, 0.4, by=0.1), limits=c(0,0.46)) + 
  scale_x_continuous(breaks = seq(0, 2.5, by=0.5), limits=c(0,2.7)) 

# Female achromatic rates and breeding sympatry
f_1a_lum_phylolm <- phylolm(lumrate.female~symp_breed_4thTran, data=sex_diffs, phy = tre ,model = "lambda")
summary(f_1a_lum_phylolm)
f_1a_lum_lm <- lm(lumrate.female~symp_breed_4thTran, data=sex_diffs)
summary(f_1a_lum_lm)

FigS9C <- ggplot(data=sex_diffs, aes(x = symp_breed_4thTran, y = lumrate.female)) + geom_point() + 
  geom_point(aes(fill=sp_col), size=4, pch=21, stroke = 0.5)+theme(legend.position = "none")+
  ggtitle("breeding sympatry vs. female lum rate")+ geom_smooth(method = "lm", se = FALSE) +
  scale_y_continuous(breaks = seq(0, 0.5, by=0.1), limits=c(0,0.53)) + 
  scale_x_continuous(breaks = seq(0, 2.5, by=0.5), limits=c(0,2.7)) 

# Female achromatic rates and non-breeding sympatry
f_1b_lum_phylolm <- phylolm(lumrate.female~symp_all_4thTran, data=sex_diffs, phy = tre ,model = "lambda")
summary(f_1b_lum_phylolm)
f_1b_lum_lm <- lm(lumrate.female~symp_all_4thTran, data=sex_diffs)
summary(f_1b_lum_lm)

FigS9D <- ggplot(data=sex_diffs, aes(x = symp_all_4thTran, y = lumrate.female)) + geom_point() + 
  geom_point(aes(fill=sp_col), size=4, pch=21, stroke = 0.53)+theme(legend.position = "none")+
  ggtitle("year round sympatry vs. female lum rate")+ 
  scale_y_continuous(breaks = seq(0, 0.5, by=0.1), limits=c(0,0.53)) + 
  scale_x_continuous(breaks = seq(0, 2.5, by=0.5), limits=c(0,3)) 

#Making Figure S8: both sexes for color vs. sympatry 
#pdf(file="FigS8_colrates-sympatry.pdf", width=8, height=6)
grid.arrange(FigS8A, FigS8B, FigS8C,FigS8D, ncol=2 )
#dev.off()

#Making Figure S9: both sexes for luminance vs. sympatry 
#pdf(file="FigS9_lumrates-sympatry.pdf", width=8, height=6)
grid.arrange(FigS9A, FigS9B, FigS9C,FigS9D, ncol=2 )
#dev.off()

#2. SEXUAL DICHROMATISM 
# Male chromatic rate and chromatic dichromatism
m_2a_color_phylolm <- phylolm(colrate.male~AvgDicCol, data=sex_diffs, phy = tre ,model = "lambda")
summary(m_2a_color_phylolm)
m_2a_color_lm <- lm(colrate.male~AvgDicCol, data=sex_diffs)
summary(m_2a_color_lm)

FigS7A <- ggplot(data=sex_diffs, aes(x = AvgDicCol , y = colrate.male)) + geom_point() + 
  geom_point(aes(fill=sp_col), size=4, pch=21, stroke = 1)+theme(legend.position = "none")+
  ggtitle("chromatic dichromatism vs. male col rate") +
  scale_y_continuous(breaks = seq(0, 0.5, by=0.1), limits=c(0,0.52)) +
  scale_fill_discrete(name = "branch", labels = setNames(sex_diffs$branch, sex_diffs$sp_col))

# Male achromatic rate and achromatic dichromatism
m_2b_lum_phylolm <- phylolm(lumrate.male~AvgDicLum, data=sex_diffs, phy = tre ,model = "lambda")
summary(m_2b_lum_phylolm)
m_2b_lum_lm <- lm(lumrate.male~AvgDicLum, data=sex_diffs)
summary(m_2b_lum_lm)

FigS7B <- ggplot(data=sex_diffs, aes(x = AvgDicLum , y = lumrate.male)) + geom_point() + 
  geom_point(aes(fill=sp_col), size=4, pch=21, stroke = 1)+theme(legend.position = "none")+
  ggtitle("achromatic dichromatism vs. male lum rate") +
  scale_y_continuous(breaks = seq(0, 0.5, by=0.1), limits=c(0,0.52)) +
  scale_fill_discrete(name = "branch", labels = setNames(sex_diffs$branch, sex_diffs$sp_col))

# Female chromatic rate and chromatic dichromatism
f_2a_color_phylolm <- phylolm(colrate.female~AvgDicCol, data=sex_diffs, phy = tre ,model = "lambda")
summary(f_2a_color_phylolm)
f_2a_color_lm <- lm(colrate.female~AvgDicCol, data=sex_diffs)
summary(f_2a_color_lm)

FigS7C <- ggplot(data=sex_diffs, aes(x = AvgDicCol , y = colrate.female)) + geom_point() + 
  geom_point(aes(fill=sp_col), size=4, pch=21, stroke = 1)+theme(legend.position = "none")+
  ggtitle("chromatic dichromatism vs. female col rate") +
  scale_y_continuous(breaks = seq(0, 0.5, by=0.1), limits=c(0,0.52)) +
  scale_fill_discrete(name = "branch", labels = setNames(sex_diffs$branch, sex_diffs$sp_col))

# Female achromatic rate and achromatic dichromatism
f_2b_lum_phylolm <- phylolm(lumrate.female~AvgDicLum, data=sex_diffs, phy = tre ,model = "lambda")
summary(f_2b_lum_phylolm)
f_2b_lum_lm <- lm(lumrate.female~AvgDicLum, data=sex_diffs)
summary(f_2b_lum_lm)

FigS7D <- ggplot(data=sex_diffs, aes(x = AvgDicLum , y = lumrate.female)) + geom_point() + 
  geom_point(aes(fill=sp_col), size=4, pch=21, stroke = 1)+theme(legend.position = "none")+
  ggtitle("achromatic dichromatism vs. female lum rate") +
  scale_y_continuous(breaks = seq(0, 0.5, by=0.1), limits=c(0,0.52)) +
  scale_fill_discrete(name = "branch", labels = setNames(sex_diffs$branch, sex_diffs$sp_col))

# Making Figure S7 
#pdf(file="FigS7_dichromatism-rates.pdf", width=8, height=6)
grid.arrange(FigS7A,FigS7B, FigS7C, FigS7D, ncol=2)
#dev.off()  

# ISLAND SIZE 
# Male chromatic rate and total island area
m_3a_color_phylolm <- phylolm(colrate.male~log_totalarea, data=sex_diffs, phy = tre ,model = "lambda", lower.bound = 0.0000000001)
summary(m_3a_color_phylolm)
m_3a_color_lm <- lm(colrate.male~log_totalarea, data=sex_diffs)
summary(m_3a_color_lm)

FigS11A <- ggplot(data=sex_diffs, aes(x = log_totalarea, y = colrate.male)) + geom_point() + 
  geom_point(aes(fill=sp_col), size=3, pch=21)+ theme(legend.position = "none")+
  theme(plot.title = element_text(size = 9))

# Male chromatic rate and mean island area
m_3b_color_phylolm <- phylolm(colrate.male~log_meanarea, data=sex_diffs, phy = tre ,model = "lambda", lower.bound = 0.0000000001)
summary(m_3b_color_phylolm)
m_3b_color_lm <- lm(colrate.male~log_meanarea, data=sex_diffs)
summary(m_3b_color_lm)

FigS10A <- ggplot(data=sex_diffs, aes(x = log_meanarea, y = colrate.male)) + geom_point() + 
  geom_point(aes(fill=sp_col), size=3, pch=21)+ theme(legend.position = "none")+
  theme(plot.title = element_text(size = 9))

# Male achromatic rate and total island area
m_3a_lum_phylolm <- phylolm(lumrate.male~log_totalarea, data=sex_diffs, phy = tre ,model = "lambda")
summary(m_3a_lum_phylolm)
m_3a_lum_lm <- lm(lumrate.male~log_totalarea, data=sex_diffs)
summary(m_3a_lum_lm)

FigS11B <- ggplot(data=sex_diffs, aes(x = log_totalarea, y = lumrate.male)) + geom_point() + 
  geom_point(aes(fill=sp_col), size=3, pch=21)+ theme(legend.position = "none")+
  theme(plot.title = element_text(size = 9))

# Male achromatic rate and mean island area
m_3b_lum_phylolm <- phylolm(lumrate.male~log_meanarea, data=sex_diffs, phy = tre ,model = "lambda")
summary(m_3b_lum_phylolm)
m_3b_lum <- lm(lumrate.male~log_meanarea, data=sex_diffs)
summary(m_3b_lum)

FigS10B <- ggplot(data=sex_diffs, aes(x = log_meanarea, y = lumrate.male)) + geom_point() + 
  geom_point(aes(fill=sp_col), size=3, pch=21)+ theme(legend.position = "none")+
  theme(plot.title = element_text(size = 9))

# Female chromatic rate and total island area
f_3a_color_phylolm <- phylolm(colrate.female~log_totalarea, data=sex_diffs, phy = tre ,model = "lambda", lower.bound = 0.0000000001)
summary(f_3a_color_phylolm)
f_3a_color_lm <- lm(colrate.female~log_totalarea, data=sex_diffs)
summary(f_3a_color_lm)

FigS11C <- ggplot(data=sex_diffs, aes(x = log_totalarea, y = colrate.female)) + geom_point() + 
  geom_point(aes(fill=sp_col), size=3, pch=21)+ theme(legend.position = "none")+
  theme(plot.title = element_text(size = 9))

# Female chromatic rate and mean island area
f_3b_color_phylolm <- phylolm(colrate.female~log_meanarea, data=sex_diffs, phy = tre ,model = "lambda", lower.bound = 0.0000000001)
summary(f_3b_color_phylolm)
f_3b_color_lm <- lm(colrate.female~log_meanarea, data=sex_diffs)
summary(f_3b_color_lm)

FigS10C <- ggplot(data=sex_diffs, aes(x = log_meanarea, y = colrate.female)) + geom_point() + 
  geom_point(aes(fill=sp_col), size=3, pch=21)+ theme(legend.position = "none")+
  theme(plot.title = element_text(size = 9))

# Female achromatic rate and total island area
f_3a_lum_phylolm <- phylolm(lumrate.female~log_totalarea, data=sex_diffs, phy = tre ,model = "lambda")
summary(f_3a_lum_phylolm)
f_3a_lum_lm <- lm(lumrate.female~log_totalarea, data=sex_diffs)
summary(f_3a_lum_lm)

FigS11D <- ggplot(data=sex_diffs, aes(x = log_meanarea, y = lumrate.female)) + geom_point() + 
  geom_point(aes(fill=sp_col), size=3, pch=21)+ theme(legend.position = "none")+
  theme(plot.title = element_text(size = 9))

# Female achromatic rate and mean island area
f_3b_lum_phylolm <- phylolm(lumrate.female~log_meanarea, data=sex_diffs, phy = tre ,model = "lambda")
summary(f_3b_lum_phylolm)
f_3b_lum_lm <- lm(lumrate.female~log_meanarea, data=sex_diffs)
summary(f_3b_lum_lm)

FigS10D <- ggplot(data=sex_diffs, aes(x = log_meanarea, y = lumrate.female)) + geom_point() + 
  geom_point(aes(fill=sp_col), size=3, pch=21)+ theme(legend.position = "none")+
  theme(plot.title = element_text(size = 9))

# Making Figure S10
#pdf(file="FigS10_mean_island-size-rates.pdf", width=8, height=6)
grid.arrange(FigS10A,FigS10B, FigS10C, FigS10D, ncol=2)
#dev.off()  

# Making Figure S11
#pdf(file="FigS11_total_island-size-rates.pdf", width=8, height=6)
grid.arrange(FigS11A,FigS11B, FigS11C, FigS11D, ncol=2)
#dev.off()  

# HETEROZYGOSITY 
# Male chromatic rate and heterozygosity (all samples)
m_3c_color_phylolm <- phylolm(colrate.male~het, data=sex_diffs, phy = tre ,model = "lambda", lower.bound = 0.0000000001)
summary(m_3c_color_phylolm)
coef(m_3c_color_phylolm)
m_3c_color_lm <- lm(colrate.male~het, data=sex_diffs)
summary(m_3c_color_lm)

FigS12A <- ggplot(data=sex_diffs, aes(y = colrate.male, x = het)) + geom_point() + 
  geom_point(aes(fill=sp_col), size=4, pch=21, color=sex_diffs$tp, stroke = 1)+
  geom_abline(aes(slope=-21.2287128, intercept= 0.2103221)) +theme(legend.position = "none")+
  ggtitle("Heterozygosity vs. male color rate with toepads (red border)") + 
  scale_y_continuous(breaks = seq(0, 0.5, by=0.1), limits=c(0,0.52))

# Male chromatic rate and heterozygosity (tissue samples only)
m_3d_color_noTP_phylolm <- phylolm(colrate.male~het, data=sex_diffs_noTP, phy = tre ,model = "lambda", lower.bound = 0.0000000001)
summary(m_3d_color_noTP_phylolm)
coef(m_3d_color_noTP_phylolm)
m_3d_pixy2_color_noTP_lm <- lm(colrate.male~het, data=sex_diffs_noTP)
summary(m_3d_color_noTP_lm)
coef(m_3d_color_noTP_lm)

FigS12B <- ggplot(data=sex_diffs_noTP, aes(x = het, y = colrate.male)) + geom_point() + 
  geom_point(aes(fill=sp_col), size=4, pch=21, stroke = 1)+theme(legend.position = "none")+
  geom_abline(aes(slope=-22.7350764, intercept= 0.1924968)) +
  ggtitle("Heterozygosity vs. male color rate with no toepads") + 
  scale_y_continuous(breaks = seq(0, 0.5, by=0.1), limits=c(0,0.52))

# Male achromatic rate and heterozygosity (all samples)
m_3c_lum_phylolm <- phylolm(lumrate.male~het, data=sex_diffs, phy = tre ,model = "lambda")
summary(m_3c_lum_phylolm)
coef(m_3c_lum_phylolm)
m_3c_lum_lm <- lm(lumrate.male~het, data=sex_diffs)
summary(m_3c_lum_lm)
coef(m_3c_lum_lm)

FigS13A <- ggplot(data=sex_diffs, aes(y = lumrate.male, x = het)) + geom_point() + 
  geom_point(aes(fill=sp_col), size=4, pch=21, color=sex_diffs$tp, stroke = 1)+
  geom_abline(aes(slope=-22.6738053, intercept= 0.2377594)) +theme(legend.position = "none")+
  ggtitle("Heterozygosity vs. male lum rate, all (toepads red border)") + 
  scale_y_continuous(breaks = seq(0, 0.5, by=0.1), limits=c(0,0.52))

# Male achromatic rate and heterozygosity (tissue samples only)
m_3d_lum_noTP_phylolm <- phylolm(lumrate.male~het, data=sex_diffs_noTP, phy = tre ,model = "lambda")
summary(m_3d_lum_noTP_phylolm)
coef(m_3d_lum_noTP_phylolm)
m_3d_lum_noTP_lm <- lm(lumrate.male~het, data=sex_diffs_noTP)
summary(m_3d_lum_noTP_lm)
coef(m_3d_lum_noTP_lm)

FigS13B <- ggplot(data=sex_diffs_noTP, aes(x = het, y = lumrate.male)) + geom_point() + 
  geom_point(aes(fill=sp_col), size=4, pch=21, stroke = 1)+theme(legend.position = "none")+
  geom_abline(aes(slope=-26.0390547, intercept= 0.2226315)) +
  ggtitle("Heterozygosity vs. male lum rate with no toepads") + 
  scale_y_continuous(breaks = seq(0, 0.5, by=0.1), limits=c(0,0.52)) 

# Female chromatic rate and heterozygosity (all samples)
f_3c_color_phylolm <- phylolm(colrate.female~het, data=sex_diffs, phy = tre ,model = "lambda", lower.bound = 0.0000000001)
summary(f_3c_color_phylolm)
coef(f_3c_color_phylolm)
f_3c_color_lm <- lm(colrate.female~het, data=sex_diffs)
summary(f_3c_color_lm)

FigS12C <- ggplot(data=sex_diffs, aes(x = het, y = colrate.female)) + geom_point() + 
  geom_point(aes(fill=sp_col), size=4, pch=21, color=sex_diffs$tp, stroke = 1)+
  geom_abline(aes(slope=-16.2644779, intercept= 0.1748061)) +theme(legend.position = "none")+
  ggtitle("Pixy heterozygosity vs. female color rate, all (toepads red border)") + 
  scale_y_continuous(breaks = seq(0, 0.5, by=0.1), limits=c(0,0.52))# + 

# Female chromatic rate and heterozygosity (tissue samples only)
f_3d_color_noTP_phylolm <- phylolm(colrate.female~het, data=sex_diffs_noTP, phy = tre ,model = "lambda", lower.bound = 0.0000000001)
summary(f_3d_color_noTP_phylolm)
coef(f_3d_color_noTP_phylolm)
f_3d_color_noTP_lm <- lm(colrate.female~het, data=sex_diffs_noTP)
summary(f_3d_color_noTP_lm)
coef(f_3d_color_noTP_lm)

FigS12D <- ggplot(data=sex_diffs_noTP, aes(x = het, y = colrate.female)) + geom_point() + 
  geom_point(aes(fill=sp_col), size=4, pch=21, stroke = 1)+
  geom_abline(aes(slope=-16.1692332, intercept= 0.1616226))+ theme(legend.position = "none")+
  ggtitle("Heterozygosity vs. female color rate with no toepads") + 
  scale_y_continuous(breaks = seq(0, 0.5, by=0.1), limits=c(0,0.52))# + 

# Female achromatic rate and heterozygosity (all samples)
f_3c_lum_phylolm <- phylolm(lumrate.female~het, data=sex_diffs, phy = tre ,model = "lambda", lower.bound = 0.0000000001)
summary(f_3c_lum_phylolm)
coef(f_3c_lum_phylolm)
f_3c_lum_lm <- lm(lumrate.female~het, data=sex_diffs)
summary(f_3c_lum_lm)
coef(f_3c_lum_lm)

FigS13C <- ggplot(data=sex_diffs, aes(x = het, y = lumrate.female)) + geom_point() + 
  geom_point(aes(fill=sp_col), size=4, pch=21, color=sex_diffs$tp, stroke = 1)+theme(legend.position = "none")+
  geom_abline(aes(slope=-18.594426, intercept= 0.230352)) +
  ggtitle("Pixy heterozygosity vs. female lum rate") + 
  scale_y_continuous(breaks = seq(0, 0.5, by=0.1), limits=c(0,0.52))  


# Female achromatic rate and heterozygosity (tissue samples only)
f_3c_lum_noTP_phylolm <- phylolm(lumrate.female~het, data=sex_diffs_noTP, phy = tre ,model = "lambda", lower.bound = 0.0000000001)
summary(f_3c_lum_noTP_phylolm)
coef(f_3c_lum_noTP_phylolm)
f_3c_lum_noTP_lm <- lm(lumrate.female~het, data=sex_diffs_noTP)
summary(f_3c_lum_noTP_lm)
coef(f_3c_lum_noTP_lm)

FigS13D <- ggplot(data=sex_diffs_noTP, aes(x = het, y = lumrate.female)) + geom_point() + 
  geom_point(aes(fill=sp_col), size=4, pch=21, stroke = 1)+theme(legend.position = "none")+
  geom_abline(aes(slope=-27.8462392, intercept= 0.2301356)) +
  ggtitle("Heterozygosity vs. female lum rate with no toepads") + 
  scale_y_continuous(breaks = seq(0, 0.5, by=0.1), limits=c(0,0.52)) 

#Making Figure S12
#pdf(file="FigS12_chromatic-rates-heterozygosity.pdf", width=8, height=6)
grid.arrange(FigS12A,FigS12B, FigS12C, FigS12D, ncol=2)
#dev.off()  

#Making Figure S13 
#pdf(file="FigS13_achromatic-rates-heterozygosity.pdf", width=8, height=6)
grid.arrange(FigS13A,FigS13B, FigS13C, FigS13D, ncol=2)
#dev.off()  

### Multiple variable linear models 

# For achromatic rates, the linear models showed sympatry has also having some
# predictive power. Now we will compare them. This is related to Table S6. 

#male luminance rates for linear models only 
# heterozygosity  
summary(m_3c_lum_lm)
#sympatry in breeding season 
summary(m_1a_lum_lm)
#sympatry year-round
summary(m_1b_lum_lm)

lum_breedingsympatry_het_m <- lm(lumrate.male~symp_breed_4thTran+het, data = sex_diffs)
summary(lum_breedingsympatry_het_m)

lum_allsympatry_het_m <- lm(lumrate.male~het+symp_all_4thTran, data = sex_diffs)
summary(lum_allsympatry_het_m)

AIC(m_3c_lum_lm,m_1a_lum_lm,m_1b_lum_lm,lum_breedingsympatry_het_m, lum_allsympatry_het_m)

# female luminance rates
# heterozygosity  
summary(f_3c_lum_lm)
#explained by sympatry in breeding season 
summary(f_1a_lum_lm)

# multiple linear regressions 
# heterozygosity and sympatry 
lum_sympatry_het_f <- lm(lumrate.female~het+symp_breed_4thTran, data = sex_diffs)
summary(lum_sympatry_het_f)

AIC(f_3c_lum_lm, f_1a_lum_lm,lum_sympatry_het_f)

################################################
# PART 6: Sensitivity tests for heterozygosity #
################################################ 

# There is a known bias in the estimation of phenotypic rates for short 
# branches in phylogenetic trees. The crux is that many trait evolutionary
# models tend to produce inflated (i.e. artificially fast) rate estimates for
# short branches in the tree, which can cause problems for testing downstream
# relationships with variables that are also calculate from (or related) to 
# phylogenetic branch length.This was illustrated by Cooney and Thomas 2021
# Nat. Ecol. & Evol. 

# Much of the results relating to this section is dicussed in the additional
# results document that is available on dryad/Zenodo. 

# Now, we want to double-check our main result that only heterozygosity.
# We have three sensitivity tests that we are going to do in this section. 
# our questions are: 

#1. Are heterozygosity estimates are uncorrelated with taxon age (i.e. terminal
# branch length in the tree)? 

#2. Do relationships persist when using colour rates derived from lambda-scaled
# trees? 

#3. Do relationships persist when including taxon age (i.e. terminal branch 
# length in the tree) as an additional predictor in the multipredictor models
# to account for the fact that short branches are invariably inferred to have
# faster rates?  

#1. RELATIONSHIP OF HETEROZYGOSITY AND TERMINAL BRANCH LENGTH 

# using our sexs_diffs dataset, I've already added the terminal branchs lengths
# in the dataset from the timetree in the column sex_diffs$divergence 
cor.test(sex_diffs$het,sex_diffs$divergence,method = "pearson")

#lets plot this relationship  
FigS17 <- ggplot(data=sex_diffs, aes(x = het, y =divergence )) +
  geom_smooth(method = "lm", se = TRUE) + 
  geom_point(aes(fill=islands), size=4, pch=21, stroke = 1,color=sex_diffs$tp) + 
  geom_point(aes(fill=islands), size=4, pch=21, stroke = 1)+theme(legend.position = "none")+
  ggtitle("heterozygosity vs. divergence") 

#pdf(file="FigS17_heterozygosity-divergence.pdf",width=6.5, height=4)
FigS17
#dev.off()

## Lets do this for the tissue-only dataset as well. 
cor.test(sex_diffs_noTP$het,sex_diffs_noTP$divergence,method = "pearson")

#is pyrrhopygius causing this issue? 
sex_diffs_noTP[which.max(sex_diffs_noTP$divergence),c("branch", "het", "divergence")]
no_pyrr <- subset(sex_diffs_noTP, branch != "Todiramphus_pyrrhopygius")

cor.test(no_pyrr$het,no_pyrr$divergence,method = "pearson")

# These results are reflected in Table S12: 
#With T. pyrrhopygius included: 
#r = 0.519
#p = 0.000263
#n = 45

#Without it:
#r = 0.156
#p = 0.311
#n = 44
#95% CI = −0.147 to 0.433

#does that change the whole dataset? 
no_pyrr_full <- subset(sex_diffs,branch != "Todiramphus_pyrrhopygius")

cor.test(no_pyrr_full$het,no_pyrr_full$divergence,method = "pearson")

cor.test(sex_diffs$het,sex_diffs$divergence,method = "spearman",exact = FALSE)

cor.test(sex_diffs_noTP$het,sex_diffs_noTP$divergence,method = "spearman",exact = FALSE)

# so, yes, these two aspects are correlated. But lets see if rescaling the
# phylogeny based on lambda will wholesale change our results. 

# 2. COLOR EVOLUTIONARY RATES ESTIMATED FROM λ-SCALED TREES. 
library(mvMORPH)

mal96_14_xyz <- read.csv("color-data/tax96-m-chrom-avg-specs3.xyz.nolores.csv")
fem95_14_xyz <- read.csv("color-data/tax95-f-chrom-avg-specs3.xyz.nolores.csv")
fem95_14_lum <- read.csv("color-data/tax95-f-achrom-avg-specs3-nolores.lum.csv")
mal96_14_lum <- read.csv("color-data/tax96-m-achrom-avg-specs3-nolores.lum.csv")

#the code below for RRphylo will overwrite each other so make sure you 
# are running the lines from the top of each section (i.e., do all of 
# male chromatic before moving on to male achromatic because it will overwrite). 
# I'm a lazy coder sometimes. 

setdiff(tre$tip.label, rownames(mal96_14_xyz))
setdiff(rownames(mal96_14_xyz), tre$tip.label)
mal96_14_xyz <- mal96_14_xyz[tre$tip.label, ]
identical(rownames(mal96_14_xyz), tre$tip.label)
xyz <- as.matrix(mal96_14_xyz)

#we are going to estimate lambda for these traits. 
fit_lambda <- mvgls(
  xyz ~ 1,
  tree = Tre,
  model = "lambda",
  method = "LOO",
  penalty = "RidgeArch"
)

fit_lambda$param
# lambda 0.9913

tre_lambda <- fit_lambda$corrSt$phy
class(tre_lambda)
head(tre$edge.length)
head(tre_lambda$edge.length)
plot(tre_lambda)
lambda <- fit_lambda$param
tre_lambda <- fit_lambda$corrSt$phy

## Ok now we can redo RRphylo, this again is for males and XYZ coordinates 
RR <- RRphylo(tree=tre_lambda, y=mal96_14_xyz)
pRR <- plotRR(RR,y=mal96_14_xyz)
id <- c((Ntip(tre_lambda)+1):(Ntip(tre_lambda)+Nnode(tre)), 1:Ntip(tre_lambda))
erates <- RR$rates[,1][match(tre_lambda$edge[,2], id)]
res <- search.shift(RR, status.type="clade")
res
pRR$plotRRrates(tree.args=list(edge.width=4),color.pal=viridis,
                colorbar.args = list(x="bottomleft",direction="h",labs.adj=0.7,xpd=TRUE))
#write.csv(erates, "14patch.RESCALED.m.xyzcoord.erates.csv")

# Because there is no aneityum in the female dataset, we need to drop it from
# the tree before we estimate lambda. 
tre_f95 <- drop.tip(tre,"Todiramphus_sacer_aneityum", rooted = is.rooted(tre))
rownames(fem95_14_xyz) <- fem95_14_xyz$spp
fem95_14_xyz = subset(fem95_14_xyz, select = -c(spp) )

xyz <- as.matrix(fem95_14_xyz)
fit_lambda <- mvgls(
  xyz ~ 1,
  tree = tre_f95,
  model = "lambda",
  method = "LOO",
  penalty = "RidgeArch"
)
summary(fit_lambda)

fit_lambda$param
# lambda 0.997751

tre_f95_lambda <- fit_lambda$corrSt$phy
lambda_f95 <- fit_lambda$param
coltre_f95_lambda <- fit_lambda$corrSt$phy
str(coltre_f95_lambda)

#redo RRphylo
RR <- RRphylo(tree=coltre_f95_lambda, y=fem95_14_xyz)
pRR <- plotRR(RR,y=fem95_14_xyz)
id <- c((Ntip(coltre_f95_lambda)+1):(Ntip(coltre_f95_lambda)+Nnode(coltre_f95_lambda)), 1:Ntip(coltre_f95_lambda))

erates <- RR$rates[,1][match(coltre_f95_lambda$edge[,2], id)]
res <- search.shift(RR, status.type="clade")
res
pRR$plotRRrates(tree.args=list(edge.width=4),color.pal=viridis,
                colorbar.args = list(x="bottomleft",direction="h",labs.adj=0.7,xpd=TRUE))
#write.csv(erates, "14patch.RESCALED.f.xyzcoord.erates.csv")

# Now let's do the rescaled tree 
# remember I dropped a tip from the timetree like this. 
#tre_f95 <- drop.tip(tre,"Todiramphus_sacer_aneityum", rooted = is.rooted(tre))
rownames(fem95_14_lum) <- fem95_14_lum$spp
fem95_14_lum = subset(fem95_14_lum, select = -c(spp) )
xyz <- as.matrix(fem95_14_lum)
fit_lambda <- mvgls(
  xyz ~ 1,
  tree = tre_f95,
  model = "lambda",
  method = "LOO",
  penalty = "RidgeArch"
)
summary(fit_lambda)

fit_lambda$param
# lambda 0.7852

lumtre_f95_lambda <- fit_lambda$corrSt$phy
str(lumtre_f95_lambda)
fem95_14_xyz <- fem95_14_xyz[lumtre_f95_lambda$tip.label, ]

RR <- RRphylo(tree=lumtre_f95_lambda, y=fem95_14_lum)
pRR <- plotRR(RR,y=fem95_14_lum)
id <- c((Ntip(lumtre_f95_lambda)+1):(Ntip(lumtre_f95_lambda)+Nnode(lumtre_f95_lambda)), 1:Ntip(lumtre_f95_lambda))
erates <- RR$rates[,1][match(lumtre_f95_lambda$edge[,2], id)]
res <- search.shift(RR, status.type="clade")
res
pRR$plotRRrates(tree.args=list(edge.width=4),color.pal=viridis,
                colorbar.args = list(x="bottomleft",direction="h",labs.adj=0.7,xpd=TRUE))
#write.csv(erates, "14patch.RESCALED.f.luminance.erates.csv")

#males
rownames(mal96_14_lum) <- mal96_14_lum$spp
mal96_14_lum = subset(mal96_14_lum, select = -c(spp) )

RR <- RRphylo(tree=tre, y=mal96_14_lum)
pRR <- plotRR(RR,y=mal96_14_lum)
id <- c((Ntip(tre)+1):(Ntip(tre)+Nnode(tre)), 1:Ntip(tre))
erates <- RR$rates[,1][match(tre$edge[,2], id)]
res <- search.shift(RR, status.type="clade")
res
pRR$plotRRrates(tree.args=list(edge.width=4),color.pal=viridis,
                colorbar.args = list(x="bottomleft",direction="h",labs.adj=0.7,xpd=TRUE))
#write.csv(erates, "14patch.m.luminance.erates.csv")

xyz <- as.matrix(mal96_14_lum)
fit_lambda <- mvgls(
  xyz ~ 1,
  tree = tre,
  model = "lambda",
  method = "LOO",
  penalty = "RidgeArch"
)
summary(fit_lambda)

fit_lambda$param
# lambda 0.9287354

lumtre_m_lambda <- fit_lambda$corrSt$phy
str(lumtre_m_lambda)
mal96_14_lum <- mal96_14_lum[lumtre_m_lambda$tip.label, ]

RR <- RRphylo(tree=lumtre_m_lambda, y=mal96_14_lum)
pRR <- plotRR(RR,y=mal96_14_lum)
id <- c((Ntip(lumtre_m_lambda)+1):(Ntip(lumtre_m_lambda)+Nnode(lumtre_m_lambda)), 1:Ntip(lumtre_m_lambda))
erates <- RR$rates[,1][match(lumtre_m_lambda$edge[,2], id)]
res <- search.shift(RR, status.type="clade")
res
pRR$plotRRrates(tree.args=list(edge.width=4),color.pal=viridis,
                colorbar.args = list(x="bottomleft",direction="h",labs.adj=0.7,xpd=TRUE))
#write.csv(erates, "14patch.RESCALED.m.luminance.erates.csv")

# All of those outputs are already added to sex_diffs as four different columns. 
str(sex_diffs)
# they have "rescaled.XXXX" format. 

#So, does rescaling the tree cause the signal to disappear? 
# This is recreating Table S6. 

lm_rescaled_col_m_het <- lm(rescaled.colrate.male ~ het,data = sex_diffs)
summary(lm_rescaled_col_m_het)
coef(lm_rescaled_col_m_het)

lm_rescaled_col_m_het_noTP <- lm(rescaled.colrate.male ~ het,data = sex_diffs_noTP)
summary(lm_rescaled_col_m_het_noTP)
coef(lm_rescaled_col_m_het_noTP)

lm_rescaled_col_f_het <- lm(rescaled.colrate.female ~ het,data = sex_diffs)
summary(lm_rescaled_col_f_het)
coef(lm_rescaled_col_f_het)

lm_rescaled_col_f_het_noTP <- lm(rescaled.colrate.female ~ het,data = sex_diffs_noTP)
summary(lm_rescaled_col_f_het_noTP)
coef(lm_rescaled_col_f_het_noTP)

lm_rescaled_lum_m_het <- lm(rescaled.lumrate.male ~ het,data = sex_diffs)
summary(lm_rescaled_lum_m_het)
coef(lm_rescaled_lum_m_het)

lm_rescaled_lum_m_het_noTP <- lm(rescaled.lumrate.male ~ het,data = sex_diffs_noTP)
summary(lm_rescaled_lum_m_het_noTP)
coef(lm_rescaled_lum_m_het_noTP)

lm_rescaled_lum_f_het <- lm(rescaled.lumrate.female ~ het, data = sex_diffs)
summary(lm_rescaled_lum_f_het)
coef(lm_rescaled_lum_f_het)

lm_rescaled_lum_f_het_noTP <- lm(rescaled.lumrate.female ~ het, data = sex_diffs_noTP)
summary(lm_rescaled_lum_f_het_noTP)
coef(lm_rescaled_lum_f_het_noTP)

#These analyses recovered the same general relationship as our original 
# analyses: lower heterozygosity was associated with faster color evolutionary
# rates. Female achromatic evolution had the lowest phylogenetic signal (λ = 
# 0.79, compared with λ = 0.92–0.99 for the other three datasets). Importantly,
# heterozygosity significantly predicted female achromatic evolutionary rates in
#the λ-scaled analyses for both the complete and modern-tissue-only datasets 
# (p=0.005–0.017).

#3. MODELS ACCOUNTING DIRECTLY FOR TERMINAL BRANCH LENGTH
#Does heterozygosity predict color evolutionary rate after controlling for 
# taxon age?

# We additionally repeated the heterozygosity models from the original 
# non-λ-scaled trees with terminal branch length included as a second predictor.
# this is recreating table S5.  

#male chromatic rate and heterozygosity for all samples
m_col_het <- lm(colrate.male~het, data = sex_diffs)
summary(m_col_het)

#male chromatic rate, heterozygosity, and divergence time for all samples
m_col_het_age <- lm(colrate.male~het+divergence, data = sex_diffs)
summary(m_col_het_age)

AIC(m_col_het, m_col_het_age)

#male chromatic rate, heterozygosity, and divergence time for just tissue-based
m_col_het_noTP <- lm(colrate.male~het, data = sex_diffs_noTP)
summary(m_col_het_noTP)

#male chromatic rate and heterozygosity for just tissue-based samples 
m_col_het_age_noTP <- lm(colrate.male~het+divergence, data = sex_diffs_noTP)
summary(m_col_het_age_noTP)

AIC(m_col_het_noTP, m_col_het_age_noTP)

#male achromatic rate, heterozygosity, and divergence time for all samples
m_lum_het_age <- lm(lumrate.male~het+divergence, data = sex_diffs)
summary(m_lum_het_age)

#male achromatic rate and heterozygosity for all samples
m_lum_het <- lm(lumrate.male~het, data = sex_diffs)
summary(m_lum_het)

AIC(m_lum_het, m_lum_het_age)


#male achromatic rate, heterozygosity, and divergence time for just tissue-based
m_lum_het_noTP <- lm(lumrate.male~het, data = sex_diffs_noTP)
summary(m_lum_het_noTP)

#male achromatic rate and heterozygosity for just tissue-based samples 
m_lum_het_age_noTP <- lm(lumrate.male~het+divergence, data = sex_diffs_noTP)
summary(m_lum_het_age_noTP)

AIC(m_lum_het_noTP,m_lum_het_age_noTP)

#female chromatic rate and heterozygosity for all samples
f_col_het <- lm(colrate.female~het, data = sex_diffs)
summary(f_col_het)

#female chromatic rate, heterozygosity, and divergence time for all samples
f_col_het_age <- lm(colrate.female~het+divergence, data = sex_diffs)
summary(f_col_het_age)

AIC(f_col_het, f_col_het_age)


#female chromatic rate, heterozygosity, and divergence time for just 
# tissue-based
f_col_het_noTP <- lm(colrate.female~het, data = sex_diffs_noTP)
summary(f_col_het_noTP)

#female chromatic rate and heterozygosity for just tissue-based samples 
f_col_het_age_noTP <- lm(colrate.female~het+divergence, data = sex_diffs_noTP)
summary(f_col_het_age_noTP)

AIC(f_col_het_age_noTP,f_col_het_noTP)

#female achromatic rate and heterozygosity for all samples
f_lum_het <- lm(lumrate.female~het, data = sex_diffs)
summary(f_lum_het)

#female achromatic rate, heterozygosity, and divergence time for all samples
f_lum_het_age <- lm(lumrate.female~het+divergence, data = sex_diffs)
summary(f_lum_het_age)

AIC(f_lum_het, f_lum_het_age)


#female achromatic rate, heterozygosity, and divergence time for just
# tissue-based
f_lum_het_noTP <- lm(lumrate.female~het, data = sex_diffs_noTP)
summary(f_lum_het_noTP)
#female achromatic rate and heterozygosity for just tissue-based samples 
f_lum_het_age_noTP <- lm(lumrate.female~het+divergence, data = sex_diffs_noTP)
summary(f_lum_het_age_noTP)

AIC(f_lum_het_noTP,f_lum_het_age_noTP)

anova(m_lum_het, m_lum_het_age)
anova(f_col_het, f_col_het_age)
anova(f_lum_het, f_lum_het_age)

# When all taxa were included, terminal branch length was significantly 
# associated with evolutionary rate in all four analyses (p = 0.0036–0.0178), 
# with shorter branches associated with higher inferred evolutionary rates, 
# consistent with the potential short-branch effect. Including terminal branch
# length also improved model support relative to heterozygosity-only models 
# (ΔAIC = 3.83–6.81). Importantly, however, heterozygosity remained a 
# significant negative predictor of male chromatic (p = 0.0141), female
# chromatic (p = 0.0076), and male achromatic evolutionary rates (p = 0.0417)
# after accounting for terminal branch length. Only for female achromatic 
# evolution was heterozygosity no longer significant (p = 0.175). But when 
# analyses were restricted to heterozygosity estimates derived from high-quality
#tissue samples, terminal branch length was not significant in any of the four
# models (p = 0.274–0.721). Models with and without terminal branch length 
# received similar AIC support (ΔAIC = 0.70–1.86), with the heterozygosity-only
# model having the lower AIC in all four comparisons. Heterozygosity was 
# significantly associated with all four evolutionary-rate measures in the 
# single-predictor models (p = 0.0037–0.0082). After terminal branch length was
# included, heterozygosity remained significant for male chromatic evolution 
# (p = 0.0201) and approached significance for male achromatic, female 
# chromatic, and female achromatic evolution (p = 0.051–0.065) in this 
# substantially reduced dataset.

