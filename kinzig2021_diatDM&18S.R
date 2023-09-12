# SFB RESIST
# Kinzig 2021
# data analysed analysed here are from OMNIDIA output
# 14.08.2023
# Serge Mayombo
# digital microscopy and diatom subset of MPB biofilm 18SV9 amplicon sequencing
# data analysis

# load libraries
library(tidyverse)
#install.packages("vegan")
library(vegan)

# read environmental dataset

data <- field <- read.csv("kinzig2021_fieldData.csv", sep = ",",
                          header = T)

# discard FAL2 from env. data
data <- data[-18,]

# change class of surrounding to factor variable with 2 levels
data$surrounding <- as.factor(data$surrounding) 
str(data)
######
# load data on R
diat <- read.csv("kinzig_paper.csv", sep = ",",
                 header = T)
# select only morphology (DM) data
diat_dm <- diat %>% filter(method == "MIC")
# select only 18SV9 data
diat_18S <- diat %>% filter(method == "DNA")
# split the site_name column in 2 columns to separate site_name and replicate
# We use separate_wider_delim()
#diat_18SV9 <- diat_18S %>% separate_wider_delim(site_name, "_", names = c("site_name", "replicate"))

# Taxa plot
# working with microscopy dataset
# select all columns except "method"
spe <- diat_dm %>% dplyr::select(-method)

# make the dataset tidy

spe1 <- spe %>% 
  pivot_longer(!site_name, names_to = "taxa", values_to = "proportion")

# Keep only taxa with proportion >= 5, or 5% of the total relative abundance per sample.
spe2 <- subset(spe1, proportion >= 3)
# Now we can spread the dataset of the most dominant taxa with at least 5% of proportion per sample

spe3 <- spe2 %>%
  pivot_wider(names_from = taxa, values_from = proportion, values_fn = sum, values_fill = 0)

# make site_name column as rownames
spe4 <- spe3 %>% column_to_rownames(var="site_name") # change row names

# plotting relative abundances of dominant taxa

########
# taxa plot
########
set.seed(123)
x<-spe4/rowSums(spe4)
x<-x[,order(colSums(x),decreasing=TRUE)]
head(x)

#Extract list of top N Taxa
N<-32
taxa_list<-spe4
N<-length(taxa_list)

#Create a custom color scale
library(RColorBrewer)
# Define the number of colors you want
nb.cols <- 32
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)

# other colours
library(RColorBrewer)
#install.packages("randomcoloR")
library(randomcoloR)
n <- 32
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))
palette <- distinctColorPalette(n)

# Make a new dataframe
df<-NULL
for (i in 1:dim(x)[2]){
  tmp<-data.frame(row.names=NULL,Sample=rownames(x), surrounding=data$surrounding,  
                  Taxa=rep(colnames(x)[i],dim(x)[1]),Value=x[,i])
  if(i==1){df<-tmp} else {df<-rbind(df,tmp)}
}
# Taxa plot based on surrounding 

kin2021_dominant_DM <- ggplot(df, aes(x=Sample,y=Value, fill=Taxa))+geom_bar(stat="identity")+
  facet_grid(. ~ surrounding, drop=TRUE,scale="free",space="free_x")+
  scale_fill_manual(values=palette[1:(N+1)]) + theme_bw()+ylab("Relative abundance") +
  scale_y_continuous(expand = c(0,0))+theme(strip.background = element_rect(fill="gray85"))+
  theme(panel.margin = unit(0.3, "lines")) +
  theme(axis.text.x=element_text(angle=30,hjust=1,vjust=0.5, size = 12, face = "bold"), 
        strip.text = element_text(size = 14, face = "bold"))+
  theme(legend.key.size = unit(0.6, "cm"), legend.key.width = unit(0.6,"cm"), 
        legend.text = element_text(size = 11))+
  guides(fill=guide_legend(ncol=1))

kin2021_dominant_DM
# save the plot
ggsave("kin2021_dominant_DM.tiff", units="in", width=9, height=9, dpi=300, compression = 'lzw')

##############################
# working with amplicon data
# select all columns except "method"
dna <- diat_18S %>% dplyr::select(-method)

# make the dataset tidy

dna1 <- dna %>% 
  pivot_longer(!site_name, names_to = "taxa", values_to = "proportion")

# Keep only taxa with proportion >= 3, or 3% of the total relative abundance per sample.
dna2 <- subset(dna1, proportion >= 3)
# Now we can spread the dataset of the most dominant taxa with at least 5% of proportion per sample

dna3 <- dna2 %>%
  pivot_wider(names_from = taxa, values_from = proportion, values_fn = sum, values_fill = 0)


# make site_name column as rownames
dna4 <- dna3 %>% column_to_rownames(var="site_name") # change row names

dna4
# plotting relative abundances of dominant taxa based on reads

########
# taxa plot
########
set.seed(123)
x1<-dna4/rowSums(dna4)
x1<-x1[,order(colSums(x1),decreasing=TRUE)]
head(x1)
#Extract new list of top N Taxa
N<-20
#taxa_list<-dna4
taxa_list<-colnames(x1)[1:N]
N<-length(taxa_list)
#N1<-length(taxa_list1)

#Create a custom color scale
library(RColorBrewer)
# Define the number of colors you want
nb.cols <- 20
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)

# other colours
library(RColorBrewer)
#install.packages("randomcoloR")
library(randomcoloR)
set.seed(123)
n <- 20
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))
palette <- distinctColorPalette(n)
#Colours <- brewer.pal(50,"Set")
#names(Colours) <- levels(diatom_ab$Taxa)
#colScale <- scale_colour_manual(name = "Taxa",values = Colours)

#Generate a new table with everything added to Others
#new_x<-data.frame(x[,colnames(x) %in% taxa_list],Others=rowSums(x[,!colnames(x) %in% taxa_list]))
#You can change the Type=grouping_info[,1] should you desire any other grouping of panels

df1<-NULL
for (i in 1:dim(x1)[2]){
  tmp<-data.frame(row.names=NULL,Sample=rownames(x1), surrounding=data$surrounding,  
                  Taxa=rep(colnames(x1)[i],dim(x1)[1]),Value=x1[,i])
  if(i==1){df1<-tmp} else {df1<-rbind(df1,tmp)}
}

# Based on site location along the river system (upstream, midstream or downstream)


kin2021_dominant_dna <- ggplot(df1, aes(x=Sample,y=Value, fill=Taxa))+geom_bar(stat="identity")+
  facet_grid(. ~ surrounding, drop=TRUE,scale="free",space="free_x")+
  scale_fill_manual(values=palette[1:(N+1)]) + theme_bw()+ylab("Relative abundance") +
  scale_y_continuous(expand = c(0,0))+theme(strip.background = element_rect(fill="gray85"))+
  theme(panel.margin = unit(0.3, "lines")) +
  theme(axis.text.x=element_text(angle=30,hjust=1,vjust=0.5, size = 12, face = "bold"), 
        strip.text = element_text(size = 14, face = "bold"))+
  theme(legend.key.size = unit(0.7, "cm"), legend.key.width = unit(0.7,"cm"), 
        legend.text = element_text(size = 12))+
  guides(fill=guide_legend(ncol=1))

kin2021_dominant_dna
# save the plot
ggsave("kin2021_dominant_dna.tiff", units="in", width=9, height=9, dpi=300, compression = 'lzw')

################################
# nMDS --------------------------------------------------------------------
# digital microscopy data
# microscopy
spp.log_dm <- decostand(spe4, method = "log")
# presence/absence
spp.pa_dm <- decostand(spe4, method = "pa")

#spp.log1 <- decostand(spe5, method = "log")
spp.log.dis_dm <- vegdist(spp.log_dm, method = "bray")
#spp.log.dis1 <- vegdist(spp.log, method = "bray")
spp.log_dm
spp.log.dis_dm
# check
decorana (spp.log_dm) # The length of the first DCA axis is 2.7 S.D. (i.e. < 3. S.D.), 
# and these data are thus suitable for linear ordination methods.
#library(vegan)

# Logarithmic transformation as suggested by Anderson et al. (2006): log_b (x) + 1 for x > 0, where b is the base of 
# the logarithm; zeros are left as zeros. 



# betadisper --------------------------------------------------------------

# Before doing the PERMANOVA, first we check to see if the dispersion is the same
# Homogeneity of groups
# betadisper studies the differences in group homogeneities
# analogous to Levene's test of the equality of variances
# can only use one factor as an independent variable

# surrounding
(mod.surrounding <- with(data, betadisper(spp.log.dis_dm, data$surrounding)))
plot(mod.surrounding, sub = NULL) 
boxplot(mod.surrounding)
anova(mod.surrounding)   # this says that within-group variances can be considered homogeneous - in fact a good news!
# next step could be to compare the location of centroids and see if there is
# significant difference between groups in that respect!
permutest(mod.surrounding)
# PERMANOVA ---------------------------------------------------------------

# Permutational multivariate analysis of variance using distance matrices
# (Bray-Curtis similarities by default). ANOSIM uses only ranks of Bray-Curtis,
# so the former preserves more information.

(perm.1 <- adonis2(spp.log.dis_dm~(surrounding*pH*water_temp.*oxygen_full*dissolved_oxygen*oxygen_saturation*conductivity*average_velocity*average_depth*discharge*water_level_at_pole),
                   method = perm, data = data))

(perm.1 <- adonis2(spp.log.dis_dm~(surrounding+pH+water_temp.+oxygen_full+dissolved_oxygen+oxygen_saturation+conductivity+average_velocity+average_depth+discharge+water_level_at_pole),
                   method = perm, data = data))
## great - this tells us that salinity does have a significant effect
# but not velocity
# so the simplest model we can take would be:
(perm.2 <- adonis(spp.log.dis_dm ~ surrounding,
                  method = perm, data = data))

(perm.2 <- adonis2(spp.log.dis_dm ~ conductivity,
                  method = perm, data = data))
##################################################
# NMDS
set.seed(123)
spp.nmds <- metaMDS(spp.log_dm, k = 2,trymax = 100, permutation = 9999,
                    distance = "bray", wascores = TRUE)

spp.nmds
stressplot(spp.nmds)

##########################################
#############################
#spp.mds <- metaMDS((spe_count), distance = "bray", autotransform = F)
#spp.mds
env.fit <- envfit(spp.nmds, data, permutations = 999, na.rm = TRUE) # this fits environmental vectors
spp.fit <- envfit(spp.nmds, spe4, permutations = 999) # this fits species vectors

env.fit
spp.fit
# To plot the output from the mds using ggplot a new datasheet needs to be created which contains the x,y points for each site. 
#You can do this by calling the scores of you mds.

sample.scrs <- as.data.frame(scores(spp.nmds, display = "sites")) #save NMDS results into dataframe
sample.scrs <- cbind(sample.scrs, surrounding = data$surrounding) #add grouping variable "Management" to dataframe
#sample.scrs <- cbind(sample.scrs, year.of.restor. = metadata1$year.of.r) #add grouping variable of cluster grouping to dataframe
#site.scrs <- cbind(site.scrs, Site = rownames(site.scrs)) #add site names as variable if you want to display on plot
#sample.scrs <- cbind(sample.scrs, location = env$site)

head(sample.scrs)
#view(sample.scrs)

# A new dataset containing species data also needs to be made to look at species vectors.
#This is not necessary if you don't want to show the species on the final graph.
spp.scrs <- as.data.frame(scores(spp.fit, display = "vectors")) #save species intrinsic values into dataframe
spp.scrs <- cbind(spp.scrs, Species = rownames(spp.scrs)) #add species names to dataframe
spp.scrs <- cbind(spp.scrs, pval = spp.fit$vectors$pvals) #add pvalues to dataframe so you can select species which are significant
#spp.scrs<- cbind(spp.scrs, abrev = abbreviate(spp.scrs$Species, minlength = 6)) #abbreviate species names
sig.spp.scrs <- subset(spp.scrs, pval<=0.05) #subset data to show species significant at 0.05

head(spp.scrs)

sig.spp.scrs
# To show biologically extrinsic variables another datasheet needs to be created

env.scores <- as.data.frame(scores(env.fit, display = "vectors")) #extracts relevant scores from envifit
env.scores <- cbind(env.scores, env.variables = rownames(env.scores)) #and then gives them their names

env.scores <- cbind(env.scores, pval = env.fit$vectors$pvals) # add pvalues to dataframe
sig.env.scrs <- subset(env.scores, pval<=0.05) #subset data to show variables significant at 0.05

head(env.scores)
sig.env.scrs
env.scores
#Now we have the relevant information for plotting the ordination in ggplot! Lets get plotting!
#######################################
# add the hull
surrounding.rural <- sample.scrs[sample.scrs$surrounding == "rural", ][chull(sample.scrs[sample.scrs$surrounding == 
                                                                                           "rural", c("NMDS1", "NMDS2")]), ]  # hull values for rural areas
surrounding.urban <- sample.scrs[sample.scrs$surrounding == "urban", ][chull(sample.scrs[sample.scrs$surrounding == 
                                                                                           "urban", c("NMDS1", "NMDS2")]), ]  # hull values for urban areas


hull.data <- rbind(surrounding.rural, surrounding.urban)  #combine sampling site groups in the same dataset
hull.data

#########################################
# Basic ordination plot (restoration_date and season)
nmds.plot.diatoms <- ggplot(sample.scrs, aes(x=NMDS1, y=NMDS2), size = 7)+ #sets up the plot
  geom_point(aes(NMDS1, NMDS2, colour = factor(surrounding), shape = factor(surrounding)), size = 6)+ #adds site points to plot, shape determined season, colour determined by restorartion_date
  #geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS2,fill=surrounding,group=surrounding),alpha=0.30) +
  stat_ellipse(aes(fill=surrounding), alpha=.2,type='t',size =1, geom="polygon")+ ## add ellipses
  #geom_text(data=sample.scrs,aes(x=NMDS1,y=NMDS2,label=sample), vjust=0) +  # add the site labels
  coord_fixed()+
  theme_classic()+ 
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+
  labs(colour = "surrounding", shape = "surrounding")+ # add legend labels for Management and Landuse
  theme(legend.position = "right", legend.text = element_text(size = 14), legend.title = element_text(size = 14), axis.text = element_text(size = 16)) + # add legend at right of plot
  annotate(geom = "label", x = -1.5, y = -2, size = 5,
           label = paste("Stress: ", round(spp.nmds$stress, digits = 3)))# add stress value

nmds.plot.diatoms + labs(title = "") # displays plot

# Significant species
nmds1 <- nmds.plot.diatoms +
  geom_segment(data = sig.spp.scrs, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", lwd=0.3) + #add vector arrows of significant species
  ggrepel::geom_text_repel(data = sig.spp.scrs, aes(x=NMDS1, y=NMDS2, label = Species), cex = 3, direction = "both", segment.size = 0.25)+ #add labels for species, use ggrepel::geom_text_repel so that labels do not overlap
  labs(title = "")

print(nmds1)

# Significant abiotic variables
nmds.plot.diatoms +
  geom_segment(data = sig.env.scrs, aes(x = 0, xend=NMDS1, y = 0, yend=NMDS2), arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", lwd=0.3) + #add vector arrows of significant env variables
  ggrepel::geom_text_repel(data = sig.env.scrs, aes(x=NMDS1, y=NMDS2, label = env.variables), cex = 4, direction = "both", segment.size = 0.25)+ #add labels for env variables
  labs(title="")

# add both significant taxa and abiotic factors
kinzig2021_nmds_DM_paper <- nmds1 +
  geom_segment(data = sig.env.scrs, aes(x = 0, xend=NMDS1, y = 0, yend=NMDS2), arrow = arrow(length = unit(0.25, "cm")), colour = "red", lwd=0.3) + #add vector arrows of significant env variables
  ggrepel::geom_text_repel(data = sig.env.scrs, aes(x=NMDS1, y=NMDS2, label = env.variables), cex = 4, direction = "both", colour = "red", segment.size = 0.25)+ #add labels for env variables
  labs(title="")
kinzig2021_nmds_DM_paper

# save the plot
ggsave("kinzig2021_nmds_DM3%_paper2.tiff", units="in", width=7, height=7, dpi = 300, compression = 'lzw')

########
# nMDS --------------------------------------------------------------------
# 18SV9 amplicon sequencing data

set.seed(123)
spp.log_dna <- decostand(dna4, method = "log")
# presence/absence
spp.pa_dna <- decostand(dna4, method = "pa")

#spp.log1 <- decostand(spe5, method = "log")
spp.log.dis_dna <- vegdist(spp.log_dna, method = "bray")
# presence/absence
spp.pa.dis_dna <- vegdist(spp.pa_dna, method = "bray")

spp.log_dna
spp.pa_dna
spp.log.dis_dna
spp.pa.dis_dna
# check
decorana (spp.log_dna) # The length of the first DCA axis is 2.7 S.D. (i.e. < 3. S.D.), 
decorana(spp.pa_dna)
# and these data are thus suitable for linear ordination methods.


# betadisper --------------------------------------------------------------

# Before doing the PERMANOVA, first we check to see if the dispersion is the same
# Homogeneity of groups
# betadisper studies the differences in group homogeneities
# analogous to Levene's test of the equality of variances
# can only use one factor as an independent variable

# surrounding
(mod.surrounding <- with(data, betadisper(spp.log.dis_dna, data$surrounding)))
plot(mod.surrounding, sub = NULL) 
boxplot(mod.surrounding)
anova(mod.surrounding)   # this says that within-group variances can be considered homogeneous - in fact a good news!
# next step could be to compare the location of centroids and see if there is
# significant difference among groups in that respect!
permutest(mod.surrounding)

# for the Bray-Curtis distances of the relative abundances:
# for relative abundances, e.g. Bray-Curtis is often used:
library(MASS)
bray.dist <- vegdist(dna4)

bray.pco <- isoMDS(bray.dist)
plot(bray.pco$points[,1:2], pch = 16, col = data$surrounding, main = "PCoA of Bray distances")
legend("bottomright", legend = levels(data$surrounding), pch = 16, col = 1:5)

###########################################################################

# PERMANOVA ---------------------------------------------------------------

# Permutational multivariate analysis of variance using distance matrices
# (Bray-Curtis similarities by default). ANOSIM uses only ranks of Bray-Curtis,
# so the former preserves more information.

(perm.1 <- adonis2(spp.log.dis_dna~(surrounding*pH*water_temp.*oxygen_full*dissolved_oxygen*oxygen_saturation*conductivity*average_velocity*average_depth*discharge*water_level_at_pole),
                   method = perm, data = data))

(perm.1 <- adonis2(spp.pa.dis_dna~(surrounding+pH+water_temp.+oxygen_full+dissolved_oxygen+oxygen_saturation+conductivity+average_velocity+average_depth+discharge+water_level_at_pole),
                   method = perm, data = data))

## great - this tells us that salinity does have a significant effect
# but not velocity
# so the simplest model we can take would be:
(perm.2 <- adonis2(spp.log.dis_dna ~ surrounding,
                   method = perm, data = data))
(perm.2 <- adonis2(spp.pa.dis_dna ~ surrounding,
                   method = perm, data = data))


(perm.2 <- adonis(spp.log.dis_dna ~ surrounding,
                   method = perm, data = data))
(perm.2 <- adonis2(spp.pa.dis_dna ~ surrounding,
                   method = perm, data = data))
##################################################
# NMDS

spp.nmds_dna <- metaMDS(spp.log_dna, k = 2,trymax = 1000, permutation = 9999,
                        distance = "bray", wascores = TRUE)

spp.nmds_dna_pa <- metaMDS(spp.pa_dna, k = 2,trymax = 1000, permutation = 9999,
                           distance = "bray", wascores = TRUE)

spp.nmds_dna
stressplot(spp.nmds_dna)

# fitting environmental variables
env.fit_dna <- envfit(spp.nmds_dna, data, permutations = 999, na.rm = TRUE) # this fits environmental vectors
spp.fit_dna <- envfit(spp.nmds_dna, dna4, permutations = 999) # this fits species vectors

env.fit_dna

# To plot the output from the mds using ggplot a new datasheet needs to be created which contains the x,y points for each site. 
#You can do this by calling the scores of you mds.

sample.scrs_dna <- as.data.frame(scores(spp.nmds_dna, display = "sites")) #save NMDS results into dataframe
sample.scrs_dna <- cbind(sample.scrs_dna, surrounding = data$surrounding) #add grouping variable "Management" to dataframe
#sample.scrs <- cbind(sample.scrs, year.of.restor. = metadata1$year.of.r) #add grouping variable of cluster grouping to dataframe
#site.scrs <- cbind(site.scrs, Site = rownames(site.scrs)) #add site names as variable if you want to display on plot
#sample.scrs <- cbind(sample.scrs, location = env$site)

head(sample.scrs_dna)
#view(sample.scrs)

# A new dataset containing species data also needs to be made to look at species vectors.
#This is not necessary if you don't want to show the species on the final graph.
spp.scrs_dna <- as.data.frame(scores(spp.fit_dna, display = "vectors")) #save species intrinsic values into dataframe
spp.scrs_dna <- cbind(spp.scrs_dna, Species = rownames(spp.scrs_dna)) #add species names to dataframe
spp.scrs_dna <- cbind(spp.scrs_dna, pval = spp.fit_dna$vectors$pvals) #add pvalues to dataframe so you can select species which are significant
#spp.scrs<- cbind(spp.scrs, abrev = abbreviate(spp.scrs$Species, minlength = 6)) #abbreviate species names
sig.spp.scrs_dna <- subset(spp.scrs_dna, pval<=0.05) #subset data to show species significant at 0.05

head(spp.scrs_dna)
#view(spp.scrs)
sig.spp.scrs_dna
# To show biologically extrinsic variables another datasheet needs to be created

env.scores_dna <- as.data.frame(scores(env.fit_dna, display = "vectors")) #extracts relevant scores from envifit
env.scores_dna <- cbind(env.scores_dna, env.variables = rownames(env.scores_dna)) #and then gives them their names

env.scores_dna <- cbind(env.scores_dna, pval = env.fit_dna$vectors$pvals) # add pvalues to dataframe
sig.env.scrs_dna <- subset(env.scores_dna, pval<=0.05) #subset data to show variables significant at 0.05

head(env.scores_dna)
sig.env.scrs_dna
#Now we have the relevant information for plotting the ordination in ggplot! Lets get plotting!
#######################################
# add the hull
surrounding.rural <- sample.scrs[sample.scrs$surrounding == "rural", ][chull(sample.scrs[sample.scrs$surrounding == 
                                                                                           "rural", c("NMDS1", "NMDS2")]), ]  # hull values for rural areas
surrounding.urban <- sample.scrs[sample.scrs$surrounding == "urban", ][chull(sample.scrs[sample.scrs$surrounding == 
                                                                                           "urban", c("NMDS1", "NMDS2")]), ]  # hull values for urban areas


hull.data <- rbind(surrounding.rural, surrounding.urban)  #combine sampling site groups in the same dataset
hull.data

#########################################
# Basic ordination plot (restoration_date and season)
nmds.plot.dna <- ggplot(sample.scrs_dna, aes(x=NMDS1, y=NMDS2), size = 7)+ #sets up the plot
  geom_point(aes(NMDS1, NMDS2, colour = factor(surrounding), shape = factor(surrounding)), size = 6)+ #adds site points to plot, shape determined season, colour determined by restorartion_date
  #geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS2,fill=surrounding,group=surrounding),alpha=0.30) +
  stat_ellipse(aes(fill=surrounding), alpha=.2,type='t',size =1, geom="polygon")+ ## add ellipses
  #geom_text(data=sample.scrs,aes(x=NMDS1,y=NMDS2,label=sample), vjust=0) +  # add the site labels
  coord_fixed()+
  theme_classic()+ 
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+
  labs(colour = "surrounding", shape = "surrounding")+ # add legend labels for Management and Landuse
  theme(legend.position = "right", legend.text = element_text(size = 14), legend.title = element_text(size = 14), axis.text = element_text(size = 16)) + # add legend at right of plot
  annotate(geom = "label", x = -0.5, y = -2, size = 5,
           label = paste("Stress: ", round(spp.nmds_dna$stress, digits = 3)))# add stress value

nmds.plot.dna + labs(title = "") # displays plot

# Significant species
#tiff("EDM2023_nmdsyear", units="in", width=10, height=10, res=600)
nmds2 <- nmds.plot.dna +
  geom_segment(data = sig.spp.scrs_dna, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", lwd=0.3) + #add vector arrows of significant species
  ggrepel::geom_text_repel(data = sig.spp.scrs_dna, aes(x=NMDS1, y=NMDS2, label = Species), cex = 3, direction = "both", segment.size = 0.25)+ #add labels for species, use ggrepel::geom_text_repel so that labels do not overlap
  labs(title = "")

print(nmds2)

# Significant abiotic variables
nmds.plot.dna +
  geom_segment(data = sig.env.scrs_dna, aes(x = 0, xend=NMDS1, y = 0, yend=NMDS2), arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", lwd=0.3) + #add vector arrows of significant env variables
  ggrepel::geom_text_repel(data = sig.env.scrs_dna, aes(x=NMDS1, y=NMDS2, label = env.variables), cex = 4, direction = "both", segment.size = 0.25)+ #add labels for env variables
  labs(title="")


kinzig2021_nmds_DNA2_paper <- nmds2 +
  geom_segment(data = sig.env.scrs_dna, aes(x = 0, xend=NMDS1, y = 0, yend=NMDS2), arrow = arrow(length = unit(0.25, "cm")), colour = "red", lwd=0.3) + #add vector arrows of significant env variables
  ggrepel::geom_text_repel(data = sig.env.scrs_dna, aes(x=NMDS1, y=NMDS2, label = env.variables), cex = 4, direction = "both", colour = "red", segment.size = 0.25)+ #add labels for env variables
  labs(title="")

kinzig2021_nmds_DNA2_paper

# save high resolution ordination plot
ggsave("kinzig2021_nmds_DNA_data_paper2.tiff", units="in", width=7, height=7, dpi = 300, compression = 'lzw')
###################################################
# calculate alpha diversity
# Keep only taxa that occurred at least once

diatoms <- spe
# make site_name column as rownames
diatoms1 <- diatoms %>% column_to_rownames(var="site_name") # change row names
good.taxa <- colSums(diatoms1) > 0

diatoms2 <- diatoms1[,good.taxa]

# Compute alpha diversity indices of the diatom communities
# *******************************************************

# Get help on the diversity() function
#library(vegan)
?diversity

N0 <- rowSums(diatoms2 > 0)         # Species richness
N0a <- specnumber(diatoms2)     # Species richness (alternate)
H <- diversity(diatoms2)       # Shannon entropy (base e)
Hb2 <- diversity(diatoms2, base = 2) # Shannon entropy (base 2)
N1 <- exp(H)                   # Shannon diversity (number of abundant species) (base e)
N1b2 <- 2^Hb2                 # Shannon diversity (base 2)
N2 <- diversity(diatoms2, "inv")    # Simpson diversity (number of dominant species)
J <- H/log(N0)                 # Pielou evenness
E10 <- N1/N0                   # Shannon evenness (Hill's ratio)
E20 <- N2/N0                   # Simpson evenness (Hill's ratio)
(div <- data.frame(N0, N0a, H, Hb2, N1, N1b2, N2, E10, E20, J))
# make div dataset rownames as first column
div1 <- tibble::rownames_to_column(div, "site_name")
# assigning new column names to some columns in div dataframe
colnames(div1)[2] <- "Richness"
colnames(div1)[4] <- "Shannon"
colnames(div1)[8] <- "Simpson"
colnames(div1)[11] <- "Pielou"

write.csv(div1, file = "div_kinzig2021.csv")

# Combine this div dataset with environmental variables dataset
div2 <- left_join(div1, data)

range(div2$Richness)
range(div2$Shannon)
# select only rural sites
div_rural <- div2 %>%
  filter(surrounding == "rural")
range(div_rural$Richness)
range(div_rural$Shannon)
mean(div_rural$Richness)
mean(div_rural$Shannon)
# Select only urban sites 
div_urban <- div2 %>%
  filter(surrounding == "urban")
range(div_urban$Richness)
range(div_urban$Shannon)
mean(div_urban$Richness)
mean(div_urban$Shannon)

# remove column 3 in div2
div2 <- div2[,-3]


# run multiple comparisons for analysis of variance (t-test or ANOVA)
# Multiple comparisons
# install.packages("ggstatsplot")
library(ggstatsplot)
# Comparison between species

# edit from here
x <- "surrounding"
cols <- 2:3 # the 4 continuous dependent variables
type <- "parametric" # given the large number of observations, we use the parametric version
paired <- FALSE # FALSE for independent samples, TRUE for paired samples

# edit until here

# edit at your own risk
plotlist <-
  purrr::pmap(
    .l = list(
      data = list(as_tibble(div2)),
      x = x,
      y = as.list(colnames(div2)[cols]),
      plot.type = "box", # for boxplot
      type = type, # parametric or nonparametric
      pairwise.comparisons = FALSE, # to run post-hoc tests if more than 2 groups
      pairwise.display = "significant", # show only significant differences
      bf.message = FALSE, # remove message about Bayes Factor
      centrality.plotting = TRUE # remove central measure
    ),
    .f = ifelse(paired, # automatically use ggwithinstats if paired samples, ggbetweenstats otherwise
                ggstatsplot::ggwithinstats,
                ggstatsplot::ggbetweenstats
    )
  )

# print all plots together with statistical results
for (i in 1:length(plotlist)) {
  print(plotlist[[i]] +
          labs(caption = NULL)) # remove caption
}

#install.packages("patchwork")
library(patchwork)
tiff("alpha_div2_data_micro.tiff", units="in", width=15, height=10, res=300, compression = 'lzw')
wrap_plots(plotlist) +
  plot_layout(ncol = 2) +
  plot_annotation(
    tag_levels = "a",
    tag_suffix = ")"
  )
dev.off()

# computing welch t-test
# Richness
library(rstatix)
stat.test1 <- div2 %>%
  t_test(Richness ~ surrounding) %>%
  add_significance()
stat.test1

# Effect size: Cohen's d formula for Welch t-test
div2 %>% cohens_d(Richness ~ surrounding, var.equal = FALSE)

# Shannon
stat.test1 <- div2 %>%
  t_test(Shannon ~ surrounding) %>%
  add_significance()
stat.test1

# Effect size: Cohen's d formula for Welch t-test
div2 %>% cohens_d(Shannon ~ surrounding, var.equal = FALSE)

###########################################################
# table of diatoms >= 5%

# select all columns except "method"
#tab <- diat %>% dplyr::select(-method)
# make the dataset tidy
tab<-diat
tab1 <- tab %>% 
  pivot_longer(!c(site_name, method),  names_to = "taxa", values_to = "proportion")

# Keep only taxa with proportion >= 5, or 5% of the total relative abundance per sample.
tab2 <- subset(tab1, proportion >= 5)
# Now we can spread the dataset of the most dominant taxa with at least 5% of proportion per sample

tab3 <- tab2 %>%
  pivot_wider(names_from = taxa, values_from = proportion, values_fn = sum, values_fill = 0)

tab4 <- left_join(tab3, metadata2)
write.csv(tab4, file = "kinzig2021_table_paper.csv")

#############################################################
##### alpha div on diatom 18S data
######################################
# spread the full diatom 18S dataset
dna_div <- dna1 %>% 
  pivot_wider(names_from = taxa, values_from = proportion, values_fn = sum, values_fill = 0)
# make site_name column as rownames
dna_div1 <- dna_div %>% column_to_rownames(var="site_name") # change row names
# Keep only taxa that occurred at least once

good.taxa <- colSums(dna_div1) > 0

dna_div2 <- dna_div1[,good.taxa]


N0 <- rowSums(dna_div2 > 0)         # Species richness
N0a <- specnumber(dna_div2)     # Species richness (alternate)
H <- diversity(dna_div2)       # Shannon entropy (base e)
Hb2 <- diversity(dna_div2, base = 2) # Shannon entropy (base 2)
N1 <- exp(H)                   # Shannon diversity (number of abundant species) (base e)
N1b2 <- 2^Hb2                 # Shannon diversity (base 2)
N2 <- diversity(dna_div2, "inv")    # Simpson diversity (number of dominant species)
J <- H/log(N0)                 # Pielou evenness
E10 <- N1/N0                   # Shannon evenness (Hill's ratio)
E20 <- N2/N0                   # Simpson evenness (Hill's ratio)
(div_dna <- data.frame(N0, N0a, H, Hb2, N1, N1b2, N2, E10, E20, J))
# make div dataset rownames as first column
div_dna1 <- tibble::rownames_to_column(div_dna, "site_name")

# assigning new column names to some columns in div dataframe
colnames(div_dna1)[2] <- "Richness"
colnames(div_dna1)[4] <- "Shannon"
colnames(div_dna1)[8] <- "Simpson"
colnames(div_dna1)[11] <- "Pielou"
# Combine this div dataset with envir
div_dna2 <- left_join(div_dna1, data)
range(div_dna2$Richness)
range(div_dna2$Shannon)
# remove column 3 in div2
div_dna2 <- div_dna2[,-3]
# select only rural sites
div_rural <- div2 %>%
  filter(surrounding == "rural")
range(div_rural$Richness)
range(div_rural$Shannon)
mean(div_rural$Richness)
mean(div_rural$Shannon)
# Select only urban sites 
div_urban <- div2 %>%
  filter(surrounding == "urban")
range(div_urban$Richness)
range(div_urban$Shannon)
mean(div_urban$Richness)
mean(div_urban$Shannon)

# run multiple comparisons for analysis of variance
# Multiple comparisons
# install.packages("ggstatsplot")
#library(ggstatsplot)
# Comparison between species

# edit from here
x <- "surrounding"
cols <- 2:3 # the 4 continuous dependent variables
type <- "parametric" # given the large number of observations, we use the parametric version
paired <- FALSE # FALSE for independent samples, TRUE for paired samples

# edit until here

# edit at your own risk
plotlist1 <-
  purrr::pmap(
    .l = list(
      data = list(as_tibble(div_dna2)),
      x = x,
      y = as.list(colnames(div_dna2)[cols]),
      plot.type = "box", # for boxplot
      type = type, # parametric or nonparametric
      pairwise.comparisons = FALSE, # to run post-hoc tests if more than 2 groups
      pairwise.display = "significant", # show only significant differences
      bf.message = FALSE, # remove message about Bayes Factor
      centrality.plotting = TRUE # remove central measure
    ),
    .f = ifelse(paired, # automatically use ggwithinstats if paired samples, ggbetweenstats otherwise
                ggstatsplot::ggwithinstats,
                ggstatsplot::ggbetweenstats
    )
  )

# print all plots together with statistical results
for (i in 1:length(plotlist1)) {
  print(plotlist1[[i]] +
          labs(caption = NULL)) # remove caption
}

#library(patchwork)

wrap_plots(plotlist) +
  plot_layout(ncol = 2) +
  plot_annotation(
    tag_levels = "a",
    tag_suffix = ")"
  )

library(ggpubr)
alpha_div <- ggarrange(plotlist[[1]],plotlist[[2]], plotlist1[[1]],plotlist1[[2]], 
                       labels = c("a)", "b)", "c)", "d)"),
                       ncol = 2, nrow = 2)

alpha_div
ggsave("alpha_div.tiff", units="in", width=10, height=7,  dpi=300, compression = 'lzw')
?ggarrange
################################################
# indices using presence/absence data
ind <- read.csv("kinzig2021_indices_rerun.csv", sep=",", 
                header=T)
# combine indices dataset and env. data
ind <- left_join(ind, data)
# filter microscopy
ind_mic <- ind %>%
  filter(method == "MIC")
# filter dna
ind_dna <- ind %>%
  filter(method == "DNA")

## Method == microscopy
# select only rural sites
ind_mic_rural <- ind_mic %>%
  filter(surrounding == "rural")
range(ind_mic_rural$IPS)
mean(ind_mic_rural$IPS)
# Select only urban sites 
ind_mic_urban <- ind_mic %>%
  filter(surrounding == "urban")
range(ind_mic_urban$IPS)
mean(ind_mic_urban$IPS)

## Method == 18SV9 DNA
# select only rural sites
ind_dna_rural <- ind_dna %>%
  filter(surrounding == "rural")
range(ind_dna_rural$IPS)
mean(ind_dna_rural$IPS)
# Select only urban sites 
ind_dna_urban <- ind_dna %>%
  filter(surrounding == "urban")
range(ind_dna_urban$IPS)
mean(ind_dna_urban$IPS)

##############################
# MIC
x <- "surrounding"
cols <- 8 # the 4 continuous dependent variables
type <- "parametric" # given the large number of observations, we use the parametric version
paired <- FALSE # FALSE for independent samples, TRUE for paired samples

# edit until here

# edit at your own risk
plotlist_a <-
  purrr::pmap(
    .l = list(
      data = list(as_tibble(ind_mic)),
      x = x,
      y = as.list(colnames(ind_mic)[cols]),
      plot.type = "box", # for boxplot
      type = type, # parametric or nonparametric
      pairwise.comparisons = FALSE, # to run post-hoc tests if more than 2 groups
      pairwise.display = "significant", # show only significant differences
      bf.message = FALSE, # remove message about Bayes Factor
      centrality.plotting = TRUE # remove central measure
    ),
    .f = ifelse(paired, # automatically use ggwithinstats if paired samples, ggbetweenstats otherwise
                ggstatsplot::ggwithinstats,
                ggstatsplot::ggbetweenstats
    )
  )

# print all plots together with statistical results
for (i in 1:length(plotlist_a)) {
  print(plotlist_a[[i]] +
          labs(caption = NULL)) # remove caption
}
# 18SV9
# edit from here
x <- "surrounding"
cols <- 8 # the 4 continuous dependent variables
type <- "parametric" # given the large number of observations, we use the parametric version
paired <- FALSE # FALSE for independent samples, TRUE for paired samples

# edit until here

# edit at your own risk
plotlist_b <-
  purrr::pmap(
    .l = list(
      data = list(as_tibble(ind_dna)),
      x = x,
      y = as.list(colnames(ind_dna)[cols]),
      plot.type = "box", # for boxplot
      type = type, # parametric or nonparametric
      pairwise.comparisons = FALSE, # to run post-hoc tests if more than 2 groups
      pairwise.display = "significant", # show only significant differences
      bf.message = FALSE, # remove message about Bayes Factor
      centrality.plotting = TRUE # remove central measure
    ),
    .f = ifelse(paired, # automatically use ggwithinstats if paired samples, ggbetweenstats otherwise
                ggstatsplot::ggwithinstats,
                ggstatsplot::ggbetweenstats
    )
  )

# print all plots together with statistical results
for (i in 1:length(plotlist_b)) {
  print(plotlist_b[[i]] +
          labs(caption = NULL)) # remove caption
}

#install.packages("patchwork")
#library(patchwork)
library(ggpubr)
IPS_plot <- ggarrange(plotlist_a[[1]], plotlist_b[[1]], 
                      labels = c("a)", "b)"),
                      ncol = 2, nrow = 1)

IPS_plot
ggsave("IPS_plot.tiff", units="in", width=10, height=7, dpi=300, compression = 'lzw')

#wrap_plots(plotlist, plotlist1) +
# plot_layout(ncol = 2) +
#plot_annotation(
# tag_levels = "a",
#tag_suffix = ")"
#)
#dev.off()

