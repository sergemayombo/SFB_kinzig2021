# Analysis of raw amplicon sequencing data
# Kinzig 2021 sampling campaign
# Natrix pipeline output data
# Multivariate statistical analysis
# 30 August 2023
# Serge Mayombo

# load data in R

#library(readr)
library(vegan)
library(tidyverse)

# Sequencing data
# @Serge: I simplified the below part a bit
# I hope I got the important steps!
protists <- read.table("kinzig2021_18SV9_filtered.csv", sep=",", 
                  header=T, row.names = 1)

view(mol[grep("Metazoa", mol$Phylum),])


########
# remove metazoa
protists <- mol[mol$Phylum != "Metazoa",]

# diatom subset
dna_diat <- subset(protists, Class == "Bacillariophyta")
# write.csv(protists, file = "diatoms_kinzig2021.csv")
# after formating dna_diat_kinzig2021.csv on Excel
#diat_list <- read.csv("dna_diatom_kinzig2021_list.csv")
###### persence/absence
spp_dna.pa <- decostand(dna_diat[,1:19], method = "pa")
#write.csv(spp.pa, file = "kinzig2021_dna_pa")
# summary pa
# summary
sum_dna.pa <- colSums(spp_dna.pa)
view(sum_dna.pa)
sum(sum_dna.pa)
range(sum_dna.pa)
mean(sum_dna.pa)


diat_count <- colSums(dna_diat[,1:19])
view(diat_count)
sum(diat_count)

# readcounts full protists dataset

readcounts <- protists[,1:19]
dim(readcounts)

# transpose dataset

dna_spe_trans <- t(readcounts)
total_reads_count <- rowSums(dna_spe_trans)
view(total_reads_count)
sum(total_reads_count)
range(total_reads_count)
mean(total_reads_count)

## diatom reads
diat_count <- t(diat_count)
view(diat_count)
sum(diat_count)
range(diat_count)
mean(diat_count)
# 18SV9
# Keep only diatom OTUs with readcount > 0


reads <- function(x){
  if(is.numeric(x)){
    sum(x) > 0
  } else {
    TRUE
  }
}

diat_good_18S <- diat_count[, sapply(diat_count,  reads)]

diat_good <- as.data.frame(diat_good_18S)

# arranges the rows in reverse order in the dna_spe_trans dataframe


#dna_spe_tidy <- readcounts %>% gather(key = sample, value = abundance, -species)

sapply(dna_spe_trans, class)

#taxa
spec <- protits[,20:28]
spec_trans <- t(spec)
sapply(spec_trans, class)

# environmental data

data <- as_tibble(read_csv("kinzig2021_fieldData.csv"))
str(data)
data <-  data[-18,]
# change env variable classes

data$surrounding <- as.factor(data$surrounding)

# change row names
#data <- data %>% remove_rownames %>% column_to_rownames(var="site_name") # change row names

# Test for the effect of salinity and velocity at species level

#######################################################################
#install.packages("vegan")

#library(vegan)

# Logarithmic transformation as suggested by Anderson et al. (2006): log_b (x) + 1 for x > 0, where b is the base of 
# the logarithm; zeros are left as zeros. 
# Higher bases give less weight to quantities and more to presences.
spp.log_alg <- decostand(dna_spe_trans, method = "log")
spp.log.dis_alg <- vegdist(spp.log_alg, method = "bray")


# betadisper --------------------------------------------------------------

# Before doing the PERMANOVA, first we check to see if the dispersion is the same
# Homogeneity of groups
# betadisper studies the differences in group homogeneities
# analogous to Levene's test of the equality of variances
# can only use one factor as an independent variable

# surrounding
(mod.surrounding <- with(data, betadisper(spp.log.dis_alg, data$surrounding)))
plot(mod.surrounding, sub = NULL) 
boxplot(mod.surrounding)
anova(mod.surrounding)   # this says that within-group variances can be considered homogeneous - in fact a good news!
# next step could be to compare the location of centroids and see if there is
# significant difference among groups in that respect!
permutest(mod.surrounding)

# PERMANOVA ---------------------------------------------------------------

# Permutational multivariate analysis of variance using distance matrices
# (Bray-Curtis similarities by default). ANOSIM uses only ranks of Bray-Curtis,
# so the former preserves more information. 

(perm.1 <- adonis2(spp.log.dis_alg~(surrounding+pH+water_temp.+oxygen_full+dissolved_oxygen+oxygen_saturation+conductivity+average_velocity+average_depth+discharge+water_level_at_pole),
                   method = perm, data = data))

# so the simplest model we can take would be:
(perm.2 <- adonis2(spp.log.dis ~ surrounding,
                   method = perm, data = data))

# the distribution if it is a normal one
plot (density (spp.log.dis_alg)) 
# surrounding has no effect on protist community structure
# to visualise that we run the NMDS
# nMDS --------------------------------------------------------------------

spp.nmds_alg <- metaMDS(spp.log_alg, k = 2,trymax = 1000,
                        distance = "bray", wascores = TRUE, permutations = 999,
                        plot=T)
# I increased the number of permutations because it did not converge
spp.nmds_alg
stressplot(spp.nmds_alg) # Large scatter around the line indicates that original dissimilarities 
# are not well preserved in the reduced number of dimensions. But it looks pretty good in our case.
# nmds plot

#############################

env.fit_alg <- envfit(spp.nmds_alg, data, permutations = 999, na.rm = TRUE) # this fits environmental vectors
spp.fit_alg <- envfit(spp.nmds_alg, dna_spe_trans, permutations = 999) # this fits species vectors

env.fit_alg
spp.fit_alg
# To plot the output from the mds using ggplot a new datasheet needs to be created which contains the x,y points for each site. 
#You can do this by calling the scores of you mds.
sample.scrs_alg <- as.data.frame(scores(spp.nmds_alg, display = "sites")) #save NMDS results into dataframe
sample.scrs_alg <- cbind(sample.scrs_alg, surrounding = data$surrounding) #add grouping variable "Management" to dataframe

head(sample.scrs_alg)

#Now we have the relevant information for plotting the ordination in ggplot! Lets get plotting!
#######################################
# add the hull
surrounding.rural <- sample.scrs_alg[sample.scrs_alg$surrounding == "rural", ][chull(sample.scrs_alg[sample.scrs_alg$surrounding == 
                                                                                                       "rural", c("NMDS1", "NMDS2")]), ]  # hull values for rural areas
surrounding.urban <- sample.scrs_alg[sample.scrs_alg$surrounding == "urban", ][chull(sample.scrs_alg[sample.scrs_alg$surrounding == 
                                                                                                       "urban", c("NMDS1", "NMDS2")]), ]  # hull values for grp B

hull.data <- rbind(surrounding.rural, surrounding.urban)  #combine grp.a and grp.b
hull.data

# A new dataset containing species data also needs to be made to look at species vectors.
#This is not necessary if you don't want to show the species on the final graph.
# can't fit all OTUs here, the plot will be too messy
spp.scrs_alg <- as.data.frame(scores(spp.fit_alg, display = "vectors")) #save species intrinsic values into dataframe
spp.scrs_alg <- cbind(spp.scrs_alg, Species = rownames(spp.scrs_alg)) #add species names to dataframe
spp.scrs_alg <- cbind(spp.scrs_alg, pval = spp.fit_alg$vectors$pvals) #add pvalues to dataframe so you can select species which are significant
#spp.scrs<- cbind(spp.scrs, abrev = abbreviate(spp.scrs$Species, minlength = 6)) #abbreviate species names
sig.spp.scrs_alg <- subset(spp.scrs_alg, pval<=0.05) #subset data to show species significant at 0.05

head(spp.scrs_alg)
sig.spp.scrs_alg

# To show biologically extrinsic variables another datasheet needs to be created

env.scores_alg <- as.data.frame(scores(env.fit_alg, display = "vectors")) #extracts relevant scores from envifit
env.scores_alg <- cbind(env.scores_alg, env.variables = rownames(env.scores_alg)) #and then gives them their names

env.scores_alg <- cbind(env.scores_alg, pval = env.fit_alg$vectors$pvals) # add pvalues to dataframe
sig.env.scrs_alg <- subset(env.scores_alg, pval<=0.05) #subset data to show variables significant at 0.05

head(env.scores_alg)
sig.env.scrs_alg
#Now we have the relevant information for plotting the ordination in ggplot! Lets get plotting!
#dev.off()
# Basic ordination plot (restoration_date and season)
nmds.plot.dna_alg <- ggplot(sample.scrs_alg, aes(x=NMDS1, y=NMDS2), size = 7)+ #sets up the plot
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
           label = paste("Stress: ", round(spp.nmds_alg$stress, digits = 3)))# add stress value

nmds.plot.dna_alg + labs(title = "") # displays plot

# Significant species
nmds2_alg <- nmds.plot.dna_alg +
  geom_segment(data = sig.spp.scrs_alg, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", lwd=0.3) + #add vector arrows of significant species
  ggrepel::geom_text_repel(data = sig.spp.scrs_alg, aes(x=NMDS1, y=NMDS2, label = Species), cex = 3, direction = "both", segment.size = 0.25)+ #add labels for species, use ggrepel::geom_text_repel so that labels do not overlap
  labs(title = "")

print(nmds2_alg)

# Significant abiotic variables
nmds.plot.dna_alg +
  geom_segment(data = sig.env.scrs_alg, aes(x = 0, xend=NMDS1, y = 0, yend=NMDS2), arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", lwd=0.3) + #add vector arrows of significant env variables
  ggrepel::geom_text_repel(data = sig.env.scrs_alg, aes(x=NMDS1, y=NMDS2, label = env.variables), cex = 4, direction = "both", segment.size = 0.25)+ #add labels for env variables
  labs(title="")


kinzig2021_nmds_DNA_alg <- nmds.plot.dna_alg +
  geom_segment(data = sig.env.scrs_alg, aes(x = 0, xend=NMDS1, y = 0, yend=NMDS2), arrow = arrow(length = unit(0.25, "cm")), colour = "red", lwd=0.3) + #add vector arrows of significant env variables
  ggrepel::geom_text_repel(data = sig.env.scrs_alg, aes(x=NMDS1, y=NMDS2, label = env.variables), cex = 4, direction = "both", colour = "red", segment.size = 0.25)+ #add labels for env variables
  labs(title="")

kinzig2021_nmds_DNA_alg

ggsave("kinzig2021_nmds_DNA_protists_paper.tiff", units="in", width=7, height=7, dpi = 300, compression = 'lzw')
#dev.new()


# So we found no evidence that surrounding land use affect the abundance of protists
# as shown in PERMANOVA result


#########################################################################
############
# plotting taxa proportions
############
# Remove samples with 0 diatom abundances
# Remove samples with 0 diatom abundances
# make rownames as first column

dna <- protists %>% rownames_to_column(var="seq_IDs") # change row names


set.seed(123)
read_ab <- protists[as.logical(rowSums(protists[,1:19]!= 0)), ]
colSums(read_ab[,1:19])
x<-read_ab[,1:19]
x <- t(x)
x<-x/rowSums(x)
x<-x[,order(colSums(x),decreasing=TRUE)]
colSums(x)
#view(x)
#colSums(x1)
head(x)

#Extract list of top N Taxa
N<-40
taxa_list<-protists[,24]
#taxa_list1<-colnames(x1)[1:N]
N<-length(taxa_list)
#N1<-length(taxa_list1)

#Create a custom color scale
library(RColorBrewer)
# Define the number of colors you want
nb.cols <- 40
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)

# other colours
library(RColorBrewer)
#install.packages("randomcoloR")
library(randomcoloR)
n <- 40
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

df<-NULL
for (i in 1:dim(x)[2]){
  tmp<-data.frame(row.names=NULL,Sample=rownames(data), surrounding=data$surrounding,  
                  Class=rep(protists$Class[i],dim(x)[1]),Value=x[,i])
  if(i==1){df<-tmp} else {df<-rbind(df,tmp)}
}

# We select only OTUs with read abundances 3% of the total counts per sample.

df_prot <- subset(df, Value >= 0.03)

# Based on site location along the river system: rural vs. urban

txp_dna_algae2_kinzig2021 <- ggplot(df_prot, aes(x=Sample,y=Value, fill=Class))+geom_bar(position="fill", stat="identity")+
  facet_grid(. ~ surrounding, drop=TRUE,scale="free",space="free_x")+
  scale_fill_manual(values=palette[1:(N+1)]) + theme_bw()+ylab("Relative abundance") +
  scale_y_continuous(expand = c(0,0))+theme(strip.background = element_rect(fill="gray85"))+
  theme(panel.margin = unit(0.3, "lines")) +
  theme(axis.text.x=element_text(angle=30,hjust=1,vjust=0.5, size = 9, face = "bold"), 
        strip.text = element_text(size = 14, face = "bold"))+
  theme(legend.key.size = unit(0.7, "cm"), legend.key.width = unit(0.7,"cm"), 
        legend.text = element_text(size = 12))+
  guides(fill=guide_legend(ncol=1))
txp_dna_algae2_kinzig2021
ggsave("txp_dna_protists_kinzig2021.tiff", units="in", width=15, height=10, dpi=300, compression = 'lzw')
###########

# Compute alpha diversity indices of the diatom communities
# *******************************************************

# Get help on the diversity() function
#library(vegan)
?diversity

N0 <- rowSums(dna_spe_trans > 0)         # Species richness
N0a <- specnumber(dna_spe_trans)     # Species richness (alternate)
H <- diversity(dna_spe_trans)       # Shannon entropy (base e)
Hb2 <- diversity(dna_spe_trans, base = 2) # Shannon entropy (base 2)
N1 <- exp(H)                   # Shannon diversity (number of abundant species) (base e)
N1b2 <- 2^Hb2                 # Shannon diversity (base 2)
N2 <- diversity(dna_spe_trans, "inv")    # Simpson diversity (number of dominant species)
J <- H/log(N0)                 # Pielou evenness
E10 <- N1/N0                   # Shannon evenness (Hill's ratio)
E20 <- N2/N0                   # Simpson evenness (Hill's ratio)
(div <- data.frame(N0, N0a, H, Hb2, N1, N1b2, N2, E10, E20, J))
# make div dataset rownames as first column
div1 <- tibble::rownames_to_column(div, "site_name")
# make also metadata1 dataset rownames as first column
meta2 <- tibble::rownames_to_column(meta1, "site_name")
# assigning new column names to some columns in div dataframe
colnames(div1)[2] <- "Richness"
colnames(div1)[4] <- "Shannon"
colnames(div1)[8] <- "Simpson"
colnames(div1)[11] <- "Pielou"
# select all columns except "method"
spe <- diat_dm %>% dplyr::select(-method)
range(div2$Richness)
range(div2$Shannon)
# remove column 3 in div2
div2 <- div2[,-3]
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
library(ggstatsplot)
# Comparison between species

# edit from here
x <- "surrounding"
cols <- 2:3 # the continuous dependent variables
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

library(patchwork)
plots <- list(p1,p2,p3,p4)

tiff("alpha_micro_dna.tiff", units="in", width=15, height=10, res=300)

wrap_plots(plots) +
  plot_layout(ncol = 2, nrow = 2) +
  plot_annotation(
    tag_levels = "a",
    tag_suffix = ")"
  )

wrap_plots(p1,p2,p3,p4,
           ncol = 2,
           nrow = 2,
           byrow = NULL,
           widths = NULL,
           heights = NULL,
           guides = NULL,
           tag_level = "a",
           tag_suffix = ")",
           design = NULL
)
dev.off()
dev.new()

?wrap_plots
patched <- plotlist1/plotlist2
patched + plot_annotation(tag_levels = 'a',
                          tag_suffix = ")")
#Dissolved Oxygen#######################
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

