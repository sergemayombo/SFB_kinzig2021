# CRC_RESIST
# Data analysis
# Benthotorch data measured during 2021 sampling efforts
# Tuesday 15th June 2023
# Serge Mayombo
# AG Phykologie
# University of Duisburg-Essen

# Read data into R
library(tidyverse)

data <- read.csv("kinzig2021_benthotorch.csv")
View(data)
str(data)

head(data, 10) # First five lines
tail(data, 10) # Last two lines
glimpse(data) # A more thorough summary
names(data) # THe names of the columns

summary(data)
str(data)
# read environmental data

meta_K <- read.csv("kinzig2021_fieldData.csv")
# remove observation 18
meta_K <- meta_K[-18,]
# combine dataframes
data_K <- left_join(meta_K, data_Kinzig)

##########################################
# Multiple comparisons
# install.packages("ggstatsplot")
library(ggstatsplot)
# Comparison between species
# Kinzig

data_K$surrounding <- as.factor(data_K$surrounding)
data_K$Diatoms <- as.numeric(data_K$Diatoms)
data_K$Cyano <- as.numeric(data_K$Cyano)
data_K$Green.Algae <- as.numeric(data_K$Green.Algae)
data_K$Total.Conc. <- as.numeric(data_K$Total.Conc.)
str(data_K)


# edit from here
x <- "surrounding"
cols <- c(18,19,21,20) # the 4 continuous dependent variables
type <- "parametric" # given the large number of observations, we use the parametric version
paired <- FALSE # FALSE for independent samples, TRUE for paired samples
# edit until here

# edit at your own risk
plotlist <-
  purrr::pmap(
    .l = list(
      data = list(as_tibble(data_K)),
      x = x,
      y = as.list(colnames(data_K)[cols]),
      plot.type = "box", # for boxplot
      type = type, # parametric or nonparametric
      pairwise.comparisons = TRUE, # to run post-hoc tests if more than 2 groups
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

library(ggpubr)
bentho_kinzig2021 <- ggarrange(plotlist[[1]],plotlist[[2]], plotlist[[3]],plotlist[[4]], 
                               labels = c("a)", "b)", "c)", "d)"),
                               ncol = 2, nrow = 2)

bentho_kinzig2021
ggsave("bentho_kinzig2021.tiff", units="in", width=11, height=7,  dpi=300, compression = 'lzw')
#####################################################
# test for all other continual variables with lm
# using ggstatsplot
library(BayesFactor)
ggstatsplot::ggscatterstats(
  data = data_K, 
  x = water_level_at_pole, 
  y = Cyano,
  title = "",
  messages = FALSE
)
y
ggstatsplot::ggscatterstats(
  data = data_K, 
  x = conductivity, 
  y = Green.Algae,
  title = "",
  messages = FALSE
)
y
ggstatsplot::ggscatterstats(
  data = data_K, 
  x = conductivity, 
  y = Diatoms,
  title = "",
  messages = FALSE
)
#################################################
y

# computing welch t-test
# photosynthetic biomass
library(rstatix)
bentho.test1 <- data_K %>%
  t_test(Total.Conc. ~ surrounding) %>%
  add_significance()
bentho.test1
# Effect size: Cohen's d formula for Welch t-test
data_K %>% cohens_d(Total.Conc. ~ surrounding, var.equal = FALSE)

bentho.test2 <- data_K %>%
  t_test(Diatoms ~ surrounding) %>%
  add_significance()
bentho.test2
# Effect size: Cohen's d formula for Welch t-test
data_K %>% cohens_d(Diatoms ~ surrounding, var.equal = FALSE)

bentho.test3 <- data_K %>%
  t_test(Cyano ~ surrounding) %>%
  add_significance()
bentho.test3
# Effect size: Cohen's d formula for Welch t-test
data_K %>% cohens_d(Cyano ~ surrounding, var.equal = FALSE)

bentho.test4 <- data_K %>%
  t_test(Green.Algae ~ surrounding) %>%
  add_significance()
bentho.test4
# Effect size: Cohen's d formula for Welch t-test
data_K %>% cohens_d(Green.Algae ~ surrounding, var.equal = FALSE)

#############
#############
#install.packages("patchwork")
library(patchwork)

dev.new()
guilds + plotlist
#dev.off()
dev.new()
wrap_plots(plotlist) +
  plot_layout(ncol = 2) +
  plot_annotation(
    tag_levels = "A",
    tag_suffix = ")"
  )


#######################################################
# Assuming your data is in a dataframe called df
# Here's a mock dataframe as an example:
set.seed(123)  # for reproducibility
#df <- data.frame(
# site_type = c(rep("rural", 100), rep("urban", 100)),
#chlorophyll_biomass = c(rnorm(100, 5, 1), rnorm(100, 4.5, 1))  # example data
#)
df <- data_K
# Perform a t-test to test for differences
t_result <- t.test(Diatoms ~ surrounding, data = df)

# Calculate the percentage difference
rural_mean <- mean(df$Diatoms[df$surrounding == "rural"])
urban_mean <- mean(df$Diatoms[df$surrounding == "urban"])
percent_diff <- ((rural_mean - urban_mean) / urban_mean) * 100

# Display the results
cat(sprintf("The rural sites had on average %.2f%% higher Diatoms (p = %.3f) than the urban ones.\n", 
            percent_diff, t_result$p.value))



# computing welch t-test
# Richness
library(rstatix)
bentho.test1 <- data_K %>%
  t_test(Total.Conc. ~ surrounding) %>%
  add_significance()
bentho.test1
# Effect size: Cohen's d formula for Welch t-test
data_K %>% cohens_d(Total.Conc. ~ surrounding, var.equal = FALSE)

bentho.test2 <- data_K %>%
  t_test(Diatoms ~ surrounding) %>%
  add_significance()
bentho.test2
# Effect size: Cohen's d formula for Welch t-test
data_K %>% cohens_d(Diatoms ~ surrounding, var.equal = FALSE)

bentho.test3 <- data_K %>%
  t_test(Cyano ~ surrounding) %>%
  add_significance()
bentho.test3
# Effect size: Cohen's d formula for Welch t-test
data_K %>% cohens_d(Cyano ~ surrounding, var.equal = FALSE)

bentho.test4 <- data_K %>%
  t_test(Green.Algae ~ surrounding) %>%
  add_significance()
bentho.test4
# Effect size: Cohen's d formula for Welch t-test
data_K %>% cohens_d(Green.Algae ~ surrounding, var.equal = FALSE)

