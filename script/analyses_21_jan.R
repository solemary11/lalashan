#new analyses
#overarching aim: investigate recruitment and growth patterns in a Montane cloud forest of Taiwan
# in relation to species assemblages (BPR, sr), species functional distinctiveness and niche breadth

#1. hypothesis: we expect a positive relationship between sr and
   #a) growth at: 1. individual level, 2. subplot level
   #b) recruitment at a plot level

#1. method:
    #a)linear mixed model relating DBH_growth at individual level with sr + light + soil resources + windwardness as covariates + species and subplot as random effect
    #b)linear model relating DBH_growth at subplot level with sr + light + soil resources + windwardness as covariates
    #c)glm with poisson family relating total recruitment per subplot to sr + light + soil resources + windwardness as covariates

#2. we expect a positive relationships between species functional distinctiveness and trait variation
#and light and soil resources availability

#2. method:
    #a)linear mixed model relating species distinctiveness at individual level with sr + light + soil resources + windwardness as covariates + subplot as random effect
    #b)single linear mixed model relating species traits variation at individual level with sr + light + soil resources + windwardness as covariates + subplot as random effect

#3. hypothesis: we expect a positive relationship between growth and recruitment with
#mean species distinctiveness and trait variations or plasticity at the species level
# and we expected that species distinctiveness and trait variation are mediating factors
#of how tree species richness affects tree growth and recruitment (with environmental covariates as above)

#3. method
    #a)1 SEM connecting all components for growth
    #b)1 SEM connecting all components for recruitment

######LOADING AND CLEANING######
#########libraries and data

library(readxl)
library(funrar)
library(dplyr)
library(tidyverse)
library(purrr)
library(lme4)
library(lmerTest)

#import lalashan data
setwd("~/Documents/git/lalashan/all data compiled")
#tree comm
tree <- as.data.frame(read_excel("LFDP_master_20221027.xlsx", sheet = 1)) #importance values 0-1, could be treated as relative abundances
#tree traits
tree_trait <- read_excel("LFDP_master_20221027.xlsx", sheet = 6)
tree_trait <- tree_trait[,c(1,2,3,4,5,6,7,8)]
colnames(tree_trait) <- c("subplot", "species", "ind","LA","Lth","Chl","SLA","LDMC")

leaf_trait <- read_excel("leaf_traits_qnap_16-1-26.xlsx", sheet = 1)
leaf_trait <-leaf_trait[,c(1,2,3,4,5,6,7,8,9,10,11,12,13)]
colnames(leaf_trait) <- c("subplot","ind","species","leaf_type","leaf","fresh_weight",
                          "dry_weight","Lt1","Lt2","Lt3","Lt4","LA","Chl")

#dbh resample from me:
growth <- as.data.frame(read_excel("woody_species_records_retypeonly_Sole.20251219.xlsx", sheet = 1))
unique(growth$Note)# remove all obs with notes!
growth <- growth[is.na(growth$Note),]
growth <- growth[,c(1,2,3,7,8,9)]
colnames(growth) <- c("subplot","ind","branch","species","DBH_o","DBH_n")
#new recuitments from luke:
recruit <- read_excel("retyped_by_luke_new_recruits.xlsx", sheet = 5)
#move old trees from recruit to growth db
old_tree <- recruit[recruit$Note == "old tree",]
old_tree <- old_tree[rowSums(is.na(old_tree)) != ncol(old_tree),]
old_tree <- old_tree[,c(1,2,3,7,8,9)]
colnames(old_tree) <- c("subplot","ind","branch","species","DBH_o","DBH_n")
growth <- rbind(growth, old_tree)
#find duplicated ind+branch names to remove in case
dup <- duplicated(growth[,2:3])
growth[,dup]#no duplicates, good!

#further clean growth data
growth$D_DBH <- growth$DBH_n-growth$DBH_o #get difference, negative values are dead branches or measurements errors, remove only dead and keep errors as 0
growth$dead[growth$DBH_n == 0] <- 'yes'

#remove dead
growth <- growth %>%
   filter(is.na(dead))

#change other negative values to 0
growth$D_DBH[!growth$D_DBH >= 0] <- 0

#also put a category of diameter for further analysis to discriminate mature individuals and saplings that have higher growth rates!
growth$D_cat <- ""
growth$D_cat[growth$DBH_o >= 10] <- "mat" #10 cm threshold
growth$D_cat[!growth$DBH_o >= 10] <- "sap"

#keep other records as recruitment
recruit <- recruit[!recruit$Note %in% c("old tree"),]
recruit <- recruit[,c(1,2,3,7,8,9)]
colnames(recruit) <- c("subplot","ind","branch","species","DBH_o","DBH_n")
#add old tree to growth df

#plot data - environmental
topo <- read_excel("LFDP_master_20221027.xlsx", sheet = 9)
soil_subp <- read_excel("LFDP_master_20221027.xlsx", sheet = 10)
soil_analysis <- read_excel("LFDP_master_20221027.xlsx", sheet = 11)
light <- read_excel("LFDP_master_20221027.xlsx", sheet = 13)#use rel total light

#merge all plot metdatata
metadata <- merge(topo, soil_subp, by = "subplot", all =T)
metadata <- merge(metadata, soil_analysis, by = "subplot", all = T)
metadata <- merge(metadata, light, by = "subplot")


######DATA PREPARATION######

#individual DBH as sum keeping the cat info from the biggest branch
individual_DBH <- growth %>%
   group_by(ind, subplot, species) %>%
   summarise(
      DBH_growth = sum(D_DBH, na.rm = TRUE),
      D_cat = D_cat[which.max(DBH_o)],
      .groups = "drop"
   )

#DBH increase by species by subplot and by mat or sap
species_subplot_DBH_cat <- individual_DBH %>%
   group_by(subplot, species, D_cat) %>%
   summarise(DBH_growth = sum(DBH_growth), .groups = 'drop')

#DBH increase by species by subplot
species_subplot_DBH <- individual_DBH %>%
   group_by(subplot, species) %>%
   summarise(DBH_growth = sum(DBH_growth), .groups = 'drop')

#DBH increase by subplot
subplot_DBH <- growth %>%
   group_by(subplot) %>%
   summarise(total_D_DBH = sum(D_DBH, na.rm = TRUE))

#DBH increase by species
species_DBH <-growth %>%
   group_by(species) %>%
   summarise(total_D_DBH = sum(D_DBH, na.rm = TRUE))

#recruitment
# Sum recruitment per species across all subplots
subplot_recruitment <- recruit %>%
   group_by(subplot) %>%
   summarise(total_recruits = n(), .groups = 'drop')

species_recruitment <- recruit %>%
   group_by(species) %>%
   summarise(total_recruits = n(), .groups = 'drop')

#species richness by subplot
#Calculate species richness by subplot
tree_PA <- tree
tree_PA[!tree_PA == 0] <- 1
sr <- as.data.frame(rowSums(tree_PA[,-1]))
sr$subplot <- tree$subplot
colnames(sr) <- c("sr","subplot")

#add to metadata sr
metadata <- left_join(metadata, sr, by = "subplot")

#traits

#add SLA to leaf_ trait in mm/mg
leaf_trait$SLA <- (leaf_trait$LA/(leaf_trait$dry_weight*1000))
#add LDMC mg/g
leaf_trait$LDMC <- ((leaf_trait$dry_weight*1000)/leaf_trait$fresh_weight)
#add mean LT
leaf_trait$LT <- rowMeans(leaf_trait[, c(8:11)], na.rm = TRUE)
leaf_trait <- leaf_trait[,c(1,2,3,4,5,12:16)]

#now compute individual coefficients of variations

cv_individuals <- leaf_trait %>%
   group_by(ind, species) %>%
   summarise(
      across(
         .cols = LA:LT,                  # numeric traits
         .fns  = ~sd(.x, na.rm = TRUE) / mean(.x, na.rm = TRUE),  # CV
         .names = "{.col}_CV"
      ),
      n = n(),   # number of individuals per species
      .groups = "drop"
   )

#merge with individual traits
tree_trait_cv <- merge(tree_trait, cv_individuals, by = "ind")

#then standardize traits to log10 and scale
log_transform_scale <- function(x) {
   if (all(is.na(x))) return(x)

   min_x <- min(x, na.rm = TRUE)
   if (min_x <= 0) {
      x <- x + abs(min_x) + 1
   }

   as.numeric(scale(log10(x)))
}

#add to traits also scaled values
traits_to_transform <- c("LA", "Lth", "Chl", "SLA", "LDMC")
tree_trait_cv <- tree_trait_cv %>%
   mutate(across(
      all_of(traits_to_transform),
      log_transform_scale,
      .names = "{.col}_logz"
   ))
tree_trait_cv <- tree_trait_cv[,-c(9,15)]
names(tree_trait_cv)[names(tree_trait_cv) == 'species.x'] <- 'species'

#traits, logz traits and cv species average by plot
tree_trait_species_subplot <- tree_trait_cv %>%
   group_by(subplot, species) %>%
   summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
   ungroup()
tree_trait_species_subplot <- tree_trait_species_subplot[,-3]

######DISTINCTIVENESS######
#distinctiveness by subplot and by species
trait_cols <- c("LA_logz", "Lth_logz", "Chl_logz", "SLA_logz", "LDMC_logz")
# split by subplot and compute distance matrices

dist_matrices <- tree_trait_species_subplot %>%
   split(.$subplot) %>%
   lapply(function(sub_df) {

      # build trait table
      trait_table <- sub_df %>%
         select(species, all_of(trait_cols)) %>%
         mutate(across(-species, ~ as.numeric(.))) %>%
         as.data.frame()

      rownames(trait_table) <- trait_table$species
      trait_table$species <- NULL

      #compute distance matrix
      compute_dist_matrix(trait_table[,-1], metric = 'euclidean')
   })

#remove elements that have 1x1 (only one species)
dist_matrices_clean <- dist_matrices[!sapply(dist_matrices, function(x) {
   is.matrix(x) && all(dim(x) == c(1, 1))
})]

#split also community data into subplots as elements of list
species_cols <- setdiff(names(tree), c("subplot"))

#split into a list of matrices, one per subplot
subplot_matrices <- tree%>%
   split(.$subplot) %>%
   lapply(function(x) {
      mat <- make_relative(as.matrix(x[, species_cols, drop = FALSE]))#turn also into relative abundances
      rownames(mat) <- x$subplot
      mat
   })

#######distinctiveness with funrar package######

#transform each matrix into a daframe with columns needed for the distinctiveness function
subplot_stacked <- lapply(subplot_matrices, function(mat) {

   dat <- matrix_to_stack(mat, "value", "site", "species")

   dat$site <- as.character(dat$site)
   dat$species <- as.character(dat$species)

   dat
})

#apply this function to all elements (subplot communities and subplot functional distances)
#keep only subplots present in BOTH lists
#keep only elements of subplot_stacked that are in dist_matrices
subplot_stacked_matched <- subplot_stacked[names(dist_matrices_clean)]
#ensure the order matches dist_matrices exactly
subplot_stacked_matched <- subplot_stacked_matched[names(dist_matrices_clean)]

#try
d1 <- distinctiveness_stack(subplot_stacked_matched[[1]], "species","site","value",dist_matrices[[1]], relative = T)

#apply in parallel
result_list <- map2(
   subplot_stacked_matched,
   dist_matrices_clean,
   ~ distinctiveness_stack(.x, "species", "site", "value", .y, relative = TRUE)
)

#now restructure to plot
#expand the list and clean coordinates
df_list <- bind_rows(result_list) %>%
   mutate(
      #extract coords
      coords = str_match(site, "(\\d+)[^0-9]+(\\d+)"),
      x = as.numeric(coords[,2]),
      y = as.numeric(coords[,3])
   ) %>%
   select(-coords) #remove

distinctiveness <- df_list %>%
   complete(
      species,
      x = 0:9,
      y = 0:9,
      fill = list(Di = 0)
   )

#bind to species average traits and cv
tree_trait_species_subplot$col <- paste(tree_trait_species_subplot$subplot,tree_trait_species_subplot$species, sep = "")
distinctiveness$col <- paste("(",distinctiveness$x, ",",distinctiveness$y,")" , distinctiveness$species, sep = "")
tree_trait_species_subplot <- merge(tree_trait_species_subplot, distinctiveness[,c(6,7)], by = "col")
tree_trait_species_subplot <-tree_trait_species_subplot[,-1]

######ANALYSES######

####H1####
#1. hypothesis: we expect a positive relationship between sr and
#a) growth at:
#1. individual level
#individual_DBH, add sr, rel_light and soilPH, windwardness
hypo1_df <- left_join(individual_DBH, metadata, by = "subplot")

hypo1_df <- hypo1_df %>%
   mutate(
      DBH_log = log1p(DBH_growth),          # log(x + 1) to handle zeros
      DBH_log = as.numeric(scale(DBH_log)) # scale after log
   )

hypo1_df  <- hypo1_df %>%
   mutate(across(
      c(DBH_log, DBH_growth, sr, rel_total_light, windwardness, pH),
      ~ as.numeric(scale(.))
   ))

#fit model for saplings
mod_h1a_sap <- lmer(
   DBH_log ~ sr + rel_total_light + windwardness + soil_pH +
      (1 | species) + (1 | subplot), data = hypo1_df[hypo1_df$D_cat == "sap",])

summary(mod_h1a_sap)

# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: DBH_log ~ sr + rel_total_light + windwardness + soil_pH + (1 |
#     species) + (1 | subplot)
#    Data: hypo1_df[hypo1_df$D_cat == "sap", ]
#
# REML criterion at convergence: 7716.4
#
# Scaled residuals:
#     Min      1Q  Median      3Q     Max
# -2.3305 -0.6644 -0.2217  0.4693  7.6803
#
# Random effects:
#  Groups   Name        Variance Std.Dev.
#  subplot  (Intercept) 0.01380  0.1175
#  species  (Intercept) 0.04389  0.2095
#  Residual             0.56867  0.7541
# Number of obs: 3336, groups:  subplot, 100; species, 58
#
# Fixed effects:
#                   Estimate Std. Error         df t value Pr(>|t|)
# (Intercept)      -0.663688   0.325765 100.526085  -2.037  0.04425 *
# sr                0.089487   0.029314  74.824893   3.053  0.00314 **
# rel_total_light   0.003798   0.021088  75.202215   0.180  0.85756
# windwardness      0.060231   0.028038  77.768963   2.148  0.03481 *
# soil_pH           0.101268   0.088796  98.116843   1.140  0.25687
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Correlation of Fixed Effects:
#             (Intr) sr     rl_tt_ wndwrd
# sr           0.226
# rl_ttl_lght  0.148 -0.122
# windwardnss -0.009 -0.632 -0.145
# soil_pH     -0.992 -0.225 -0.151  0.017

#fit model for mature
mod_h1a_mat <- lmer(
   DBH_growth ~ sr + rel_total_light + windwardness + soil_pH +
      (1 | species) + (1 | subplot), data = hypo1_df[hypo1_df$D_cat == "mat",])

summary(mod_h1a_mat)

# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: DBH_growth ~ sr + rel_total_light + windwardness + soil_pH +
#     (1 | species) + (1 | subplot)
#    Data: hypo1_df[hypo1_df$D_cat == "mat", ]
#
# REML criterion at convergence: 3217.2
#
# Scaled residuals:
#     Min      1Q  Median      3Q     Max
# -1.6636 -0.5632 -0.2097  0.3333  9.2080
#
# Random effects:
#  Groups   Name        Variance Std.Dev.
#  subplot  (Intercept) 0.1032   0.3212
#  species  (Intercept) 0.2988   0.5466
#  Residual             1.8087   1.3449
# Number of obs: 913, groups:  subplot, 100; species, 36
#
# Fixed effects:
#                 Estimate Std. Error       df t value Pr(>|t|)
# (Intercept)     -1.44413    0.99513 89.56975  -1.451   0.1502
# sr              -0.02233    0.08801 73.76481  -0.254   0.8004
# rel_total_light  0.10450    0.06321 76.87443   1.653   0.1024
# windwardness    -0.04912    0.08469 75.46793  -0.580   0.5637
# soil_pH          0.51879    0.27071 86.50524   1.916   0.0586 .
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Correlation of Fixed Effects:
#             (Intr) sr     rl_tt_ wndwrd
# sr           0.205
# rl_ttl_lght  0.170 -0.147
# windwardnss  0.041 -0.590 -0.110
# soil_pH     -0.991 -0.201 -0.173 -0.026


#2. subplot level
hypo1b_df <- left_join(subplot_DBH, metadata, by = "subplot")

hypo1b_df <- hypo1b_df %>%
   mutate(
      total_DBH_log = log1p(total_D_DBH),          # log(x + 1) to handle zeros
      total_DBH_log = as.numeric(scale(total_DBH_log)) # scale after log
   )

hypo1b_df  <- hypo1b_df %>%
   mutate(across(
      c(total_DBH_log, total_D_DBH, sr, rel_total_light, windwardness, pH),
      ~ as.numeric(scale(.))
   ))

#fit model
mod_h1b <- lm(total_DBH_log ~ sr + rel_total_light + windwardness + soil_pH, data = hypo1b_df)

summary(mod_h1b)
# Call:
#    lm(formula = total_DBH_log ~ sr + rel_total_light + windwardness +
#          soil_pH, data = hypo1b_df)
#
# Residuals:
#    Min      1Q  Median      3Q     Max
# -1.6751 -0.2811  0.0497  0.3285  1.2318
#
# Coefficients:
#    Estimate Std. Error t value Pr(>|t|)
# (Intercept)      0.30976    1.00466   0.308    0.759
# sr               0.46724    0.08600   5.433 4.27e-07 ***
#    rel_total_light -0.05811    0.06797  -0.855    0.395
# windwardness     0.42006    0.08371   5.018 2.43e-06 ***
#    soil_pH         -0.08621    0.27908  -0.309    0.758
# ---
#    Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 0.6295 on 95 degrees of freedom
# Multiple R-squared:  0.6197,	Adjusted R-squared:  0.6037
# F-statistic: 38.71 on 4 and 95 DF,  p-value: < 2.2e-16

plot(mod_h1b)

#now with sapling and mature separatedly
#2. subplot level
hypo1b_cat_df <- left_join(species_subplot_DBH_cat, metadata, by = "subplot")

hypo1b_cat_df <- hypo1b_cat_df %>%
   mutate(
      DBH_log = log1p(DBH_growth),          # log(x + 1) to handle zeros
      DBH_log = as.numeric(scale(DBH_log)) # scale after log
   )

hypo1b_cat_df  <- hypo1b_cat_df %>%
   mutate(across(
      c(DBH_log, DBH_growth, sr, rel_total_light, windwardness, pH),
      ~ as.numeric(scale(.))
   ))

#fit model
mod_h1b_sap <- lm(DBH_log ~ sr + rel_total_light + windwardness + soil_pH, data = hypo1b_cat_df[hypo1b_cat_df$D_cat == "sap",])

summary(mod_h1b_sap)
# Call:
#    lm(formula = DBH_log ~ sr + rel_total_light + windwardness +
#          soil_pH, data = hypo1b_cat_df[hypo1b_cat_df$D_cat == "sap",
#          ])
#
# Residuals:
#    Min      1Q  Median      3Q     Max
# -1.5963 -0.6320 -0.2503  0.4211  4.1984
#
# Coefficients:
#    Estimate Std. Error t value Pr(>|t|)
# (Intercept)     -0.57853    0.46057  -1.256 0.209335
# sr               0.13139    0.03889   3.378 0.000755 ***
#    rel_total_light -0.02090    0.03013  -0.694 0.488110
# windwardness     0.15736    0.03796   4.146 3.65e-05 ***
#    soil_pH          0.11344    0.12719   0.892 0.372616
# ---
#    Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 0.9162 on 1114 degrees of freedom
# Multiple R-squared:  0.08181,	Adjusted R-squared:  0.07851
# F-statistic: 24.81 on 4 and 1114 DF,  p-value: < 2.2e-16

mod_h1b_mat <- lm(DBH_log ~ sr + rel_total_light + windwardness + soil_pH, data = hypo1b_cat_df[hypo1b_cat_df$D_cat == "mat",])

summary(mod_h1b_mat)
# Call:
#    lm(formula = DBH_log ~ sr + rel_total_light + windwardness +
#          soil_pH, data = hypo1b_cat_df[hypo1b_cat_df$D_cat == "mat",
#          ])
#
# Residuals:
#    Min      1Q  Median      3Q     Max
# -1.7021 -0.7977 -0.1346  0.6991  4.3454
#
# Coefficients:
#    Estimate Std. Error t value Pr(>|t|)
# (Intercept)      0.42384    0.73409   0.577   0.5640
# sr              -0.13171    0.06433  -2.047   0.0411 *
#    rel_total_light  0.05302    0.04839   1.096   0.2737
# windwardness     0.14878    0.06320   2.354   0.0190 *
#    soil_pH         -0.02194    0.20235  -0.108   0.9137
# ---
#    Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 1.01 on 505 degrees of freedom
# Multiple R-squared:  0.01493,	Adjusted R-squared:  0.00713
# F-statistic: 1.914 on 4 and 505 DF,  p-value: 0.1068

#b) recruitment at a plot level
hypo1bb_df <- left_join(subplot_recruitment, metadata, by = "subplot")

hypo1bb_df  <- hypo1bb_df %>%
   mutate(across(
      c(sr, rel_total_light, windwardness, pH),
      ~ as.numeric(scale(.))
   ))

#fit model
mod_h1bb <- glm(total_recruits ~ sr + rel_total_light + windwardness + soil_pH, family = "poisson", data = hypo1bb_df)

summary(mod_h1bb)

# Call:
#    glm(formula = total_recruits ~ sr + rel_total_light + windwardness +
#           soil_pH, family = "poisson", data = hypo1bb_df)
#
# Coefficients:
#    Estimate Std. Error z value Pr(>|z|)
# (Intercept)      0.95025    1.09105   0.871   0.3838
# sr               0.49270    0.07000   7.038 1.95e-12 ***
#    rel_total_light  0.10429    0.05929   1.759   0.0786 .
# windwardness     0.13619    0.07139   1.908   0.0564 .
# soil_pH          0.11219    0.30118   0.372   0.7095
# ---
#    Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# (Dispersion parameter for poisson family taken to be 1)
#
# Null deviance: 313.84  on 72  degrees of freedom
# Residual deviance: 133.97  on 68  degrees of freedom
# AIC: 366.08
#
# Number of Fisher Scoring iterations: 5

plot(mod_h1bb)

####H2####
#2. we expect a positive relationships between species functional distinctiveness and trait variation
#and light and soil resources availability

#2. method:
#a)linear mixed model relating species distinctiveness at individual level with sr + light + soil resources + windwardness as covariates + subplot as random effect
hypo2_df <- merge(tree_trait_cv, metadata, by = "subplot")
hypo2_df$col <- paste(hypo2_df$subplot,hypo2_df$species, sep = "")
hypo2_df <- merge(hypo2_df, distinctiveness[!distinctiveness$Di == 0,], by = "col")
hypo2_df <- hypo2_df[,-61]
#rename species col species
names(hypo2_df)[names(hypo2_df) == 'species.x'] <- 'species'

hypo2_df  <- hypo2_df %>%
   mutate(across(
      c(Di, rel_total_light, windwardness, pH),
      ~ as.numeric(scale(.))
   ))

#fit model for distinctiveness
mod_h2a <- lmer(
   Di ~ rel_total_light + windwardness + soil_pH +
      (1 | species) + (1 | subplot), data = hypo2_df)
summary(mod_h2a)

# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: Di ~ rel_total_light + windwardness + soil_pH + (1 | species) +      (1 | subplot)
#    Data: hypo2_df
#
# REML criterion at convergence: 1216.9
#
# Scaled residuals:
#     Min      1Q  Median      3Q     Max
# -3.1472 -0.6331 -0.0536  0.5195  3.0355
#
# Random effects:
#  Groups   Name        Variance Std.Dev.
#  species  (Intercept) 0.6987   0.8359
#  subplot  (Intercept) 1.1566   1.0755
#  Residual             0.2233   0.4726
# Number of obs: 690, groups:  species, 52; subplot, 32
#
# Fixed effects:
#                 Estimate Std. Error      df t value Pr(>|t|)
# (Intercept)       4.3144     2.9008 23.2348   1.487   0.1504
# rel_total_light  -0.4062     0.1961 22.8695  -2.072   0.0497 *
# windwardness      0.1985     0.2159 23.0494   0.919   0.3674
# soil_pH          -0.9881     0.7907 23.1260  -1.250   0.2239
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Correlation of Fixed Effects:
#             (Intr) rl_tt_ wndwrd
# rl_ttl_lght  0.222
# windwardnss  0.253  0.051
# soil_pH     -0.997 -0.215 -0.243

#b)single linear mixed model relating species traits variation at individual level with sr + light + soil resources + windwardness as covariates + subplot as random effect
hypo2_df  <- hypo2_df %>%
   mutate(across(
      c(LA_CV, Chl_CV, SLA_CV, LT_CV, LDMC_CV),
      ~ as.numeric(scale(.))
   ))

#fit one model for each traits CV
#LA
mod_h2b_LA <- lmer(
   LA_CV ~ rel_total_light + windwardness + soil_pH +
      (1 | species) + (1 | subplot), data = hypo2_df)
summary(mod_h2b_LA)

# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: LA_CV ~ rel_total_light + windwardness + soil_pH + (1 | species) +      (1 | subplot)
#    Data: hypo2_df
#
# REML criterion at convergence: 1951.7
#
# Scaled residuals:
#     Min      1Q  Median      3Q     Max
# -1.9990 -0.7324 -0.1211  0.6031  3.4236
#
# Random effects:
#  Groups   Name        Variance Std.Dev.
#  species  (Intercept) 0.06471  0.25437
#  subplot  (Intercept) 0.00642  0.08012
#  Residual             0.93211  0.96546
# Number of obs: 690, groups:  species, 52; subplot, 32
#
# Fixed effects:
#                   Estimate Std. Error         df t value Pr(>|t|)
# (Intercept)      0.4230408  0.6341785 21.8014081   0.667    0.512
# rel_total_light  0.0070628  0.0405239 19.8015046   0.174    0.863
# windwardness     0.0006351  0.0435998 15.7772420   0.015    0.989
# soil_pH         -0.0999597  0.1705186 20.9098949  -0.586    0.564
#
# Correlation of Fixed Effects:
#             (Intr) rl_tt_ wndwrd
# rl_ttl_lght  0.046
# windwardnss  0.309  0.038
# soil_pH     -0.995 -0.046 -0.312

#SLA
mod_h2b_SLA <- lmer(
   SLA_CV ~ rel_total_light + windwardness + soil_pH +
      (1 | species) + (1 | subplot), data = hypo2_df)
summary(mod_h2b_SLA)

# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: SLA_CV ~ rel_total_light + windwardness + soil_pH + (1 | species) +      (1 | subplot)
#    Data: hypo2_df
#
# REML criterion at convergence: 1921.6
#
# Scaled residuals:
#     Min      1Q  Median      3Q     Max
# -1.7853 -0.6493 -0.2213  0.4486  6.6256
#
# Random effects:
#  Groups   Name        Variance Std.Dev.
#  species  (Intercept) 0.15494  0.3936
#  subplot  (Intercept) 0.02068  0.1438
#  Residual             0.85640  0.9254
# Number of obs: 690, groups:  species, 52; subplot, 32
#
# Fixed effects:
#                 Estimate Std. Error       df t value Pr(>|t|)
# (Intercept)      0.75097    0.71153 23.36455   1.055   0.3020
# rel_total_light -0.03775    0.04571 20.71262  -0.826   0.4183
# windwardness     0.10941    0.04995 18.63483   2.190   0.0414 *
# soil_pH         -0.17921    0.19161 22.44856  -0.935   0.3596
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Correlation of Fixed Effects:
#             (Intr) rl_tt_ wndwrd
# rl_ttl_lght  0.078
# windwardnss  0.300  0.034
# soil_pH     -0.994 -0.078 -0.301

#LT
mod_h2b_LT <- lmer(
   LT_CV ~ rel_total_light + windwardness + soil_pH +
      (1 | species) + (1 | subplot), data = hypo2_df)
summary(mod_h2b_LT)

# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: LT_CV ~ rel_total_light + windwardness + soil_pH + (1 | species) +      (1 | subplot)
#    Data: hypo2_df
#
# REML criterion at convergence: 1910.1
#
# Scaled residuals:
#     Min      1Q  Median      3Q     Max
# -1.8080 -0.6402 -0.1882  0.3844  5.5282
#
# Random effects:
#  Groups   Name        Variance Std.Dev.
#  species  (Intercept) 0.13107  0.3620
#  subplot  (Intercept) 0.02185  0.1478
#  Residual             0.84641  0.9200
# Number of obs: 690, groups:  species, 52; subplot, 32
#
# Fixed effects:
#                 Estimate Std. Error       df t value Pr(>|t|)
# (Intercept)      0.31285    0.71463 24.20192   0.438    0.665
# rel_total_light -0.07523    0.04597 21.47717  -1.636    0.116
# windwardness     0.08241    0.05024 19.45375   1.640    0.117
# soil_pH         -0.06388    0.19255 23.31534  -0.332    0.743
#
# Correlation of Fixed Effects:
#             (Intr) rl_tt_ wndwrd
# rl_ttl_lght  0.081
# windwardnss  0.300  0.034
# soil_pH     -0.994 -0.081 -0.300

#LDMC
mod_h2b_LDMC <- lmer(
   LDMC_CV ~ rel_total_light + windwardness + soil_pH +
      (1 | species) + (1 | subplot), data = hypo2_df)
summary(mod_h2b_LDMC)

# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: LDMC_CV ~ rel_total_light + windwardness + soil_pH + (1 | species) +      (1 | subplot)
#    Data: hypo2_df
#
# REML criterion at convergence: 1954.9
#
# Scaled residuals:
#     Min      1Q  Median      3Q     Max
# -1.2715 -0.5695 -0.2429  0.1919  7.8344
#
# Random effects:
#  Groups   Name        Variance Std.Dev.
#  species  (Intercept) 0.05299  0.2302
#  subplot  (Intercept) 0.01431  0.1196
#  Residual             0.93616  0.9676
# Number of obs: 690, groups:  species, 52; subplot, 32
#
# Fixed effects:
#                  Estimate Std. Error        df t value Pr(>|t|)
# (Intercept)      0.011355   0.686357 21.043588   0.017    0.987
# rel_total_light -0.030011   0.044066 18.864859  -0.681    0.504
# windwardness     0.053957   0.047780 16.155325   1.129    0.275
# soil_pH         -0.003372   0.184863 20.312203  -0.018    0.986
#
# Correlation of Fixed Effects:
#             (Intr) rl_tt_ wndwrd
# rl_ttl_lght  0.066
# windwardnss  0.306  0.037
# soil_pH     -0.996 -0.066 -0.306

#Chl
mod_h2b_Chl <- lmer(
   Chl_CV ~ rel_total_light + windwardness + soil_pH +
      (1 | species) + (1 | subplot), data = hypo2_df)
summary(mod_h2b_Chl)

# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: Chl_CV ~ rel_total_light + windwardness + soil_pH + (1 | species) +      (1 | subplot)
#    Data: hypo2_df
#
# REML criterion at convergence: 1964
#
# Scaled residuals:
#     Min      1Q  Median      3Q     Max
# -1.1418 -0.5118 -0.1880  0.2875 17.9250
#
# Random effects:
#  Groups   Name        Variance Std.Dev.
#  species  (Intercept) 0.022172 0.14890
#  subplot  (Intercept) 0.005545 0.07446
#  Residual             0.970387 0.98508
# Number of obs: 690, groups:  species, 52; subplot, 32
#
# Fixed effects:
#                 Estimate Std. Error       df t value Pr(>|t|)
# (Intercept)      0.27700    0.63308 19.81587   0.438   0.6664
# rel_total_light -0.08633    0.04055 18.17517  -2.129   0.0472 *
# windwardness     0.01557    0.04347 14.21521   0.358   0.7255
# soil_pH         -0.07478    0.17031 19.09460  -0.439   0.6655
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Correlation of Fixed Effects:
#             (Intr) rl_tt_ wndwrd
# rl_ttl_lght  0.045
# windwardnss  0.315  0.040
# soil_pH     -0.997 -0.045 -0.316

#3. hypothesis: we expect a positive relationship between growth and recruitment with
# trait variations and distinctiveness at the species by suplot level

#individual growth merge with individual traits
hypo3_df <- left_join(tree_trait_cv, hypo1_df, by = "ind")
#add distinctiveness
hypo3_df$col <- paste(hypo3_df$subplot.x,hypo3_df$species.x, sep = "")
#join with distinctiveness
hypo3_df <- merge(hypo3_df, distinctiveness, by = "col", all = F)

#use CV directly
mod_h3b_sap <- lmer(
   DBH_log ~ SLA_CV + LA_CV + LDMC_CV + LT_CV + Chl_CV + Di + windwardness + soil_pH +
      (1 | species.x) + (1 | subplot.x), data = hypo3_df[hypo3_df$D_cat == "sap",])
summary(mod_h3b_sap)

# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: DBH_log ~ SLA_CV + LA_CV + LDMC_CV + LT_CV + Chl_CV + Di + windwardness +
#     soil_pH + (1 | species.x) + (1 | subplot.x)
#    Data: hypo3_df[hypo3_df$D_cat == "sap", ]
#
# REML criterion at convergence: 959.2
#
# Scaled residuals:
#     Min      1Q  Median      3Q     Max
# -2.1586 -0.6965 -0.1822  0.5103  3.5190
#
# Random effects:
#  Groups    Name        Variance Std.Dev.
#  species.x (Intercept) 0.04987  0.2233
#  subplot.x (Intercept) 0.02920  0.1709
#  Residual              0.53123  0.7289
# Number of obs: 420, groups:  species.x, 49; subplot.x, 33
#
# Fixed effects:
#               Estimate Std. Error        df t value Pr(>|t|)
# (Intercept)    0.13766    0.81050  31.91687   0.170   0.8662
# SLA_CV         0.83337    0.96249 410.63555   0.866   0.3871
# LA_CV          1.00325    0.42226 408.10024   2.376   0.0180 *
# LDMC_CV        0.06499    1.09040 405.55450   0.060   0.9525
# LT_CV          0.35999    0.87942 404.67359   0.409   0.6825
# Chl_CV        -2.04546    0.90458 408.01666  -2.261   0.0243 *
# Di            -0.28033    0.29108 222.53122  -0.963   0.3365
# windwardness   0.15310    0.06729  22.08347   2.275   0.0330 *
# soil_pH       -0.08224    0.20575  27.36261  -0.400   0.6925
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Correlation of Fixed Effects:
#             (Intr) SLA_CV LA_CV  LDMC_C LT_CV  Chl_CV Di     wndwrd
# SLA_CV      -0.068
# LA_CV       -0.047  0.033
# LDMC_CV      0.050 -0.558 -0.039
# LT_CV        0.006 -0.380 -0.235  0.148
# Chl_CV      -0.098 -0.163 -0.040 -0.205  0.030
# Di          -0.314 -0.024 -0.058 -0.069 -0.068  0.180
# windwardnss  0.358 -0.087  0.027  0.038 -0.014 -0.031 -0.131
# soil_pH     -0.973  0.057 -0.012 -0.043 -0.033  0.024  0.140 -0.312

library(performance)
r2(mod_h3b_sap)
# Conditional R2: 0.173
# Marginal R2: 0.050

mod_h3b_mat <- lmer(
   DBH_log ~ SLA_CV + LA_CV + LDMC_CV + LT_CV + Chl_CV + Di + windwardness + soil_pH +
      (1 | species.x) + (1 | subplot.x), data = hypo3_df[hypo3_df$D_cat == "mat",])
summary(mod_h3b_mat)

# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: DBH_log ~ SLA_CV + LA_CV + LDMC_CV + LT_CV + Chl_CV + Di + windwardness +
#     soil_pH + (1 | species.x) + (1 | subplot.x)
#    Data: hypo3_df[hypo3_df$D_cat == "mat", ]
#
# REML criterion at convergence: 416.8
#
# Scaled residuals:
#     Min      1Q  Median      3Q     Max
# -1.8906 -0.6518 -0.1223  0.5611  2.9940
#
# Random effects:
#  Groups    Name        Variance Std.Dev.
#  species.x (Intercept) 0.3438   0.5864
#  subplot.x (Intercept) 0.1474   0.3839
#  Residual              1.0368   1.0182
# Number of obs: 141, groups:  species.x, 26; subplot.x, 25
#
# Fixed effects:
#              Estimate Std. Error       df t value Pr(>|t|)
# (Intercept)    1.8441     2.1129  15.6632   0.873    0.396
# SLA_CV         2.1996     2.7419 129.7876   0.802    0.424
# LA_CV         -0.7579     1.1081 122.2081  -0.684    0.495
# LDMC_CV       -2.2911     4.3031 130.6293  -0.532    0.595
# LT_CV          2.0446     1.7496 131.2118   1.169    0.245
# Chl_CV        -0.6649     0.9071 126.1519  -0.733    0.465
# Di            -1.1974     1.0379  72.3883  -1.154    0.252
# windwardness  -0.1226     0.1631  13.5820  -0.752    0.465
# soil_pH       -0.2559     0.5343  14.2187  -0.479    0.639
#
# Correlation of Fixed Effects:
#             (Intr) SLA_CV LA_CV  LDMC_C LT_CV  Chl_CV Di     wndwrd
# SLA_CV      -0.122
# LA_CV       -0.179  0.121
# LDMC_CV     -0.012 -0.640 -0.106
# LT_CV        0.082 -0.499 -0.234  0.259
# Chl_CV      -0.011  0.092 -0.109 -0.131 -0.020
# Di          -0.248  0.171  0.092 -0.039 -0.267 -0.050
# windwardnss  0.319 -0.026 -0.048  0.089 -0.066 -0.007  0.024
# soil_pH     -0.962  0.051  0.096  0.018 -0.041 -0.003  0.017 -0.295

r2(mod_h3b_mat)
# Conditional R2: 0.351
# Marginal R2: 0.043

#now same but with number or recruits per species CV and by species distinctiveness
#generate species level traits db
hypo3b_df <- hypo3_df %>%
   group_by(species) %>%
   summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
   ungroup()

hypo3b_df <- left_join(species_recruitment, hypo3b_df, by = "species")

#model
mod_h3b_2 <- glm(
   total_recruits ~ SLA_CV + LA_CV + LDMC_CV + LT_CV + Chl_CV + Di, family = "poisson",data = hypo3b_df)
summary(mod_h3b_2)

#Call:
#   glm(formula = total_recruits ~ SLA_CV + LA_CV + LDMC_CV + LT_CV +
#          Chl_CV + Di, family = "poisson", data = hypo3b_df)

#Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
#(Intercept)   4.2896     0.4159  10.314  < 2e-16 ***
#   SLA_CV      -28.2653     4.8180  -5.867 4.45e-09 ***
#   LA_CV        -0.6410     1.8654  -0.344 0.731116
#LDMC_CV      34.9161     5.9284   5.890 3.87e-09 ***
#   LT_CV         4.3040     3.9178   1.099 0.271950
#Chl_CV       11.5567     3.0252   3.820 0.000133 ***
#   Di           -4.9418     0.6436  -7.679 1.61e-14 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#(Dispersion parameter for poisson family taken to be 1)

#Null deviance: 592.64  on 31  degrees of freedom
#Residual deviance: 424.48  on 25  degrees of freedom
#(5 observations deleted due to missingness)
#AIC: 546.73

#Number of Fisher Scoring iterations: 6

#####SEM part
#4. # and we expected that species distinctiveness and trait variation are mediating factors
#of how tree species richness affects tree growth and recruitment (with environmental covariates as above)

#3. method
#a)1 SEM connecting all components for growth
#b)1 SEM connecting all components for recruitment

