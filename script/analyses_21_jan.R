#new analyses
#overarching aim: investigate recruitment and growth patterns in a Montane cloud forest of Taiwan
# in relation to species assemblages (BPR, sr), species functional distinctiveness and niche breadth

#hypotheses:
#H1A (growth in relation to sr and environmental covariates (light, ph, windwardness) at individual and subplot level)

#H1B (number or recruitment by species by subplot in relation to sr and environmental covariates)

#H2A (relationships between species functional distinctiveness and trait variation and light and soil resources availability)

#H2B (relationships between species variation in traits (separatedly, i also try pca but no results) and environmental covariates

#H3A ( relationship between growth with species trait variability and species distinctiveness at the species level)

#H3B ( relationship between total number of recruits per species with species trait variability and species distinctiveness at the species level)


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

growth$BA <- (growth$D_DBH*growth$D_DBH)*0.00007854 #returns BA in squared meters


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

#BA
#individual BA as sum keeping the cat info from the biggest branch
individual_BA <- growth %>%
   group_by(ind, subplot, species) %>%
   summarise(
      BA_growth = sum(BA, na.rm = TRUE),
      D_cat = D_cat[which.max(DBH_o)],
      .groups = "drop"
   )

#DBH increase by species by subplot and by mat or sap
species_subplot_BA_cat <- individual_BA %>%
   group_by(subplot, species, D_cat) %>%
   summarise(BA_growth = sum(BA_growth), .groups = 'drop')

#DBH increase by species by subplot
species_subplot_BA <- individual_BA %>%
   group_by(subplot, species) %>%
   summarise(BA_growth = sum(BA_growth), .groups = 'drop')

#DBH increase by subplot
subplot_BA <- growth %>%
   group_by(subplot) %>%
   summarise(total_D_BA = sum(BA, na.rm = TRUE))

#BA increase by species
species_BA <-growth %>%
   group_by(species) %>%
   summarise(total_D_BA = sum(BA, na.rm = TRUE))

#recruitment
# Sum recruitment per species across all subplots

#recruitment by species by subplot
species_subplot_recruitment <- recruit %>%
   group_by(subplot, species) %>%
   summarise(total_recruits = n(), .groups = 'drop')


#subplot
subplot_recruitment <- recruit %>%
   group_by(subplot) %>%
   summarise(total_recruits = n(), .groups = 'drop')

#species
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

#calculate subplot BA to add as covariate to consider structure
#reload raw data to calculate initial (t0) BA
DBH <- as.data.frame(read_excel("woody_species_records_retypeonly_Sole.20251219.xlsx", sheet = 1))
DBH$BA <- (DBH$DBH_o*DBH$DBH_o)*0.00007854 #returns BA in squared meters
names(DBH)[names(DBH) == "10 x 10"] <- 'subplot'

BA_subplot <-DBH%>%
   group_by(subplot) %>%
   summarise(
      BA = sum(BA, na.rm = TRUE),
      .groups = "drop"
   )

metadata$subplot_BA <- BA_subplot$BA

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

###IMPORTANT, remove species taken outside the 25 subplots?
#table to understand how many individuals by subplot by species
species_counts_df <- tree_trait_cv %>%
   group_by(subplot, species) %>%
   summarize(Count = n())

#remove that plot data
tree_trait_cv_clean <- tree_trait_cv[!tree_trait_cv$subplot %in% c("(3,8)", "(5,7)", "(7,1)",
                                                             "(9,0)","(5,0)","(8,1)","(7,0)",
                                                             "(4,3)","(7,3)","(9,9)","(6,7)"),]


#traits, logz traits and cv species average by plot
tree_trait_species_subplot <- tree_trait_cv_clean %>%
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
individual_DBH$BA_growth <- individual_BA$BA_growth
hypo1_df <- left_join(individual_DBH, metadata, by = "subplot")

hypo1_df <- hypo1_df %>%
   mutate(
      DBH_log = log1p(DBH_growth),          # log(x + 1) to handle zeros
      DBH_log = as.numeric(scale(DBH_log)), # scale after log
      BA_log = log1p(BA_growth),          # log(x + 1) to handle zeros
      BA_log = as.numeric(scale(BA_log)) # scale after log

   )

hypo1_df  <- hypo1_df %>%
   mutate(across(
      c(sr, rel_total_light, windwardness, pH, subplot_BA),
      ~ as.numeric(scale(.))
   ))

#plot to visualize growth and variables throughout the subplots
library(tidyr)
library(dplyr)
plot_df <- individual_DBH %>%
   mutate(subplot = gsub("[()]", "", subplot)) %>%
   separate(subplot, into = c("x", "y"), sep = ",")
subplot_summary <- plot_df %>%
   group_by(x, y) %>%
   summarize(
      total_growth = sum(DBH_growth, na.rm = TRUE),
      avg_growth = mean(DBH_growth, na.rm = TRUE),
      tree_count = n(), # Optional: to see how many trees are in the subplot
      .groups = "drop"
   )

#plot
library(ggplot2)
library(viridis) # For colorblind-friendly scales

d1 <- ggplot(subplot_summary, aes(x = x, y = y, fill = total_growth)) +
   geom_tile(color = "white") + # White border makes the grid clear
   geom_text(aes(label = round(total_growth, 1)), color = "white") +
   scale_fill_viridis_c() +
   coord_fixed() +
   theme_minimal() +
   labs(title = "Heatmap of Aggregate DBH Growth",
        fill = "Total Growth in cm")+theme(legend.position = "bottom")

d2 <- ggplot(subplot_summary, aes(x = x, y = y, fill = avg_growth)) +
   geom_tile(color = "white") + # White border makes the grid clear
   geom_text(aes(label = round(avg_growth, 1)), color = "white") +
   scale_fill_viridis_c() +
   coord_fixed() +
   theme_minimal() +
   labs(title = "Heatmap of Aggregate DBH Growth",
        fill = "Average Growth in cm")+theme(legend.position = "bottom")

#sr
metadata_plot <- metadata%>%
   mutate(subplot = gsub("[()]", "", subplot)) %>%
   separate(subplot, into = c("x", "y"), sep = ",")

s <- ggplot(metadata_plot, aes(x = x, y = y, fill = sr)) +
   geom_tile(color = "white") + # White border makes the grid clear
   geom_text(aes(label = round(sr, 1)), color = "white") +
   scale_fill_viridis_c(option = "magma") +
   coord_fixed() +
   theme_minimal() +
   labs(title = "Heatmap of SR",
        fill = "number of species") +theme(legend.position = "bottom")

s1 <- ggplot(metadata_plot, aes(x = x, y = y, fill = subplot_BA)) +
   geom_tile(color = "white") + # White border makes the grid clear
   geom_text(aes(label = round(subplot_BA, 1)), color = "white") +
   scale_fill_viridis_c(option = "magma") +
   coord_fixed() +
   theme_minimal() +
   labs(title = "Heatmap of Basal Area",
        fill = "number of species") +theme(legend.position = "bottom")

p <- ggplot(metadata_plot, aes(x = x, y = y, fill = soil_pH)) +
   geom_tile(color = "white") + # White border makes the grid clear
   geom_text(aes(label = round(soil_pH, 1)), color = "white") +
   scale_fill_viridis_c(option = "magma") +
   coord_fixed() +
   theme_minimal() +
   labs(title = "Heatmap of soil pH",
        fill = "number of species") +theme(legend.position = "bottom")

l <- ggplot(metadata_plot, aes(x = x, y = y, fill = rel_total_light)) +
   geom_tile(color = "white") + # White border makes the grid clear
   geom_text(aes(label = round(rel_total_light, 1)), color = "white") +
   scale_fill_viridis_c(option = "magma") +
   coord_fixed() +
   theme_minimal() +
   labs(title = "Heatmap of relative total light",
        fill = "number of species") +theme(legend.position = "bottom")

w <- ggplot(metadata_plot, aes(x = x, y = y, fill = windwardness)) +
   geom_tile(color = "white") + # White border makes the grid clear
   geom_text(aes(label = round(windwardness, 1)), color = "white") +
   scale_fill_viridis_c(option = "magma") +
   coord_fixed() +
   theme_minimal() +
   labs(title = "Heatmap of windwardness",
        fill = "number of species") +theme(legend.position = "bottom")

#distinciveness now
dist_summary <- distinctiveness %>%
   group_by(x, y) %>%
   summarize(
      total_Di = sum(Di, na.rm = TRUE),
      tree_count = n(), # Optional: to see how many trees are in the subplot
      .groups = "drop"
   )
dist_summary$x <- as.character(dist_summary$x)
dist_summary$y <- as.character(dist_summary$y)

di <- ggplot(dist_summary, aes(x = x, y = y, fill = total_Di)) +
   geom_tile(color = "white") + # White border makes the grid clear
   geom_text(aes(label = round(total_Di, 1)), color = "white") +
   scale_fill_viridis_c(option = "turbo") +
   coord_fixed() +
   theme_minimal() +
   labs(title = "Heatmap of Di",
        fill = "Distinctiveness")+theme(legend.position = "bottom")

di

library(ggpubr)
ggarrange(s,d,di, ncol =3)


#fit model for saplings and DBH adding plot BA
library(MuMIn)
library(performance)
library(sjPlot)
# mod_h1a_sap <- lmer(
#    DBH_log ~ sr + rel_total_light + windwardness + soil_pH + subplot_BA +
#       (1 | species) + (1 | subplot), data = hypo1_df[hypo1_df$D_cat == "sap",])
#
# summary(mod_h1a_sap)
# r2_nakagawa(mod_h1a_sap)
# plot_model(mod_h1a_sap)
#
# library(DHARMa)
# res <- simulateResiduals(mod_h1a_sap)
# plot(res)
#
# library(car)
# vif(mod_h1a_sap)


###USE TWEEDIE
#a specialized Generalized Linear Model (GLM) designed for data that is non-negative,
#heavily right-skewed, and contains a cluster of exact zero values

# library(glmmTMB)
#
# # Use RAW DBH_growth (NOT scaled), but scale your predictors
# mod_h1a_sap_tweedie <- glmmTMB(DBH_growth ~ sr + rel_total_light +
#                           windwardness + soil_pH + subplot_BA +
#                           (1 | species) + (1 | subplot),
#                        data = hypo1_df[hypo1_df$D_cat == "sap", ],
#                        family = tweedie(link = "log"))
#
# summary(mod_h1a_sap_tweedie)
# r2_nakagawa(mod_h1a_sap_tweedie)
# plot_model(mod_h1a_sap_tweedie)
#
# library(DHARMa)
# # For glmmTMB, it's best to use at least 1000 simulations
# res_tweedie <- simulateResiduals(mod_h1a_sap_tweedie, n = 1000)
# plot(res_tweedie)
#
#
# #BA
# mod_h1a_sap_tweedie_BA <- glmmTMB(BA_growth ~ sr + rel_total_light +
#                                   windwardness + soil_pH + subplot_BA +
#                                   (1 | species) + (1 | subplot),
#                                data = hypo1_df[hypo1_df$D_cat == "sap", ],
#                                family = tweedie(link = "log"))
#
#
# summary(mod_h1a_sap_tweedie_BA)
# r2_nakagawa(mod_h1a_sap_tweedie_BA)
# plot_model(mod_h1a_sap_tweedie_BA)
# plot(simulateResiduals(mod_h1a_sap_tweedie_BA))
#
# #plot both
# library(broom.mixed)
#
# # Extract results
# res_dbh <- tidy(mod_h1a_sap_tweedie, conf.int = TRUE) %>% mutate(model = "DBH Growth")
# res_ba  <- tidy(mod_h1a_sap_tweedie_BA, conf.int = TRUE) %>% mutate(model = "BA Growth")
#
# # Combine and clean (remove random effects and intercepts)
# all_res <- bind_rows(res_dbh, res_ba) %>%
#    filter(effect == "fixed", term != "(Intercept)") %>%
#    mutate(is_significant = ifelse(p.value < 0.05, "Yes", "No"))
#
# ggplot(all_res, aes(x = estimate, y = term, color = is_significant)) +
#    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
#    geom_point(size = 3) +
#    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
#    facet_wrap(~model) +
#    # Custom colors: Blue for significant, Grey for not
#    scale_color_manual(values = c("Yes" = "steelblue", "No" = "grey70")) +
#    theme_minimal() +
#    labs(title = "Comparing Predictors of Growth",
#         subtitle = "Significant results (p < 0.05) highlighted in blue",
#         x = "Estimate (Log Scale)",
#         y = "Predictor Variable",
#         color = "Significant?")
#
#
# #fit model for mature
# mod_h1a_mat_tweedie <- glmmTMB(DBH_growth ~ sr + rel_total_light +
#                           windwardness + soil_pH + subplot_BA +
#                           (1 | species) + (1 | subplot),
#                        data = hypo1_df[hypo1_df$D_cat == "mat", ],
#                        family = tweedie(link = "log"))
#
# summary(mod_h1a_mat_tweedie)
# r2_nakagawa(mod_h1a_mat_tweedie)
# plot_model(mod_h1a_mat_tweedie)
# plot(simulateResiduals(mod_h1a_mat_tweedie))
#
# #BA
# mod_h1a_mat_tweedie_BA <- glmmTMB(BA_growth ~ sr + rel_total_light +
#                                   windwardness + soil_pH + subplot_BA +
#                                   (1 | species) + (1 | subplot),
#                                data = hypo1_df[hypo1_df$D_cat == "mat", ],
#                                family = tweedie(link = "log"))
# summary(mod_h1a_mat_tweedie_BA)
# r2_nakagawa(mod_h1a_mat_tweedie_BA)
# plot_model(mod_h1a_mat_tweedie_BA)
# plot(simulateResiduals(mod_h1a_mat_tweedie_BA))


#### try with spatial structure
# 1. Clean the subplot strings into X and Y
library(tidyr)
hypo1_df <- hypo1_df %>%
   mutate(subplot_xy = gsub("[()]", "", subplot)) %>%
   separate(subplot_xy, into = c("x", "y"), sep = ",", convert = TRUE)

# 2. Create the pos factor and the group_id column
library(glmmTMB)
hypo1_df$pos <- numFactor(hypo1_df$x, hypo1_df$y)
hypo1_df$group_id <- factor(1) # This is the missing step!

# 3. Filter your data for the model
sap_data <- hypo1_df[hypo1_df$D_cat == "sap", ]

# 4. Run the model using the new dataframe
mod_h1a_sap_tweedie_spatial <- glmmTMB(DBH_growth ~ sr + rel_total_light + windwardness +
                          soil_pH + subplot_BA +
                          (1 | species) + (1 | subplot) +
                          exp(pos + 0 | group_id),
                       data = sap_data,
                       family = tweedie(link = "log"))



summary(mod_h1a_sap_tweedie_spatial)
plot_model(mod_h1a_sap_tweedie_spatial)
plot(simulateResiduals(mod_h1a_sap_tweedie_spatial))

#mature
mat_data <- hypo1_df[hypo1_df$D_cat == "mat", ]

# 4. mature
# mod_h1a_mat_tweedie_spatial <- glmmTMB(DBH_growth ~ sr + rel_total_light + windwardness +
#                                           soil_pH + subplot_BA +
#                                           (1 | species) + (1 | subplot) +
#                                           exp(pos + 0 | group_id),
#                                        data = mat_data,
#                                        family = tweedie(link = "log"))
#
# summary(mod_h1a_mat_tweedie_spatial)
# #r2(mod_h1a_mat_tweedie_spatial)
# plot_model(mod_h1a_mat_tweedie_spatial)
# plot(simulateResiduals(mod_h1a_mat_tweedie_spatial))

#simplify since we have NA in AIC and no spatial corrlations seemingly
mod_h1a_mat_tweedie_spatial <- glmmTMB(DBH_growth ~ sr + rel_total_light + windwardness +
                            soil_pH + subplot_BA +
                            (1 | species) + (1 | subplot),
                         data = mat_data,
                         family = tweedie(link = "log"))

summary(mod_h1a_mat_tweedie_spatial)
plot_model(mod_h1a_mat_tweedie_spatial)
plot(simulateResiduals(mod_h1a_mat_tweedie_spatial))

library(ggeffects)
# Predict effects for the two main drivers
eff_mat_ph <- ggpredict(mod_h1a_mat_tweedie_spatial, terms = "soil_pH")
eff_mat_light <- ggpredict(mod_h1a_mat_tweedie_spatial, terms = "rel_total_light")

# Plot pH
plot(eff_mat_ph) +
   labs(title = "Mature Tree Growth vs Soil pH", x = "Soil pH (Scaled)", y = "DBH Growth")
# Plot light
plot(eff_mat_light) +
   labs(title = "Mature Tree Growth vs Light availability", x = "Relative light available (Scaled)", y = "DBH Growth")


# #different try
# library(mgcv)

# # A Generalized Additive Mixed Model (GAMM)
# # 's(x, y)' creates a flexible spatial surface
# library(mgcv)
# sap_data$species <- as.factor(sap_data$species)
# mod_h1a_sap_tweedie_gam <- gam(DBH_growth ~ sr + rel_total_light + windwardness +
#                   soil_pH + subplot_BA +
#                   s(species, bs="re") +
#                   s(x, y), # The spatial 'map'
#                data = sap_data,
#                family = tw(link="log"))
#
# summary(mod_h1a_sap_tweedie_gam)
# r2(mod_h1a_sap_tweedie_gam)
# plot_model(mod_h1a_sap_tweedie_gam)
# plot(simulateResiduals(mod_h1a_sap_tweedie_gam))
#
# mat_data$species <- as.factor(mat_data$species)
# mod_h1a_mat_tweedie_gam <- gam(DBH_growth ~ sr + rel_total_light + windwardness +
#                                   soil_pH + subplot_BA +
#                                   s(species, bs="re") +
#                                   s(x, y), # The spatial 'map'
#                                data = mat_data,
#                                family = tw(link="log"))
#
# summary(mod_h1a_mat_tweedie_gam)
# r2(mod_h1a_mat_tweedie_gam)
# plot_model(mod_h1a_mat_tweedie_gam)
# plot(simulateResiduals(mod_h1a_mat_tweedie_gam))

#
# #2. subplot level
# hypo1b_df <- left_join(subplot_DBH, metadata, by = "subplot")
#
# hypo1b_df <- hypo1b_df %>%
#    mutate(
#       total_DBH_log = log1p(total_D_DBH),          # log(x + 1) to handle zeros
#       total_DBH_log = as.numeric(scale(total_DBH_log)) # scale after log
#    )
#
# hypo1b_df  <- hypo1b_df %>%
#    mutate(across(
#       c(total_DBH_log, total_D_DBH, sr, rel_total_light, windwardness, pH),
#       ~ as.numeric(scale(.))
#    ))
#
# #fit model
# mod_h1b <- lm(total_DBH_log ~ sr + rel_total_light + windwardness + soil_pH, data = hypo1b_df)

#summary(mod_h1b)
#plot(mod_h1b)

#now with sapling and mature separatedly
#2. subplot level
# hypo1b_cat_df <- left_join(species_subplot_DBH_cat, metadata, by = "subplot")
#
# hypo1b_cat_df <- hypo1b_cat_df %>%
#    mutate(
#       DBH_log = log1p(DBH_growth),          # log(x + 1) to handle zeros
#       DBH_log = as.numeric(scale(DBH_log)) # scale after log
#    )
#
# hypo1b_cat_df  <- hypo1b_cat_df %>%
#    mutate(across(
#       c(DBH_log, DBH_growth, sr, rel_total_light, windwardness, pH),
#       ~ as.numeric(scale(.))
#    ))
#
# #fit model
# mod_h1b_sap <- lm(DBH_log ~ sr + rel_total_light + windwardness + soil_pH, data = hypo1b_cat_df[hypo1b_cat_df$D_cat == "sap",])
#
# summary(mod_h1b_sap)
#mod_h1b_mat <- lm(DBH_log ~ sr + rel_total_light + windwardness + soil_pH, data = hypo1b_cat_df[hypo1b_cat_df$D_cat == "mat",])

#summary(mod_h1b_mat)

#b) recruitment at a plot level
hypo1bb_df <- left_join(subplot_recruitment, metadata, by = "subplot")
hypo1bb_df  <- hypo1bb_df %>%
   mutate(across(
      c(sr, rel_total_light, windwardness, pH),
      ~ as.numeric(scale(.))
   ))

hypo1bb_df <- hypo1bb_df %>%
   mutate(subplot_clean = gsub("[()]", "", subplot)) %>%
   separate(subplot_clean, into = c("x", "y"), sep = ",", convert = TRUE)

# 2. Create the spatial factor
hypo1bb_df$pos <- numFactor(hypo1bb_df$x, hypo1bb_df$y)
hypo1bb_df$group_id <- factor(1) # Still used to define the single study area
# # 2. Fit the Spatial Poisson Model
# mod_recruitment <- glmmTMB(total_recruits ~ sr + rel_total_light +
#                               soil_pH + subplot_BA + windwardness +
#                               exp(pos + 0 | group_id),
#                            data = hypo1bb_df,
#                            family = poisson(link = "log"))
#
# summary(mod_recruitment)

# Switch family to nbinom2
mod_rec_final <- glmmTMB(total_recruits ~ sr + rel_total_light + soil_pH +
                            subplot_BA + windwardness,
                         data = hypo1bb_df,
                         family = nbinom2(link = "log"))

summary(mod_rec_final)

library(DHARMa)
testDispersion(simulateResiduals(mod_rec_final))

library(ggeffects)
# Predict recruits based on SR
eff_rec <- ggpredict(mod_rec_final, terms = "sr")

# Plot it
plot(eff_rec) +
   labs(title = "Effect of Species Richness on Recruitment",
        x = "Species Richness (Scaled)",
        y = "Predicted Number of Recruits") +
   theme_minimal()

# 1. Generate predicted values for Species Richness
# This uses the 'mod_rec_final' (nbinom2) you just ran
eff_rec_sr <- ggpredict(mod_rec_final, terms = "sr")

# 2. Create a high-quality visualization
plot_recruitment <- ggplot(eff_rec_sr, aes(x = x, y = predicted)) +
   # Add the confidence interval ribbon
   geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "#2c7fb8", alpha = 0.15) +
   # Add the main prediction line
   geom_line(color = "#2c7fb8", size = 1.2) +
   # Formatting
   labs(
      title = "Impact of Diversity on Forest Recruitment",
      subtitle = "Predicted counts with 95% Confidence Intervals",
      x = "Species Richness (Standardized)",
      y = "Predicted Number of Recruits per Plot"
   ) +
   theme_minimal(base_size = 14) +
   theme(
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
   )

# plot
print(plot_recruitment)
plot_model(mod_rec_final)
plot(simulateResiduals(mod_rec_final))


#r2
# 1. Fit the Null Model (Intercept only)
mod_null <- glmmTMB(total_recruits ~ 1,
                    data = hypo1bb_df,
                    family = nbinom2(link = "log"))

# 2. Extract values
L_full <- as.numeric(logLik(mod_rec_final))
L_null <- as.numeric(logLik(mod_null))
n <- nobs(mod_rec_final)

# 3. Nagelkerke Formula
r2_nagelkerke <- (1 - exp(-2/n * (L_full - L_null))) / (1 - exp(2/n * L_null))

print(paste("Nagelkerke Pseudo-R2:", round(r2_nagelkerke, 3)))


####plot all effects in one
library(broom.mixed)

# 1. Extract Sapling results (Spatial Tweedie)
res_sap <- tidy(mod_h1a_sap_tweedie_spatial, conf.int = TRUE) %>%
   mutate(Model = "Sapling Growth (n=3336)")

# 2. Extract Mature results (Standard Tweedie)
res_mat <- tidy(mod_h1a_mat_tweedie_spatial, conf.int = TRUE) %>%
   mutate(Model = "Mature Growth (n=913)")

# 3. Extract Recruitment results (Negative Binomial)
res_rec <- tidy(mod_rec_final, conf.int = TRUE) %>%
   mutate(Model = "Recruitment (n=73 Plots)")

# 4. Combine all, filter out Intercepts, and flag Significance
all_results <- bind_rows(res_sap, res_mat, res_rec) %>%
   filter(effect == "fixed", term != "(Intercept)") %>%
   mutate(Significant = ifelse(p.value < 0.05, "Yes", "No"))

# 5. Create the 3-panel Plot
ggplot(all_results, aes(x = estimate, y = term, color = Significant)) +
   geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
   geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.3, size = 0.8) +
   geom_point(size = 3) +
   # Facet by Model, keeping the 3 lifecycle stages separate
   facet_wrap(~factor(Model, levels=c("Recruitment (n=73 Plots)",
                                      "Sapling Growth (n=3336)",
                                      "Mature Growth (n=913)")),
              scales = "free_x") +
   scale_color_manual(values = c("Yes" = "#2171b5", "No" = "#bdbdbd")) +
   theme_minimal() +
   labs(title = "Forest Dynamics Across Life Stages",
        subtitle = "Significant drivers (p < 0.05) highlighted in blue",
        x = "Coefficient Estimate (Log Scale)",
        y = "") +
   theme(legend.position = "bottom",
         strip.text = element_text(size = 15, face = "bold"),
         axis.text.y = element_text(size = 14),
         axis.text.x = element_text(size = 13))


#tables
library(broom.mixed)
library(purrr)

# Function to tidy and add model labels
get_tidy <- function(model, label) {
   tidy(model, conf.int = TRUE) %>%
      mutate(model_name = label) %>%
      filter(effect == "fixed") # Focus only on fixed effects
}

# Combine them
final_table <- bind_rows(
   get_tidy(mod_rec_final, "Recruitment"),
   get_tidy(mod_h1a_sap_tweedie_spatial, "Sapling"),
   get_tidy(mod_h1a_mat_tweedie_spatial, "Mature")
)

# Export to Excel-ready CSV
setwd("~/Documents/git/lalashan/results")
write.csv(final_table, "Growth_rec_SR_Forest_Models.csv", row.names = FALSE)


# #fit model
# mod_h1bb <- glm(total_recruits ~ sr + rel_total_light + windwardness + soil_pH, family = "poisson", data = hypo1bb_df)
#
# summary(mod_h1bb)
#
# # Call:
# #    glm(formula = total_recruits ~ sr + rel_total_light + windwardness +
# #           soil_pH, family = "poisson", data = hypo1bb_df)
# #
# # Coefficients:
# #    Estimate Std. Error z value Pr(>|z|)
# # (Intercept)      0.95025    1.09105   0.871   0.3838
# # sr               0.49270    0.07000   7.038 1.95e-12 ***
# #    rel_total_light  0.10429    0.05929   1.759   0.0786 .
# # windwardness     0.13619    0.07139   1.908   0.0564 .
# # soil_pH          0.11219    0.30118   0.372   0.7095
# # ---
# #    Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# #
# # (Dispersion parameter for poisson family taken to be 1)
# #
# # Null deviance: 313.84  on 72  degrees of freedom
# # Residual deviance: 133.97  on 68  degrees of freedom
# # AIC: 366.08
# #
# # Number of Fisher Scoring iterations: 5
#
# plot(mod_h1bb)

####H2####
#2. we expect a positive relationships between species functional distinctiveness and trait variation
#and light and soil resources availability

#2. method:
#a)linear mixed model relating species distinctiveness at individual level with sr + light + soil resources + windwardness as covariates + subplot as random effect
hypo2_df <- merge(tree_trait_cv, metadata, by = "subplot")
hypo2_df$col <- paste(hypo2_df$subplot,hypo2_df$species, sep = "")
hypo2_df <- merge(hypo2_df, distinctiveness[!distinctiveness$Di == 0,], by = "col")
hypo2_df <- hypo2_df[,-62]
#rename species col species
names(hypo2_df)[names(hypo2_df) == 'species.x'] <- 'species'

hypo2_df  <- hypo2_df %>%
   mutate(across(
      c(rel_total_light, windwardness, soil_pH, subplot_BA),
      ~ as.numeric(scale(.))
   ))

hypo2_df <- hypo2_df %>%
   mutate(subplot_clean = gsub("[()]", "", subplot)) %>%
   separate(subplot_clean, into = c("x", "y"), sep = ",", convert = TRUE)

# 2. Create the spatial factor
hypo2_df$pos <- numFactor(hypo2_df$x, hypo2_df$y)
hypo2_df$group_id <- factor(1) # Still used to define the single study area

#fit model for distinctiveness
distinctiveness_gamma <- glmmTMB(Di ~ rel_total_light + windwardness +
                        soil_pH + subplot_BA +
                        (1 | species) +
                        exp(pos + 0 | group_id),
                     data = hypo2_df,
                     family = Gamma(link = "log"))

summary(distinctiveness_gamma)
plot(simulate_residuals(distinctiveness_gamma))

# Create the plot

####plot all effects in one
library(broom.mixed)

# 1. Extract Sapling results (Spatial Tweedie)
res_di <- tidy(distinctiveness_gamma, conf.int = TRUE)%>%
   mutate(Significant = ifelse(p.value < 0.05, "Yes", "No"))%>%
   filter(effect == "fixed", term != "(Intercept)") %>%
   # Optional: Clean up term names for a prettier y-axis
   mutate(term = gsub("_", " ", term))

# 5. Create Plot
ggplot(res_di, aes(x = estimate, y = reorder(term, estimate), color = Significant)) +
   # Add a vertical line at 0 (the 'no effect' line)
   geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
   # Confidence intervals
   geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                  height = 0.2, size = 0.8) +
   # Point estimates
   geom_point(size = 3.5) +
   # Custom colors: Blue for significant, Gray for non-significant
   scale_color_manual(values = c("Yes" = "#2171b5", "No" = "#bdbdbd")) +
   theme_minimal() +
   labs(x = "Coefficient Estimate (Log Scale)",
        y = NULL,
        title = "Fixed Effects Results",
        subtitle = "Error bars represent 95% Confidence Intervals") +
   theme(legend.position = "bottom",
         strip.text = element_text(size = 15, face = "bold"),
         axis.text.y = element_text(size = 14),
         axis.text.x = element_text(size = 13),
         panel.grid.minor = element_blank())

#b)model relating species traits variation at individual level with sr + light + soil resources + windwardness as covariates + subplot as random effect
library(tidyr)
# 1. Reshape data
long_df <- hypo2_df %>%
   pivot_longer(cols = c(LA_CV, SLA_CV, LDMC_CV, LT_CV, Chl_CV),
                names_to = "Trait",
                values_to = "CV")

# 2. Model with interaction
mod_all_traits <- glmmTMB(CV ~ Trait * rel_total_light + soil_pH + windwardness + subplot_BA +
                             (1 | species) + (1 | subplot),
                          data = long_df,
                          family = tweedie(link = "log"),
                          start = list(psi = log(1.5 - 1) - log(2 - 1.5)),
                          map = list(psi = factor(NA)))

summary(mod_all_traits)
plot(simulate_residuals(mod_all_traits))

# Calculate predicted marginal effects
# 'terms = c("rel_total_light", "Trait")' tells it to plot light on X and group by Trait
pred_cv <- predict_response(mod_all_traits, terms = c("rel_total_light", "Trait"))

# 2. Create the Interaction Plot
ggplot(pred_cv, aes(x = x, y = predicted, color = group, fill = group)) +
   geom_line(size = 1.2) +
   geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1, color = NA) +
   scale_color_brewer(palette = "Set1") +
   scale_fill_brewer(palette = "Set1") +
   theme_minimal() +
   labs(
      title = "Effect of Light on Trait Variability (CV)",
      subtitle = "Lines show predicted slopes; Shaded areas are 95% CI",
      x = "Relative Total Light (Standardized)",
      y = "Predicted CV",
      color = "Leaf Trait",
      fill = "Leaf Trait"
   ) +
   theme(
      legend.position = "right",
      axis.title = element_text(face = "bold"),
      plot.title = element_text(size = 16)
   )


#facet plot
ggplot() +
   # Raw data points in the background
   geom_jitter(data = long_df,
               aes(x = rel_total_light, y = CV, color = Trait),
               alpha = 0.15, size = 0.8, width = 0.02) +
   # Model prediction lines
   geom_line(data = pred_cv,
             aes(x = x, y = predicted, color = group),
             size = 1.2) +
   # Confidence intervals
   geom_ribbon(data = pred_cv,
               aes(x = x, ymin = conf.low, ymax = conf.high, fill = group),
               alpha = 0.2) +
   # The Facet Command
   # "scales = 'free_y'" lets each trait have its own Y-axis range
   facet_wrap(~group, scales = "free_y", ncol = 3) +
   scale_color_brewer(palette = "Dark2") +
   scale_fill_brewer(palette = "Dark2") +
   theme_minimal() +
   labs(
      title = "Trait Variability (CV) across Light Gradient",
      subtitle = "Faceted by Leaf Trait (Note different Y-axis scales)",
      x = "Relative Total Light (Standardized)",
      y = "Trait CV",
      color = "Trait",
      fill = "Trait"
   ) +
   theme(
      legend.position = "none", # Legend is redundant with facet titles
      strip.text = element_text(face = "bold", size = 12),
      panel.spacing = unit(1.5, "lines")
   )

# install.packages("sjPlot")
library(sjPlot)

tab_model(mod_all_traits,
          show.re.var = TRUE,
          show.stat = TRUE,
          show.zeroinf = FALSE,
          p.style = "numeric_stars",
          title = "Mixed-Effects Model Summary (Tweedie Family)",
          dv.labels = "Trait CV (Log Link)")


#3. hypothesis: we expect a positive relationship between growth and recruitment with
# trait variations and distinctiveness at the species by suplot level

#individual growth merge with individual traits
hypo3_df <- left_join(tree_trait_cv, hypo1_df, by = "ind")
#add distinctiveness
hypo3_df$col <- paste(hypo3_df$subplot.x,hypo3_df$species.x, sep = "")
#join with distinctiveness
hypo3_df <- merge(hypo3_df, distinctiveness, by = "col", all = F)
#Create the pos factor and the group_id column

hypo3_df <- hypo3_df %>%
   mutate(subplot_clean = gsub("[()]", "", subplot.x)) %>%
   separate(subplot_clean, into = c("x", "y"), sep = ",", convert = TRUE)

hypo3_df  <- hypo3_df %>%
   mutate(across(
      c(SLA_CV, LA_CV, LDMC_CV, LT_CV, Chl_CV, Di, rel_total_light, windwardness, soil_pH, subplot_BA),
      ~ as.numeric(scale(.))
   ))

hypo3_df$pos <- numFactor(hypo3_df$x, hypo3_df$y)
hypo3_df$group_id <- factor(1) # This is the missing step!

library(glmmTMB)
library(broom.mixed)
library(dplyr)
library(purrr)

# Define your formula once to stay DRY (Don't Repeat Yourself)
model_formula <- DBH_growth ~ SLA_CV + LA_CV + LDMC_CV + LT_CV + Chl_CV + Di +
   rel_total_light + windwardness + soil_pH + subplot_BA +
   (1 | species) + (1 | subplot.x) + exp(pos + 0 | group_id)

# Split data by life stage and run models in a list
model_list <- hypo3_df %>%
   split(.$D_cat) %>%
   map(~ glmmTMB(model_formula, data = .x, family = tweedie(link = "log")))

# Extract tidied results from both models into one dataframe
plot_data_auto <- model_list %>%
   map_dfr(tidy, conf.int = TRUE, .id = "Model") %>%
   filter(effect == "fixed", term != "(Intercept)")

# Optional: Rename categories for the plot
plot_data_auto <- plot_data_auto %>%
   mutate(Model = dplyr::recode(Model, "sap" = "Sapling", "mat" = "Mature"),
          term = dplyr::recode(term, "rel_total_light" = "Light",
                        "windwardness" = "Windwardness",
                        "soil_pH" = "pH",
                        "subplot_BA" = "Subplot_BA"))

# Add a significance column to help the viewer's eye
plot_data_auto <- plot_data_auto %>%
   mutate(Significant = if_else(conf.low > 0 | conf.high < 0, "Yes", "No"))

# The Plot
ggplot(plot_data_auto, aes(x = estimate, y = reorder(term, estimate), color = Model)) +
   # Add a vertical reference line at 0
   geom_vline(xintercept = 0, linetype = "dashed", color = "gray60", size = 0.6) +

   # Add the error bars (Confidence Intervals)
   geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                  height = 0.3, size = 0.9, alpha = 0.7,
                  position = position_dodge(width = 0.6)) +

   # Add the point estimates
   # We use a different shape for significant results to make it 'colorblind friendly'
   geom_point(aes(shape = Significant), size = 3.5,
              position = position_dodge(width = 0.6)) +

   # Custom Styling
   scale_color_manual(values = c("Sapling" = "#2c7fb8", "Mature" = "#d95f02")) +
   scale_shape_manual(values = c("Yes" = 19, "No" = 1)) + # Solid for Sig, Hollow for Non-Sig

   theme_minimal(base_size = 14) +
   labs(title = "Drivers of DBH Growth by Life Stage",
        subtitle = "Points represent log-betas; error bars are 95% CIs",
        x = "Standardized Coefficient (Log Scale)",
        y = NULL,
        color = "Life Stage",
        shape = "Significant (p < 0.05)") +

   theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(face = "bold", color = "black"),
      plot.title = element_text(face = "bold", size = 16),
      strip.background = element_rect(fill = "gray95", color = NA)
   )




###RECRUITS
#now same but with number or recruits per species CV and by species distinctiveness
#generate species level traits db
hypo3_df$ind <- as.factor(hypo3_df$ind)
hypo3_df$key <- paste(hypo3_df$species, hypo3_df$subplot.x)

species_subplot_recruitment$key <- paste(species_subplot_recruitment$species, species_subplot_recruitment$subplot)

hypo3b_df <- merge(species_subplot_recruitment, hypo3_df, by = "key", all = T)
hypo3b_df$total_recruits[is.na(hypo3b_df$total_recruits)] <- 0

#model
mod_h3_simplified <- glmmTMB(
   total_recruits ~ scale(SLA_CV) + scale(Di) + scale(rel_total_light) +
      scale(windwardness) + scale(soil_pH) + scale(subplot_BA) +
      (1 | species.x) + (1 | subplot.x),
   data = hypo3b_df,
   ziformula = ~0,           # Disabled because -18.86 suggests it's not needed
   family = nbinom2
)

summary(mod_h3_simplified)
r2(mod_h3_simplified)

library(MASS)

# Fit the Negative Binomial GLM
mod_final <- glm.nb(total_recruits ~ scale(SLA_CV) + scale(Di) +
                       scale(rel_total_light) + scale(windwardness) +
                       scale(soil_pH) + scale(subplot_BA),
                    data = hypo3b_df)

summary(mod_final)
r2(mod_final)

library(DHARMa)
plot(simulateResiduals(mod_final))

#plot
# Load necessary libraries
library(sjPlot)
library(ggplot2)

# Create the plot
plot_model(mod_final,
           transform = NULL,        # Keeps estimates on the log-link scale (standardized)
           show.values = TRUE,      # Shows the estimate value on the plot
           value.offset = .3,       # Offsets the text slightly
           vline.color = "red",     # Adds a red line at 0 (null effect)
           title = "Drivers of Seedling Recruitment (Negative Binomial GLM)",
           axis.labels = c("Subplot Basal Area", "Soil pH", "Windwardness",
                           "Rel. Total Light", "Distinctiveness (Di)", "SLA CV"),
           width = 0.2) +           # Adjusts thickness of the error bars
   theme_minimal() +
   labs(y = "Standardized Coefficients (Log-scale)")

#####SEM part
#4. # and we expected that species distinctiveness and trait variation are mediating factors
#of how tree species richness affects tree growth and recruitment (with environmental covariates as above)

#3. method
#a)1 SEM connecting all components for growth
#b)1 SEM connecting all components for recruitment

