#########libraries and data

library(readxl)
library(funrar)
library(dplyr)
library(corrgram)
library(purrr)

#import lalashan data
setwd("~/Documents/git/lalashan/all data compiled")
#tree comm
tree <- as.data.frame(read_excel("LFDP_master_20221027.xlsx", sheet = 1)) #importance values 0-1, could be treated as relative abundances
#tree traits
tree_trait <- read_excel("LFDP_master_20221027.xlsx", sheet = 6)

#new recuitments from luke:
recruit <- read_excel("retyped_by_luke_new_recruits.xlsx", sheet = 5)
#dbh resample from me:
growth <- as.data.frame(read_excel("woody_species_records_retypeonly_Sole.20251219.xlsx", sheet = 1))

#plot data - environmental
topo <- read_excel("LFDP_master_20221027.xlsx", sheet = 9)
soil_subp <- read_excel("LFDP_master_20221027.xlsx", sheet = 10)
soil_analysis <- read_excel("LFDP_master_20221027.xlsx", sheet = 11)

#merge all plot metdatata
metadata <- merge(topo, soil_subp, by = "subplot", all =T)
metadata <- merge(metadata, soil_analysis, by = "subplot", all = T)

####calculate response variables first!!!
#recruitment and DBH increase by species

#recruitment
#remove old trees
unique(recruit$Note)
recruit <- recruit[!recruit$Note %in% c("old tree"),]
# Sum recruitment per species across all subplots
total_recruit <- recruit %>%
   group_by(Spe_Lat) %>%
   summarise(total_recruits = n(), .groups = 'drop')

#DBH increase
unique(growth$Note)# remove all obs with notes!

growth <- growth[is.na(growth$Note),]
growth$D_DBH <- growth$DBH_n-growth$DBH_o #negative values are dead branches, remove
growth <- growth[!growth$D_DBH < 0,]
#also put a treshold to already more mature individuals to not to count saplings that have higher growth rates!
growth <- growth[!growth$DBH_o < 10,] #use 10cm for now

#now summarize species avg increment in DBH
total_DBH_growth <- growth %>%
   group_by(`Species (Latin)`) %>%
   summarise(DBH_growth = mean(D_DBH), .groups = 'drop')

setdiff(total_DBH_growth$`Species (Latin)`, tree_trait$species)
setdiff(total_DBH_growth$`Species (Latin)`, tree_trait$species)
# [1] "Benthamidia japonica var. chinensis" "Chamaecyparis formosensis"
# [3] "Chamaecyparis obtusa var. formosana" "Ilex lonicerifolia var. matsudae"
# [5] "Ilex suzukii"                        "Machilus thunbergii"
# [7] "Pieris taiwanensis"                  "Rhododendron leptosanthum"
# [9] "Rhododendron pseudochrysanthum"      "Symplocos migoi"
# [11] "Tsuga chinensis var. formosana"      NA

setdiff(tree_trait$species, total_DBH_growth$`Species (Latin)`)
# [1] "Schima superba"      "Viburnum urceolatum"

###keep in mind that there are differences!

#####now calculate tree species functional distinctiveness and niche breadth local and plot wide


#1.first plot wide trait variations!!!

#average values for all individuals to compute functional distinctiveness
#(row names must be species names)
#group and calculate mean values per species by subplot
tree_trait_avg <- as.data.frame(tree_trait %>%
   group_by(species) %>%
   summarise(
      LA = mean(LA, na.rm = TRUE),
      Lth = mean(Lth, na.rm = TRUE),
      Chl = mean(Chl, na.rm = TRUE),
      SLA = mean(SLA, na.rm = TRUE),
      LDMC = mean(LDMC, na.rm = TRUE),
      Lsu = mean(Lsu, na.rm = TRUE),
      LCC = mean(LCC, na.rm = TRUE),
      LNC = mean(LNC, na.rm = TRUE)
   ))

#then standardize traits to log10 and scale
log_transform_scale <- function(x) {
   if (all(is.na(x))) return(x)

   min_x <- min(x, na.rm = TRUE)
   if (min_x <= 0) {
      x <- x + abs(min_x) + 1
   }

   as.numeric(scale(log10(x)))
}

tree_trait_logsc <- tree_trait_avg %>%
   mutate(
      across(
         -c(species),
         log_transform_scale,
         .names = "{.col}_log10scaled"
      )
   )

rownames(tree_trait_logsc) <- tree_trait_avg$species
rownames(tree) <- tree$subplot
#distance matrix
dist_matrix.p <- compute_dist_matrix(tree_trait_logsc[,c(10,11,12,13,14)], metric = "euclidean")

#distinctiveness
distinct_plot <- distinctiveness(as.matrix(tree[,-1]), dist_matrix.p, relative = TRUE)

t <- as.data.frame(t(distinct_plot))
t$species <- rownames(t)

#to long format
l <- t %>%
   pivot_longer(
      cols = -species,
      names_to = "subplot",
      values_to = "value"
   )

#get species distinciviness across all subplots average, min, max and sd
species_distinctiveness <- l %>%
   group_by(species) %>%
   summarise(
      mean = mean(value, na.rm = TRUE),
      min  = min(value, na.rm = TRUE),
      max  = max(value, na.rm = TRUE),
      sd   = sd(value, na.rm = TRUE),
      n    = sum(!is.na(value)),
      .groups = "drop"
   )


#now get species functional space measure
# #use mFD package
# library(mFD)
#
# # Assess the quality of spaces with different numbers of dimensions (axes)
# # This helps in choosing the optimal number of dimensions to retain
# fspace_quality <-quality.fspaces(
#    sp_dist = as.dist(dist_matrix.p),
#    maxdim_pcoa = 5)
#
# # Select the coordinates from the best functional space (e.g., 4 dimensions had the lowest 'mad' index)
# # Replace '4D' with the optimal number of dimensions for your data
# sp_faxes_coord <- fspace_quality$details_fspaces$sp_pc_coord
# sp_faxes_coord
#
# plot(fspace_quality$details_fspaces$pc_eigenvalues$Rel_corr_eig)
#
# # Plot the position of species in a 2D functional space (e.g., PC1 vs PC2)
# funct.space.plot(
#    sp_faxes_coord[, c(1, 2)], # Select the first two axes for a 2D plot
#    plot_vertices = TRUE,
#    plot_sp_nm = rownames(sp_faxes_coord),
#    nm_size = 1# You can add species names near points
# )

#calculate species coefficient of variation in traits

library(dplyr)

cv_summary <- tree_trait %>%
   group_by(species) %>%
   summarise(
      across(
         .cols = LA:DBH,                  # numeric traits
         .fns  = ~sd(.x, na.rm = TRUE) / mean(.x, na.rm = TRUE),  # CV
         .names = "{.col}_CV"
      ),
      n = n(),   # number of individuals per species
      .groups = "drop"
   )

###now relate both distinctiveness and CVs to recruitment and DBH increase
#merge in one df:

species_plot <- merge(cv_summary, species_distinctiveness, by = "species", all = T)
colnames(total_recruit) <- c("species", "recruit")
species_plot <- merge(species_plot, total_recruit, by = "species", all = T)
colnames(total_DBH_growth) <- c("species", "growth")
species_plot <- merge(species_plot, total_DBH_growth, by = "species", all = T)

species_plot$recruit[is.na(species_plot$recruit)] <- 0

# Example: recruit ~ functional distinctiveness or trait CVs
plot(recruit ~ mean, data = species_plot)

#recruit
glm_recruit <- glm(
   recruit ~ mean + LA_CV + SLA_CV + LDMC_CV + Chl_CV + Lth_CV,
   family = poisson(link = "log"),
   data = species_plot
)

summary(glm_recruit)
car::vif(glm_recruit)#not high
anova(glm_recruit, test = "Chisq")
# Calculate the pseudo R-squared manually
pseudo_r_squared_recruit <- 1 - (glm_recruit$deviance / glm_recruit$null.deviance) #0.105187....

#growth
glm_growth <- glm(
   growth ~ mean + LA_CV + SLA_CV + LDMC_CV + Chl_CV + Lth_CV,
   family = gaussian(link = "log"),
   data = species_plot
)

summary(glm_growth)
car::vif(glm_growth)#not high
pseudo_r_squared_growth <- 1 - (glm_growth$deviance / glm_growth$null.deviance) #0.137916....

