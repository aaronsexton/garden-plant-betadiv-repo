
# Repository code for "Wild plants drive biotic differentiation across urban gardens"
# Code originally archived on Jan. 15th, 2025
# by Aaron Sexton


## Code Outline:
# 1. NMDS & Vector analysis of community composition
# 2. Urbanization & Beta div regressions
# 3. Wild and Cult contribution to overall Beta div




# 1. NMDS & Vector analysis ----


library(tidyverse)
library(vegan)


### Full Comm ----

# Read in the wide df
full_wide <- read.csv("../Data/Plants/UBH_manuscript_files/full_wide_2020_2023.csv")
str(full_wide)

# Drop rows with no cover
full_wide <- filter(full_wide, total_cover > 0)

# Create a df grouped to the Garden level
full_wide_grpd <- group_by(full_wide, garden_yr)
full_wide_gyrlvl <- full_wide_grpd %>% 
  summarise(across(Amaranthus.sp.:Veronica.sp., ~ sum(.x, na.rm = T)))

# Do the same at the survey level
grpd_survey <- group_by(full_wide, survey)
full_wide_surveylvl <- grpd_survey %>% 
  summarise(across(Amaranthus.sp.:Veronica.sp., ~ sum(.x, na.rm = T)))

# Clean up the env
rm(full_wide_grpd, grpd_survey)


# Create distance matrices at both levels
gyr_dist     <- vegdist(decostand(full_wide_gyrlvl[2:1110], "norm"), "euclidean")
survey_dist  <- vegdist(decostand(full_wide_surveylvl[2:1110], "norm"), "euclidean")


# Run NMDS's

## Garden_yr level
gyr_nms <- metaMDS(gyr_dist, distance = "euclidean", k = 3, maxit = 999,
                   try = 5, trymax = 250, wascores = T)
plot(gyr_nms, "sites")   # Produces distance 
orditorp(gyr_nms, "sites")   # Gives points labels

## Survey level
survey_nms <- metaMDS(survey_dist, distance = "euclidean", k = 3, maxit = 999,
                      try = 5, trymax = 25, wascores = T)
plot(survey_nms, "sites")   # Produces distance 
orditorp(survey_nms, "sites")   # Gives points labels

# Check the histograms here
hist(rowSums(full_wide_surveylvl[2:1110]))
hist(vegan::diversity(full_wide_surveylvl[2:1110]))



# Create predictor df's so we can do vector analyses
preds <- full_wide[,1:8]

gry_preds <- full_wide_gyrlvl[,1]
preds_gyr <- preds[!duplicated(preds$garden_yr), ]
gry_preds <- left_join(gry_preds, preds_gyr)
rm(preds_gyr)
survey_preds <- full_wide_surveylvl[,1]
preds_survey <- preds[!duplicated(preds$survey), ]
survey_preds <- left_join(survey_preds, preds_survey)
rm(preds_survey)



str(gry_preds)
str(survey_preds)

gry_preds$status <- gry_preds$impervious_1000
gry_preds$status[gry_preds$status > 50] <- "Urban"
gry_preds$status[gry_preds$status < 50] <- "Rural"
gry_preds$status <- gsub("9.57", "Rural", as.character(gry_preds$status))
gry_preds$status2 <- paste(gry_preds$city, gry_preds$status)

survey_preds$status <- survey_preds$impervious_1000
survey_preds$status[survey_preds$status > 50] <- "Urban"
survey_preds$status[survey_preds$status < 50] <- "Rural"
survey_preds$status <- gsub("9.57", "Rural", as.character(survey_preds$status))
survey_preds$status2 <- paste(survey_preds$city, survey_preds$status)









# Run the vector analyses


## With City and Urban as two categorical terms
str(gry_preds)
gry_preds <- select(gry_preds, garden_yr, garden, city, status, everything())
str(survey_preds)
survey_preds <- select(survey_preds, survey, garden, city, status, everything())

en_gyr    <- envfit(gyr_nms, gry_preds[3:4], permutations = 999, na.rm = TRUE)
en_survey <- envfit(survey_nms, survey_preds[3:4], permutations = 999, na.rm = TRUE)

en_gyr
en_survey

plot(gyr_nms)
plot(en_gyr)
plot(survey_nms)
plot(en_survey)




# Now make it pretty for plotting




# Garden_yr Level
data.scores.gardenyr <- as.data.frame(gyr_nms$points)
gyr_nms$species
gyr_nms$nobj
data.scores.gardenyr$status <- gry_preds$status2
data.scores.gardenyr$city <- gry_preds$city
data.scores.gardenyr$Urb <- gry_preds$status


spp.rel2 <- (decostand(full_wide_gyrlvl[2:1110], "norm"))
spp_scrs2 <- 
  sppscores(gyr_nms) <- spp.rel2

nms_species <- as.data.frame(gyr_nms$species)
nms_species$Taxa <- rownames(nms_species)




# Survey Level
data.scores.survey <- as.data.frame(survey_nms$points)
survey_nms$species
survey_nms$nobj
data.scores.survey$status <- survey_preds$status2
data.scores.survey$city <- survey_preds$city
data.scores.survey$Urb <- survey_preds$status


spp.rel2s <- (decostand(full_wide_surveylvl[2:1110], "norm"))
spp_scrs2s <- 
  sppscores(survey_nms) <- spp.rel2s

nms_species_s <- as.data.frame(survey_nms$species)
nms_species_s$Taxa <- rownames(nms_species)



## ggplot it

# Survey-level City
b <- ggplot(data = data.scores.survey, aes(x = MDS1, y = MDS2)) +
  stat_ellipse(aes(color = city), size = 1.5) +
  geom_point(data = data.scores.survey, aes(size = MDS3, color = city), alpha = 0.6) +
  scale_color_manual(values = c("royalblue", "tomato3")) +
  labs(color = "City") +
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), 
        axis.text = element_text(size = 8), 
        legend.title = element_text(size = 26, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 18, colour = "grey30"))
b
# Survey-level City and Urban
c <- ggplot(data = data.scores.survey, aes(x = MDS1, y = MDS2)) +
  stat_ellipse(aes(color = status), size = 1.5) +
  geom_point(data = data.scores.survey, aes(size = MDS3, color = status), alpha = 0.6) +
  scale_color_manual(values = c("powderblue", "royalblue", "indianred2", "tomato4")) +
  labs(color = "City & Status") +
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), 
        axis.text = element_text(size = 8), 
        legend.title = element_text(size = 26, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 18, colour = "grey30"))
c
# Garden-level City and Urban
d <- ggplot(data = data.scores.gardenyr, aes(x = MDS1, y = MDS2)) +
  stat_ellipse(aes(color = status), size = 1.5) +
  geom_point(data = data.scores.gardenyr, aes(size = MDS3, color = status), alpha = 0.6) +
  scale_color_manual(values = c("powderblue", "royalblue", "indianred2", "tomato4")) +
  labs(color = "City & Status") +
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), 
        axis.text = element_text(size = 8), 
        legend.title = element_text(size = 26, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 18, colour = "grey30"))
d

# Garden-level City and Urban
a <- ggplot(data = data.scores.gardenyr, aes(x = MDS1, y = MDS2)) +
  stat_ellipse(aes(color = city), size = 1.5) +
  geom_point(data = data.scores.gardenyr, aes(size = MDS3, color = city), alpha = 0.6) +
  scale_color_manual(values = c("royalblue", "tomato3")) +
  labs(color = "City") +
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), 
        axis.text = element_text(size = 8), 
        legend.title = element_text(size = 26, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 18, colour = "grey30"))
a


# Put them all together

library(ggpubr)

ggarrange(a, d, b, c, 
          ncol = 2, nrow = 2)

# Just City
ggarrange(a, b, 
          nrow= 2)




# Can save the coordinates for easier re-plotting

# write.csv(data.scores.gardenyr, "../Data/Plants/UBH_manuscript_files/nms_datascores_gardenyr.csv", row.names = F)
# write.csv(data.scores.survey,   "../Data/Plants/UBH_manuscript_files/nms_datascores_survey.csv", row.names = F)







### Wild ----

# Read in the wide df
wild_wide <- read.csv("../Data/Plants/UBH_manuscript_files/full_wide_wild_2020_2023.csv")
str(wild_wide)

# Drop rows with no cover
wild_wide <- filter(wild_wide, wild_total_cover > 0)

# Create a df grouped to the Garden_yr level
wild_wide_grpd <- group_by(wild_wide, garden_yr)
wild_wide_gyrlvl <- wild_wide_grpd %>% 
  summarise(across(Amaranthus.sp.:Taraxacum.sp., ~ sum(.x, na.rm = T)))

# Do the same at the survey  level
grpd_survey <- group_by(wild_wide, survey)
wild_wide_surveylvl <- grpd_survey %>% 
  summarise(across(Amaranthus.sp.:Taraxacum.sp., ~ sum(.x, na.rm = T)))

# Clean up the env
rm(wild_wide_grpd, grpd_survey)


# Create distance matrices at both levels
gyr_dist_wild     <- vegdist(decostand(wild_wide_gyrlvl[2:670], "norm"), "euclidean")
survey_dist_wild  <- vegdist(decostand(wild_wide_surveylvl[2:670], "norm"), "euclidean")


# Run NMDS's

## Garden_yr level
gyr_nms_wild <- metaMDS(gyr_dist_wild, distance = "euclidean", k = 3, maxit = 999,
                        try = 5, trymax = 250, wascores = T)
plot(gyr_nms_wild, "sites")   # Produces distance 
orditorp(gyr_nms_wild, "sites")   # Gives points labels

## Survey level
survey_nms_wild <- metaMDS(survey_dist_wild, distance = "euclidean", k = 3, maxit = 999,
                           try = 5, trymax = 25, wascores = T)
plot(survey_nms_wild, "sites")   # Produces distance 
orditorp(survey_nms_wild, "sites")   # Gives points labels

# Need to remove row 338. sorry.
wild_wide_surveylvl <- wild_wide_surveylvl[-338, ]
survey_dist_wild    <- vegdist(decostand(wild_wide_surveylvl[2:670], "norm"), "euclidean")
survey_nms_wild     <- metaMDS(survey_dist_wild, distance = "euclidean", k = 3, maxit = 999,
                               try = 5, trymax = 25, wascores = T)
plot(survey_nms_wild, "sites")   # Produces distance 
orditorp(survey_nms_wild, "sites")   # Gives points labels






# Create predictor df's so we can do vector analyses
preds_wild <- wild_wide[,1:8]

gry_preds_wild <- wild_wide_gyrlvl[,1]
preds_gyr_wild <- preds_wild[!duplicated(preds_wild$garden_yr), ]
gry_preds_wild <- left_join(gry_preds_wild, preds_gyr_wild)
rm(preds_gyr_wild)
survey_preds_wild <- wild_wide_surveylvl[,1]
preds_survey_wild <- preds_wild[!duplicated(preds_wild$survey), ]
survey_preds_wild <- left_join(survey_preds_wild, preds_survey_wild)
rm(preds_survey_wild)



str(gry_preds_wild)
str(survey_preds_wild)

gry_preds_wild$status <- gry_preds_wild$impervious_1000
gry_preds_wild$status[gry_preds_wild$status > 50] <- "Urban"
gry_preds_wild$status[gry_preds_wild$status < 50] <- "Rural"
gry_preds_wild$status <- gsub("9.57", "Rural", as.character(gry_preds_wild$status))
gry_preds_wild$status2 <- paste(gry_preds_wild$city, gry_preds_wild$status)

survey_preds_wild$status <- survey_preds_wild$impervious_1000
survey_preds_wild$status[survey_preds_wild$status > 50] <- "Urban"
survey_preds_wild$status[survey_preds_wild$status < 50] <- "Rural"
survey_preds_wild$status <- gsub("9.57", "Rural", as.character(survey_preds_wild$status))
survey_preds_wild$status2 <- paste(survey_preds_wild$city, survey_preds_wild$status)







# Run the vector analyses

# With City and Urban as two categorical terms
str(gry_preds_wild)
gry_preds_wild <- select(gry_preds_wild, garden_yr, garden, city, status, everything())
str(survey_preds_wild)
survey_preds_wild <- select(survey_preds_wild, survey, garden, city, status, everything())

en_gyr_wild    <- envfit(gyr_nms_wild, gry_preds_wild[3:4], permutations = 999, na.rm = TRUE)
en_survey_wild <- envfit(survey_nms_wild, survey_preds_wild[3:4], permutations = 999, na.rm = TRUE)

en_gyr_wild
en_survey_wild

plot(gyr_nms_wild)
plot(en_gyr_wild)
plot(survey_nms_wild)
plot(en_survey_wild)












# Now make it pretty




# Garden_yr Level
data.scores.gardenyr_wild <- as.data.frame(gyr_nms_wild$points)
gyr_nms_wild$species
gyr_nms_wild$nobj
data.scores.gardenyr_wild$status <- gry_preds_wild$status2
data.scores.gardenyr_wild$city   <- gry_preds_wild$city
data.scores.gardenyr_wild$Urb    <- gry_preds_wild$status


spp.rel2_wild <- (decostand(wild_wide_gyrlvl[2:670], "norm"))
spp_scrs2_wild <- 
  sppscores(gyr_nms_wild) <- spp.rel2_wild

nms_species_wild <- as.data.frame(gyr_nms_wild$species)
nms_species_wild$Taxa <- rownames(nms_species_wild)



# Survey Level
data.scores.survey_wild <- as.data.frame(survey_nms_wild$points)
survey_nms_wild$species
survey_nms_wild$nobj
data.scores.survey_wild$status <- survey_preds_wild$status2
data.scores.survey_wild$city   <- survey_preds_wild$city
data.scores.survey_wild$Urb    <- survey_preds_wild$status


spp.rel2s_wild <- (decostand(wild_wide_surveylvl[2:670], "norm"))
spp_scrs2s_wild <- 
  sppscores(survey_nms_wild) <- spp.rel2s_wild

nms_species_s_wild <- as.data.frame(survey_nms_wild$species)
nms_species_s_wild$Taxa <- rownames(nms_species_wild)



# Plot it


# Survey
ggplot(data = data.scores.survey_wild, aes(x = MDS1, y = MDS2)) +
  # geom_segment(data = en_coord_survey, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
  #              size = 2)+
  stat_ellipse(aes(color = city), size = 1.5) +
  geom_point(data = data.scores.survey_wild, aes(size = MDS3, color = city), alpha = 0.6) +
  scale_color_manual(values = c("royalblue", "tomato3")) +
  labs(color = "City") +
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), 
        axis.text = element_text(size = 8), 
        legend.title = element_text(size = 26, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 18, colour = "grey30"))


# Gardenyr
ggplot(data = data.scores.gardenyr_wild, aes(x = MDS1, y = MDS2)) +
  stat_ellipse(aes(color = city), size = 1.5) +
  geom_point(data = data.scores.gardenyr_wild, aes(size = MDS3, color = city), alpha = 0.6) +
  scale_color_manual(values = c("royalblue", "tomato3")) +
  labs(color = "City") +
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), 
        axis.text = element_text(size = 8), 
        legend.title = element_text(size = 26, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 18, colour = "grey30"))





# Save the coordinates re-plotting

# write.csv(data.scores.gardenyr_wild, "../Data/Plants/UBH_manuscript_files/nms_datascores_gardenyr_wild.csv", row.names = F)
# write.csv(data.scores.survey_wild,   "../Data/Plants/UBH_manuscript_files/nms_datascores_survey_wild.csv", row.names = F)













### Cult ----

# Read in the wide df
cult_wide <- read.csv("../Data/Plants/UBH_manuscript_files/full_wide_cult_2020_2023.csv")
str(cult_wide)

# Drop rows with no cover
cult_wide <- filter(cult_wide, cult_total_cover > 0)

# Create a df grouped to the Garden_yr level
cult_wide_grpd <- group_by(cult_wide, garden_yr)
cult_wide_gyrlvl <- cult_wide_grpd %>% 
  summarise(across(Beta.vulgaris:Veronica.sp., ~ sum(.x, na.rm = T)))

# Do the same at the survey and garden_yr level
grpd_survey <- group_by(cult_wide, survey)
cult_wide_surveylvl <- grpd_survey %>% 
  summarise(across(Beta.vulgaris:Veronica.sp., ~ sum(.x, na.rm = T)))

# Clean up the env
rm(cult_wide_grpd, grpd_survey)


# Run dists at both levels
gyr_dist_cult     <- vegdist(decostand(cult_wide_gyrlvl[2:780], "norm"), "euclidean")
survey_dist_cult  <- vegdist(decostand(cult_wide_surveylvl[2:780], "norm"), "euclidean")


# Run NMDS's

## Garden_yr level
gyr_nms_cult <- metaMDS(gyr_dist_cult, distance = "euclidean", k = 3, maxit = 999,
                        try = 5, trymax = 250, wascores = T)
plot(gyr_nms_cult, "sites")   # Produces distance 
orditorp(gyr_nms_cult, "sites")   # Gives points labels

## Survey level
survey_nms_cult <- metaMDS(survey_dist_cult, distance = "euclidean", k = 3, maxit = 999,
                           try = 5, trymax = 25, wascores = T)
plot(survey_nms_cult, "sites")   # Produces distance 
orditorp(survey_nms_cult, "sites")   # Gives points labels



# Need to remove row 200
cult_wide_surveylvl <- cult_wide_surveylvl[-200, ]
survey_dist_cult  <- vegdist(decostand(cult_wide_surveylvl[2:780], "norm"), "euclidean")
survey_nms_cult <- metaMDS(survey_dist_cult, distance = "euclidean", k = 3, maxit = 999,
                           try = 5, trymax = 25, wascores = T)
plot(survey_nms_cult, "sites")   
orditorp(survey_nms_cult, "sites") 
# Need to remove row 181
cult_wide_surveylvl <- cult_wide_surveylvl[-181, ]
survey_dist_cult  <- vegdist(decostand(cult_wide_surveylvl[2:780], "norm"), "euclidean")
survey_nms_cult <- metaMDS(survey_dist_cult, distance = "euclidean", k = 3, maxit = 999,
                           try = 5, trymax = 25, wascores = T)
plot(survey_nms_cult, "sites")   
orditorp(survey_nms_cult, "sites") 
# Need to remove row 318, 183 and 320
cult_wide_surveylvl <- cult_wide_surveylvl[-320, ]
cult_wide_surveylvl <- cult_wide_surveylvl[-318, ]
cult_wide_surveylvl <- cult_wide_surveylvl[-183, ]
survey_dist_cult  <- vegdist(decostand(cult_wide_surveylvl[2:780], "norm"), "euclidean")
survey_nms_cult <- metaMDS(survey_dist_cult, distance = "euclidean", k = 3, maxit = 999,
                           try = 5, trymax = 25, wascores = T)
plot(survey_nms_cult, "sites")   
orditorp(survey_nms_cult, "sites") 

# In the end, 5 rows removed (1% of surveys)




# Create predictor df's so we can do vector analyses
preds_cult <- cult_wide[,1:8]

gry_preds_cult <- cult_wide_gyrlvl[,1]
preds_gyr_cult <- preds_cult[!duplicated(preds_cult$garden_yr), ]
gry_preds_cult <- left_join(gry_preds_cult, preds_gyr_cult)
rm(preds_gyr_cult)
survey_preds_cult <- cult_wide_surveylvl[,1]
preds_survey_cult <- preds_cult[!duplicated(preds_cult$survey), ]
survey_preds_cult <- left_join(survey_preds_cult, preds_survey_cult)
rm(preds_survey_cult)



str(gry_preds_cult)
str(survey_preds_cult)

gry_preds_cult$status <- gry_preds_cult$impervious_1000
gry_preds_cult$status[gry_preds_cult$status > 50] <- "Urban"
gry_preds_cult$status[gry_preds_cult$status < 50] <- "Rural"
gry_preds_cult$status <- gsub("9.57", "Rural", as.character(gry_preds_cult$status))
gry_preds_cult$status2 <- paste(gry_preds_cult$city, gry_preds_cult$status)

survey_preds_cult$status <- survey_preds_cult$impervious_1000
survey_preds_cult$status[survey_preds_cult$status > 50] <- "Urban"
survey_preds_cult$status[survey_preds_cult$status < 50] <- "Rural"
survey_preds_cult$status <- gsub("9.57", "Rural", as.character(survey_preds_cult$status))
survey_preds_cult$status2 <- paste(survey_preds_cult$city, survey_preds_cult$status)







# Run the vector analyses

# With City and Urban as two categorical terms
str(gry_preds_cult)
gry_preds_cult <- select(gry_preds_cult, garden_yr, garden, city, status, everything())
str(survey_preds_cult)
survey_preds_cult <- select(survey_preds_cult, survey, garden, city, status, everything())

en_gyr_cult    <- envfit(gyr_nms_cult, gry_preds_cult[3:4], permutations = 999, na.rm = TRUE)
en_survey_cult <- envfit(survey_nms_cult, survey_preds_cult[3:4], permutations = 999, na.rm = TRUE)

en_gyr_cult
en_survey_cult

plot(gyr_nms_cult)
plot(en_gyr_cult)
plot(survey_nms_cult)
plot(en_survey_cult)















# Now make it pretty




# Garden_yr Level
data.scores.gardenyr_cult <- as.data.frame(gyr_nms_cult$points)
gyr_nms_cult$species
gyr_nms_cult$nobj
data.scores.gardenyr_cult$status <- gry_preds_cult$status2
data.scores.gardenyr_cult$city   <- gry_preds_cult$city
data.scores.gardenyr_cult$Urb    <- gry_preds_cult$status


spp.rel2_cult <- (decostand(cult_wide_gyrlvl[2:670], "norm"))
spp_scrs2_cult <- 
  sppscores(gyr_nms_cult) <- spp.rel2_cult

nms_species_cult <- as.data.frame(gyr_nms_cult$species)
nms_species_cult$Taxa <- rownames(nms_species_cult)



# Survey Level
data.scores.survey_cult <- as.data.frame(survey_nms_cult$points)
survey_nms_cult$species
survey_nms_cult$nobj
data.scores.survey_cult$status <- survey_preds_cult$status2
data.scores.survey_cult$city   <- survey_preds_cult$city
data.scores.survey_cult$Urb    <- survey_preds_cult$status


spp.rel2s_cult <- (decostand(cult_wide_surveylvl[2:670], "norm"))
spp_scrs2s_cult <- 
  sppscores(survey_nms_cult) <- spp.rel2s_cult

nms_species_s_cult <- as.data.frame(survey_nms_cult$species)
nms_species_s_cult$Taxa <- rownames(nms_species_cult)



# Plot it


# Survey
ggplot(data = data.scores.survey_cult, aes(x = MDS1, y = MDS2)) +
  # geom_segment(data = en_coord_survey, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
  #              size = 2)+
  stat_ellipse(aes(color = city), size = 1.5) +
  geom_point(data = data.scores.survey_cult, aes(size = MDS3, color = city), alpha = 0.6) +
  scale_color_manual(values = c("royalblue", "tomato3")) +
  labs(color = "City") +
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), 
        axis.text = element_text(size = 8), 
        legend.title = element_text(size = 26, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 18, colour = "grey30"))


# Gardenyr
ggplot(data = data.scores.gardenyr_cult, aes(x = MDS1, y = MDS2)) +
  stat_ellipse(aes(color = city), size = 1.5) +
  geom_point(data = data.scores.gardenyr_cult, aes(size = MDS3, color = city), alpha = 0.6) +
  scale_color_manual(values = c("royalblue", "tomato3")) +
  labs(color = "City") +
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), 
        axis.text = element_text(size = 8), 
        legend.title = element_text(size = 26, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 18, colour = "grey30"))





# Save the coordinates re-plotting

# write.csv(data.scores.gardenyr_cult, "../Data/Plants/UBH_manuscript_files/nms_datascores_gardenyr_cult.csv", row.names = F)
# write.csv(data.scores.survey_cult,   "../Data/Plants/UBH_manuscript_files/nms_datascores_survey_cult.csv", row.names = F)















####### Plot it together ----

library(tidyverse)
library(ggpubr)

data.scores.gardenyr      <- read.csv("../Data/Plants/UBH_manuscript_files/nms_datascores_gardenyr.csv")
data.scores.survey        <- read.csv("../Data/Plants/UBH_manuscript_files/nms_datascores_survey.csv")
data.scores.gardenyr_wild <- read.csv("../Data/Plants/UBH_manuscript_files/nms_datascores_gardenyr_wild.csv")
data.scores.survey_wild   <- read.csv("../Data/Plants/UBH_manuscript_files/nms_datascores_survey_wild.csv")
data.scores.gardenyr_cult <- read.csv("../Data/Plants/UBH_manuscript_files/nms_datascores_gardenyr_cult.csv")
data.scores.survey_cult   <- read.csv("../Data/Plants/UBH_manuscript_files/nms_datascores_survey_cult.csv")




## Full
## Survey
f_s <- ggplot(data = data.scores.survey, aes(x = MDS1, y = MDS2)) +
  stat_ellipse(aes(color = city), size = 1.5) +
  geom_point(data = data.scores.survey, aes(size = MDS3, color = city), alpha = 0.6) +
  scale_color_manual(values = c("royalblue", "tomato3")) +
  labs(color = "City") +
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), 
        axis.text = element_text(size = 8), 
        legend.title = element_text(size = 30, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 30, colour = "grey30"),
        plot.title = element_text(size = 40, color = "grey20", face = "bold", hjust = 0.5))


### Gardenyr
f_g <- ggplot(data = data.scores.gardenyr, aes(x = MDS1, y = MDS2)) +
  stat_ellipse(aes(color = city), size = 1.5) +
  geom_point(data = data.scores.gardenyr, aes(size = MDS3, color = city), alpha = 0.6) +
  scale_color_manual(values = c("royalblue", "tomato3")) +
  labs(color = "City", title = "Full Community") +
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), 
        axis.text = element_text(size = 8), 
        legend.title = element_text(size = 30, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 30, colour = "grey30"),
        plot.title = element_text(size = 40, color = "grey20", face = "bold", hjust = 0.5))




## wild

### Survey
w_s <- ggplot(data = data.scores.survey_wild, aes(x = MDS1, y = MDS2)) +
  stat_ellipse(aes(color = city), size = 1.5) +
  geom_point(data = data.scores.survey_wild, aes(size = MDS3, color = city), alpha = 0.6) +
  scale_color_manual(values = c("royalblue", "tomato3")) +
  labs(color = "City") +
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), 
        axis.text = element_text(size = 8), 
        legend.title = element_text(size = 30, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 30, colour = "grey30"),
        plot.title = element_text(size = 40, color = "grey20", face = "bold", hjust = 0.5))
w_s

#### Gardenyr
w_g <- ggplot(data = data.scores.gardenyr_wild, aes(x = MDS1, y = MDS2)) +
  stat_ellipse(aes(color = city), size = 1.5) +
  geom_point(data = data.scores.gardenyr_wild, aes(size = MDS3, color = city), alpha = 0.6) +
  scale_color_manual(values = c("royalblue", "tomato3")) +
  labs(color = "City", title = "Wild Species") +
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), 
        axis.text = element_text(size = 8), 
        legend.title = element_text(size = 30, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 30, colour = "grey30"),
        plot.title = element_text(size = 40, color = "forestgreen", face = "bold", hjust = 0.5))



## Cult

### Survey
c_s <- ggplot(data = data.scores.survey_cult, aes(x = MDS1, y = MDS2)) +
  stat_ellipse(aes(color = city), size = 1.5) +
  geom_point(data = data.scores.survey_cult, aes(size = MDS3, color = city), alpha = 0.6) +
  scale_color_manual(values = c("royalblue", "tomato3")) +
  labs(color = "City") +
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), 
        axis.text = element_text(size = 8), 
        legend.title = element_text(size = 30, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 30, colour = "grey30"),
        plot.title = element_text(size = 40, color = "orange3", face = "bold", hjust = 0.5))


#### Gardenyr
c_g <- ggplot(data = data.scores.gardenyr_cult, aes(x = MDS1, y = MDS2)) +
  stat_ellipse(aes(color = city), size = 1.5) +
  geom_point(data = data.scores.gardenyr_cult, aes(size = MDS3, color = city), alpha = 0.6) +
  scale_color_manual(values = c("royalblue", "tomato3")) +
  labs(color = "City", title = "Cultivated Species") +
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), 
        axis.text = element_text(size = 8), 
        legend.title = element_text(size = 30, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 30, colour = "grey30"),
        plot.title = element_text(size = 40, color = "orange3", face = "bold", hjust = 0.5))


ggarrange(f_g, w_g, c_g,
          f_s, w_s, c_s,
          ncol = 3, nrow = 2,
          common.legend = T, legend = "right")


































# 2. Regressions ----


library(tidyverse)
library(sjPlot)
library(lme4)




# Survey level
survey_metrics <- read.csv("../Data/Plants/UBH_manuscript_files/metrics_survey_level.csv")
str(survey_metrics)


## Full
full_lm_survey <- lm(survey_beta ~ impervious_1000, data = survey_metrics)
summary(full_lm_survey)
plot(full_lm_survey)
plot(survey_metrics$survey_beta ~ survey_metrics$impervious_1000)


## Wild
wild_lm_survey <- lm(survey_beta_wild ~ impervious_1000, data = survey_metrics)
summary(wild_lm_survey)
plot(wild_lm_survey)
plot(survey_metrics$survey_beta_wild ~ survey_metrics$impervious_1000)


## Cult
cult_lm_survey <- lm(survey_beta_cult ~ impervious_1000, data = survey_metrics)
summary(cult_lm_survey)
plot(cult_lm_survey)
plot(survey_metrics$survey_beta_cult ~ survey_metrics$impervious_1000)





# Garden_yr level
garden_yr_metrics <- read.csv("../Data/Plants/UBH_manuscript_files/metrics_garden_yr_level.csv")
str(garden_yr_metrics)


## Full
full_lm_garden_yr <- lm(garden_yr_beta ~ impervious_1000, data = garden_yr_metrics)
summary(full_lm_garden_yr)
plot(full_lm_garden_yr)
plot(garden_yr_metrics$garden_yr_beta ~ garden_yr_metrics$impervious_1000)


## Wild
wild_lm_garden_yr <- lm(garden_yr_beta_wild ~ impervious_1000, data = garden_yr_metrics)
summary(wild_lm_garden_yr)
plot(wild_lm_garden_yr)
plot(garden_yr_metrics$garden_yr_beta_wild ~ garden_yr_metrics$impervious_1000)


## Cult
cult_lm_garden_yr <- lm(garden_yr_beta_cult ~ impervious_1000, data = garden_yr_metrics)
summary(cult_lm_garden_yr)
plot(cult_lm_garden_yr)
plot(garden_yr_metrics$garden_yr_beta_cult ~ garden_yr_metrics$impervious_1000)


# Compile model outputs
library(stargazer)

beta_lm_outs <- stargazer(f_lm_s, w_lm_s, c_lm_s,
                          f_lm_g, w_lm_g, c_lm_g,
                          type = "html",
                          single.row = TRUE,
                          digits = 3,
                          star.cutoffs = c(0.05, 0.01, 0.001),
                          digit.separator = "",
                          no.space = FALSE,
                          out = "../Data/Plants/UBH_manuscript_files/beta_lm_outputs.html")























## Plotting ----

library(ggpubr)


# Garden_yr level

## Full
gy_full <- ggplot(garden_yr_metrics, aes(x = impervious_1000, y = garden_yr_beta)) + 
  geom_point(aes(color = city, fill = city), size = 5, pch = 21, alpha = 0.5) + 
  scale_fill_manual(values = c("dodgerblue", "tomato3")) + 
  scale_color_manual(values = c("royalblue", "tomato4")) + 
  geom_smooth(method = "lm", color = "mediumpurple1", fill = "mediumpurple1", size = 1, alpha = 0.2) +
  theme_bw() +
  ylim(0.55,0.96)+
  theme(panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 40, face = "bold"),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 30),
        legend.title = element_blank(),
        plot.title = element_text(size = 38, face = "bold", hjust = 0.5, 
                                  color = "black")) +
  annotate("text", x = 80, y = 0.575, label = "p < 0.05", color = "mediumpurple1", size = 12, fontface = "bold") +
  labs(title = "Full Community", y = "Garden Beta Diversty", fill = "City", color = "City")
gy_full


## Wild
gy_wild <- ggplot(garden_yr_metrics, aes(x = impervious_1000, y = garden_yr_beta_wild)) + 
  geom_point(aes(color = city, fill = city), size = 5, pch = 21, alpha = 0.5) + 
  scale_fill_manual(values = c("dodgerblue", "tomato3")) + 
  scale_color_manual(values = c("royalblue", "tomato4")) + 
  geom_smooth(method = "lm", color = "mediumpurple1", fill = "mediumpurple1", size = 1, alpha = 0.2) +
  theme_bw() +
  ylim(0.55,0.96)+
  theme(panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 30),
        legend.title = element_blank(),
        plot.title = element_text(size = 38, face = "bold", hjust = 0.5, 
                                  color = "forestgreen")) +
  annotate("text", x = 82, y = 0.575, label = "p < 0.01", color = "mediumpurple1", size = 12, fontface = "bold") +
  labs(title = "Wild Species", fill = "City", color = "City")
gy_wild


## Cult
gy_cult <- ggplot(garden_yr_metrics, aes(x = impervious_1000, y = garden_yr_beta_cult)) + 
  geom_point(aes(color = city, fill = city), size = 5, pch = 21, alpha = 0.5) + 
  scale_fill_manual(values = c("dodgerblue", "tomato3")) + 
  scale_color_manual(values = c("royalblue", "tomato4")) + 
  geom_smooth(method = "lm", color = "mediumpurple1", fill = "mediumpurple1", size = 1, alpha = 0.2) +
  theme_bw() +
  ylim(0.55,0.96)+
  theme(panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 30),
        legend.title = element_blank(),
        plot.title = element_text(size = 38, face = "bold", hjust = 0.5, 
                                  color = "orange3")) +
  annotate("text", x = 85, y = 0.575, label = "p < 0.10", color = "grey60", size = 11) +
  labs(title = "Cultivated Species", fill = "City", color = "City")
gy_cult




gardenyr_arranged <- ggarrange(gy_full, gy_wild, gy_cult,
                               common.legend = T, legend = "right",  
                               ncol = 3)

annotate_figure(gardenyr_arranged, 
                bottom = text_grob("Urbanization", size = 40, face = "bold"))




s_full <- ggplot(survey_metrics, aes(x = impervious_1000, y = survey_beta)) + 
  geom_point(aes(color = city, fill = city), size = 4, pch = 21, alpha = 0.5) + 
  scale_fill_manual(values = c("dodgerblue", "tomato3")) + 
  scale_color_manual(values = c("royalblue", "tomato4")) + 
  geom_smooth(method = "lm", color = "mediumpurple1", fill = "mediumpurple1", size = 1, alpha = 0.2) +
  theme_bw() +
  ylim(0.1,0.96)+
  theme(panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 40, face = "bold"),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 30),
        legend.title = element_blank(),
        plot.title = element_blank()) +
  annotate("text", x = 90, y = 0.15, label = "ns", color = "grey60", size = 11) +
  labs(y = "Survey Beta Diversty", fill = "City", color = "City")
s_full

s_wild <- ggplot(survey_metrics, aes(x = impervious_1000, y = survey_beta_wild)) + 
  geom_point(aes(color = city, fill = city), size = 4, pch = 21, alpha = 0.5) + 
  scale_fill_manual(values = c("dodgerblue", "tomato3")) + 
  scale_color_manual(values = c("royalblue", "tomato4")) + 
  geom_smooth(method = "lm", color = "mediumpurple1", fill = "mediumpurple1", size = 1, alpha = 0.2) +
  theme_bw() +
  ylim(0.1,0.96)+
  theme(panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 30),
        legend.title = element_blank(),
        plot.title = element_blank()) +
  annotate("text", x = 85, y = 0.15, label = "p < 0.10", color = "grey60", size = 11) +
  labs(fill = "City", color = "City")

s_cult <- ggplot(survey_metrics, aes(x = impervious_1000, y = survey_beta_cult)) + 
  geom_point(aes(color = city, fill = city), size = 4, pch = 21, alpha = 0.5) + 
  scale_fill_manual(values = c("dodgerblue", "tomato3")) + 
  scale_color_manual(values = c("royalblue", "tomato4")) + 
  geom_smooth(method = "lm", color = "mediumpurple1", fill = "mediumpurple1", size = 1, alpha = 0.2) +
  theme_bw() +
  ylim(0.1,0.96)+
  theme(panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 30),
        legend.title = element_blank(),
        plot.title = element_blank()) +
  annotate("text", x = 90, y = 0.15, label = "ns", color = "grey60", size = 11) +
  labs(fill = "City", color = "City")
s_cult





both_arranged <- ggarrange(gy_full, gy_wild, gy_cult,
                            s_full, s_wild, s_cult,
                            common.legend = T, legend = "right",  
                            nrow = 2, ncol = 3)

annotate_figure(both_arranged, 
                bottom = text_grob("Urbanization", size = 40, face = "bold"))












# 3. Beta Div Contributions ---- 


## Gardenyr level

### Wild
summary(lm(garden_yr_metrics$garden_yr_beta ~ garden_yr_metrics$garden_yr_beta_wild))


gyr_betas_wild <- ggplot(garden_yr_metrics, aes(x = garden_yr_beta_wild, y = garden_yr_beta)) + 
  geom_point(aes(color = city, fill = city), size = 5, pch = 21, alpha = 0.5) + 
  scale_fill_manual(values = c("dodgerblue", "tomato3")) + 
  scale_color_manual(values = c("royalblue", "tomato4")) + 
  geom_smooth(method = "lm", color = "mediumpurple1", fill = "mediumpurple1", size = 1, alpha = 0.2) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 36, face = "bold", color = "darkgreen"),
        axis.title.y = element_text(size = 36, face = "bold", color = "black"),
        axis.text  = element_text(size = 16),
        legend.text = element_text(size = 30),
        legend.title = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = "Wild Beta Diversity",
       y = "Overall Garden \nBeta Diversity",
       color = "City", fill = "City") +
  annotate("text", x = 0.64, y = 0.98, label = "R2 = 0.85", size = 8, color = "mediumpurple1")
gyr_betas_wild

### Cult
summary(lm(garden_yr_metrics$garden_yr_beta ~ garden_yr_metrics$garden_yr_beta_cult))


gyr_betas_cult <- ggplot(garden_yr_metrics, aes(x = garden_yr_beta_cult, y = garden_yr_beta)) + 
  geom_point(aes(color = city, fill = city), size = 5, pch = 21, alpha = 0.5) + 
  scale_fill_manual(values = c("dodgerblue", "tomato3")) + 
  scale_color_manual(values = c("royalblue", "tomato4")) + 
  geom_smooth(method = "lm", color = "mediumpurple1", fill = "mediumpurple1", size = 1, alpha = 0.2) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 36, face = "bold", color = "orange3"),
        axis.title.y = element_blank(),
        axis.text  = element_text(size = 16),
        legend.text = element_text(size = 30),
        legend.title = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = "Cultivated Beta Diversity",
       y = "Overall Garden \nBeta Diversity") +
  annotate("text", x = 0.73, y = 0.97, label = "R2 = 0.16", size = 8, color = "mediumpurple1")
gyr_betas_cult



## Survey level

### Wild
summary(lm(survey_metrics$survey_beta ~ survey_metrics$survey_beta_wild))


s_betas_wild <- ggplot(survey_metrics, aes(x = survey_beta_wild, y = survey_beta)) + 
  geom_point(aes(color = city, fill = city), size = 4, pch = 21, alpha = 0.5) + 
  scale_fill_manual(values = c("dodgerblue", "tomato3")) + 
  scale_color_manual(values = c("royalblue", "tomato4")) + 
  geom_smooth(method = "lm", color = "mediumpurple1", fill = "mediumpurple1", size = 1, alpha = 0.2) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 36, face = "bold", color = "darkgreen"),
        axis.title.y = element_text(size = 36, face = "bold", color = "black"),
        axis.text  = element_text(size = 16),
        legend.text = element_text(size = 30),
        legend.title = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = "Wild Beta Diversity",
       y = "Overall Survey \nBeta Diversity",
       color = "City", fill = "City") +
  annotate("text", x = 0.35, y = 0.935, label = "R2 = 0.65", size = 8, color = "mediumpurple1")
s_betas_wild

### Cult
summary(lm(survey_metrics$survey_beta ~ survey_metrics$survey_beta_cult))


s_betas_cult <- ggplot(survey_metrics, aes(x = survey_beta_cult, y = survey_beta)) + 
  geom_point(aes(color = city, fill = city), size = 4, pch = 21, alpha = 0.5) + 
  scale_fill_manual(values = c("dodgerblue", "tomato3")) + 
  scale_color_manual(values = c("royalblue", "tomato4")) + 
  geom_smooth(method = "lm", color = "mediumpurple1", fill = "mediumpurple1", size = 1, alpha = 0.2) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 36, face = "bold", color = "orange3"),
        axis.title.y = element_blank(),
        axis.text  = element_text(size = 16),
        legend.text = element_text(size = 30),
        legend.title = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = "Cultivated Beta Diversity",
       y = "Overall Survey \nBeta Diversity",
       color = "City", fill = "City") +
  annotate("text", x = 0.18, y = 0.935, label = "R2 = 0.15", size = 8, color = "mediumpurple1")
s_betas_cult


## Combined


contributions_p <- ggarrange(gyr_betas_wild, gyr_betas_cult,
                             s_betas_wild, s_betas_cult,
                             common.legend = T, legend = "right",
                             widths = c(1, 0.89, 1, 0.89),
                             ncol = 2, nrow = 2)
contributions_p







# et voila
rm(list = ls())



















