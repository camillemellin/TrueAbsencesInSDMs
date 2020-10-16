##################################################################
# ZERO MODELS - CM 09/07/19                                      #
##################################################################

# Load libraries ------------
rm(list = ls())

library(dplyr)
library(stringr)
library(RLSPrivate)
library(psych)
library(vegan)
library(goeveg)
library(RColorBrewer)
library(tidyr)
library(boral)
library(corrplot)
library(broom)
library(visreg)
library(metafor)
library(mgcv)
library(ggplot2)
library(gridExtra)
library(pscl)
library(sp)
library(ks)
library(gstat)
library(PBSmapping)
library(sf)
library(cowplot)
library(psych)
library(spatialkernel)
library(lme4)
library(pROC)
library(caret)
library(ecospat)

scale2 <- function(x){
  (x-mean(x, na.rm = T))/sd(x, na.rm = T)
}

# Load and filter data --------------

setwd("~/Dropbox/My documents/Projects/UTAS/xxx MS Statistical report/modelling")
load('Zero models.RData')

aus <- importShapefile("~/Dropbox/My documents/Projects/UTAS/NESP/SoE/250K_coastline", readDBF=FALSE)
aus2 <- aus %>% dplyr::select(group=PID, POS=POS,long=X,lat=Y)
aus.sf <- st_read("~/Dropbox/My documents/Projects/UTAS/NESP/SoE/250K_coastline.shp")

data(fdat)
names(fdat)[5] <- "SPECIES_NAME"

load('rls-sites-plectropomus.rda')

GBR_fish_data <- read.table("GBR_fish_data.csv", header = TRUE, quote = "'", sep = ",")
GBR_site_data <- read.table("GBR_site_data.csv", header = TRUE, quote = "'", sep = ",")

GBR_fish_data$SurveyDate <- as.Date(GBR_fish_data$SurveyDate, "%d/%m/%y")

# Build metadata and list of sites surveyed both pre and post bleaching
GBR_metadata <- GBR_fish_data %>% 
  filter(Pre.or.post.bleach != "during" & include == "Y" & Method != 0) %>%
  group_by(SiteCode, Site.name, SiteLat, SiteLong, Reef, SurveyID, SurveyDate, Pre.or.post.bleach, Year) %>%
  summarize()

GBR_metadata_pre <- GBR_metadata %>% filter(Pre.or.post.bleach == "Pre")
GBR_metadata_post <- GBR_metadata %>% filter(Pre.or.post.bleach == "Post")

GBR_site.ls <- GBR_metadata %>% filter(SiteCode %in% GBR_metadata_pre$SiteCode & SiteCode %in% GBR_metadata_post$SiteCode) %>%
                group_by(SiteCode, SiteLat, SiteLong) %>% summarize()


# Compute pre vs. post species-site matrices based on different methods for zero insertion  ------

# 1- Ignore zeros
# 2- Add zeros for all species in the data matrix when not recorded at a site.
# 3- Add zeros for species recorded at each individual site on at least one occasion (this is what I have done for the population trend analysis).
# 4- Add zeros for species that occur within a convex hull that encompasses the site (i.e. extent of occurrence), where the kernel is calculated using all data.
# 5- Add zeros for species that occur within a convex hull that encompasses the site, where the kernel is calculated for that year/time slice only. This scenario takes into account changes in distributional ranges through time.

GBR_fish_data <- GBR_fish_data %>% filter(CLASS %in% fish_classes() & Method == 1)
GBR_fish_data$SPECIES_NAME <- factor(GBR_fish_data$SPECIES_NAME)

#1- Ignore zeros  
GBR_fish_site_1 <- GBR_fish_data %>%  
  filter(SiteCode %in% GBR_site.ls$SiteCode & Pre.or.post.bleach != "during" & include == "Y" & Method != 0) %>%
  group_by(SiteCode, Site.name, SiteLat, SiteLong, SurveyID, Year, Pre.or.post.bleach, SPECIES_NAME) %>%
  summarise(N = sum(N, na.rm=T)) %>%
  group_by(SiteCode, Site.name, SiteLat, SiteLong, Pre.or.post.bleach, Year, SPECIES_NAME) %>%
  summarise(N = mean(N, na.rm=T)) %>%
  group_by(SiteCode, Site.name, SiteLat, SiteLong, Pre.or.post.bleach, SPECIES_NAME) %>%
  summarise(N.1 = mean(N, na.rm=T))

#2- Add zeros everywhere when not recorded at a site
GBR_fish_site_2 <- GBR_fish_data %>%  
  filter(SiteCode %in% GBR_site.ls$SiteCode & Pre.or.post.bleach != "during" & include == "Y" & Method != 0) %>%
  group_by(SiteCode, Site.name, SiteLat, SiteLong, SurveyID, Year, Pre.or.post.bleach, SPECIES_NAME) %>%
  summarise(N = sum(N, na.rm=T)) %>%
  ungroup() %>% 
  complete(nesting(SiteCode, Site.name, SiteLat, SiteLong, SurveyID, Year, Pre.or.post.bleach), SPECIES_NAME, fill = list(N = 0)) %>%
  group_by(SiteCode, Site.name, SiteLat, SiteLong, Pre.or.post.bleach, Year, SPECIES_NAME) %>%
  summarise(N = mean(N, na.rm=T)) %>%
  group_by(SiteCode, Site.name, SiteLat, SiteLong, Pre.or.post.bleach, SPECIES_NAME) %>%
  summarise(Year = mean(Year), N.2 = mean(N, na.rm=T))

#3- Add zeros for species recorded at each individual site on at least one occasion, and at sites whithin their vicinity (i.e. within 1-degree radius)
# Or use KDE?

SiteSpecies.ls <- GBR_fish_site_1 %>% group_by(SiteCode, SiteLat, SiteLong, SPECIES_NAME) %>% summarise()

GBR_spp.ls <- names(table(SiteSpecies.ls$SPECIES_NAME)[table(SiteSpecies.ls$SPECIES_NAME)>0])

SiteSpecies.ls.rd <- SiteSpecies.ls
SiteSpecies.ls.rd$SiteLong <- round(SiteSpecies.ls.rd$SiteLong)
SiteSpecies.ls.rd$SiteLat <- round(SiteSpecies.ls.rd$SiteLat)

GBR_site.ls.rd <- GBR_site.ls
GBR_site.ls.rd$SiteLong <- round(GBR_site.ls.rd$SiteLong)
GBR_site.ls.rd$SiteLat <- round(GBR_site.ls.rd$SiteLat)

SiteSpecies.ls.rd <- SiteSpecies.ls.rd %>% left_join(GBR_site.ls.rd, by = c("SiteLong", "SiteLat"))

GBR_fish_site_3 <- subset(GBR_fish_site_2, paste(SiteCode, SPECIES_NAME, sep = "_") %in% with(SiteSpecies.ls.rd, paste(SiteCode.y, SPECIES_NAME, sep = "_")))
#GBR_fish_site_3 <- subset(GBR_fish_site_2, paste(SiteCode, SPECIES_NAME, sep = "_") %in% with(SiteSpecies.ls, paste(SiteCode, SPECIES_NAME, sep = "_")))

names(GBR_fish_site_3)[ncol(GBR_fish_site_3)] <- "N.3"

# 4- Add zeros for species that occur within a convex hull that encompasses the site (i.e. extent of occurrence), where the kernel is calculated using all data.

SiteSpecies.ls.4 <- data.frame(SPECIES_NAME = as.character(NA), SiteCode = NA)

GBR_spp_range.area <- data.frame(SPECIES_NAME = GBR_spp.ls, range.area = NA)

for (i in 1:length(GBR_spp.ls)) {
  hull.data <- subset(SiteSpecies.ls, SPECIES_NAME == GBR_spp.ls[i], select = c(SiteLong, SiteLat))
  hull <- chull(hull.data)
  hull <- c(hull, hull[1])
  GBR_spp_range.area$range.area[i] <- areapoly(as.matrix(hull.data[hull,]))$area
  GBR_site.in.hull <- point.in.polygon(GBR_site.ls$SiteLong, GBR_site.ls$SiteLat, hull.data$SiteLong[hull], hull.data$SiteLat[hull])
  SiteSpecies.ls.4 <- rbind(SiteSpecies.ls.4, data.frame(SPECIES_NAME = GBR_spp.ls[i], SiteCode = GBR_site.ls$SiteCode[GBR_site.in.hull %in% c(1,3)]))
}
SiteSpecies.ls.4 <- SiteSpecies.ls.4[-1,]

# Check convex hulls and inside/outside sites
plot(hull.data)
lines(hull.data[hull,])
polygon(hull.data[hull,], col = "lightgrey")
points(GBR_site.ls$SiteLong, GBR_site.ls$SiteLat, col = "blue", pch = 19)
points(hull.data, pch = 19, col = "green")
points(GBR_site.ls$SiteLong[GBR_site.in.hull %in% c(1,3)], GBR_site.ls$SiteLat[GBR_site.in.hull %in% c(1,3)], col = "red", pch = 19)

GBR_fish_site_4 <- subset(GBR_fish_site_2, paste(SiteCode, SPECIES_NAME, sep = "_") %in% with(SiteSpecies.ls.4, paste(SiteCode, SPECIES_NAME, sep = "_")))
names(GBR_fish_site_4)[ncol(GBR_fish_site_4)] <- "N.4"

# 5- Add zeros for species that occur within a convex hull that encompasses the site, where the kernel is calculated for that year/time slice only. This scenario takes into account changes in distributional ranges through time.

SiteSpecies.ls.pre.post <- GBR_fish_site_1 %>% group_by(SiteCode, SiteLat, SiteLong, Pre.or.post.bleach, SPECIES_NAME) %>% summarise()

SiteSpecies.ls.5 <- data.frame(SPECIES_NAME = as.character(NA), SiteCode = NA, Pre.or.post.bleach = NA)

for (i in 1:length(GBR_spp.ls)) {
  pre.hull.data <- subset(SiteSpecies.ls.pre.post, SPECIES_NAME == GBR_spp.ls[i] & Pre.or.post.bleach == "Pre", select = c(SiteLong, SiteLat))
  pre.hull <- chull(pre.hull.data)
  pre.hull <- c(pre.hull, pre.hull[1])
  GBR_site.in.pre.hull <- point.in.polygon(GBR_site.ls$SiteLong, GBR_site.ls$SiteLat, pre.hull.data$SiteLong[pre.hull], pre.hull.data$SiteLat[pre.hull])
  
  if(length(GBR_site.ls$SiteCode[GBR_site.in.pre.hull %in% c(1,3)]) > 0) {
  SiteSpecies.ls.5 <- rbind(SiteSpecies.ls.5, data.frame(SPECIES_NAME = GBR_spp.ls[i], SiteCode = GBR_site.ls$SiteCode[GBR_site.in.pre.hull %in% c(1,3)], Pre.or.post.bleach = "Pre"))
  }
  
  post.hull.data <- subset(SiteSpecies.ls.pre.post, SPECIES_NAME == GBR_spp.ls[i] & Pre.or.post.bleach == "Post", select = c(SiteLong, SiteLat))
  post.hull <- chull(post.hull.data)
  post.hull <- c(post.hull, post.hull[1])
  GBR_site.in.post.hull <- point.in.polygon(GBR_site.ls$SiteLong, GBR_site.ls$SiteLat, post.hull.data$SiteLong[post.hull], post.hull.data$SiteLat[post.hull])
  
  if(length(GBR_site.ls$SiteCode[GBR_site.in.post.hull %in% c(1,3)]) > 0) {
    SiteSpecies.ls.5 <- rbind(SiteSpecies.ls.5, data.frame(SPECIES_NAME = GBR_spp.ls[i], SiteCode = GBR_site.ls$SiteCode[GBR_site.in.post.hull %in% c(1,3)], Pre.or.post.bleach = "Post"))
  }
  }

SiteSpecies.ls.5 <- SiteSpecies.ls.5[-1,]


# Check convex hulls and inside/outside sites
 par(mfcol = c(2,1))
 plot(pre.hull.data)
 lines(pre.hull.data[pre.hull,])
 polygon(pre.hull.data[pre.hull,], col = "lightgrey")
 points(GBR_site.ls$SiteLong, GBR_site.ls$SiteLat, col = "blue", pch = 19)
 points(pre.hull.data, pch = 19, col = "green")
 points(GBR_site.ls$SiteLong[GBR_site.in.pre.hull %in% c(1,3)], GBR_site.ls$SiteLat[GBR_site.in.pre.hull %in% c(1,3)], col = "red", pch = 19)
 
 plot(post.hull.data)
 lines(post.hull.data[post.hull,])
 polygon(post.hull.data[post.hull,], col = "lightgrey")
 points(GBR_site.ls$SiteLong, GBR_site.ls$SiteLat, col = "blue", pch = 19)
 points(post.hull.data, pch = 19, col = "green")
 points(GBR_site.ls$SiteLong[GBR_site.in.post.hull %in% c(1,3)], GBR_site.ls$SiteLat[GBR_site.in.post.hull %in% c(1,3)], col = "red", pch = 19)

GBR_fish_site_5 <- subset(GBR_fish_site_2, paste(SiteCode, Pre.or.post.bleach, SPECIES_NAME, sep = "_") %in% with(SiteSpecies.ls.5, paste(SiteCode, Pre.or.post.bleach, SPECIES_NAME, sep = "_")))
names(GBR_fish_site_5)[ncol(GBR_fish_site_5)] <- "N.5"


# Build single table with 5 abundance estimates, one for each method

GBR_fish_site_all <- GBR_fish_site_2 %>% 
  left_join(GBR_fish_site_1, by = c("SiteCode","Site.name","SiteLat","SiteLong","SPECIES_NAME","Pre.or.post.bleach")) %>%
  left_join(GBR_fish_site_3, by = c("SiteCode","Site.name","SiteLat","SiteLong","SPECIES_NAME","Pre.or.post.bleach")) %>%
  left_join(GBR_fish_site_4, by = c("SiteCode","Site.name","SiteLat","SiteLong","SPECIES_NAME","Pre.or.post.bleach")) %>%
  left_join(GBR_fish_site_5, by = c("SiteCode","Site.name","SiteLat","SiteLong","SPECIES_NAME","Pre.or.post.bleach")) %>%
  dplyr::select("SiteCode","Site.name","SiteLat","SiteLong","Year", "Pre.or.post.bleach","SPECIES_NAME","N.1","N.2","N.3","N.4","N.5")
                    
GBR_fish_site_P <- GBR_fish_site_all  
GBR_fish_site_P[,7:11][GBR_fish_site_P[,7:11] > 0] <- 1
names(GBR_fish_site_P)[7:11] <- c("P.1", "P.2", "P.3", "P.4", "P.5")

# Illustrate the method with Ctenochaetus cyanocheilus -----

plot.data.pre <- GBR_fish_site_all %>% filter(SPECIES_NAME == "Ctenochaetus cyanocheilus" & Pre.or.post.bleach == "Pre")
plot.data.post <- GBR_fish_site_all %>% filter(SPECIES_NAME == "Ctenochaetus cyanocheilus" & Pre.or.post.bleach == "Post")

plot.data.pre[,7:11][plot.data.pre[,7:11] > 0] <- 1
plot.data.post[,7:11][plot.data.post[,7:11] > 0] <- 1

map.1_pre <- ggplot() +
  geom_polygon(data=aus2, aes(long, lat, group=group), fill="lightgray", color="darkgray") +
  coord_map(xlim=c(143,156), ylim=c(-22,-10)) +
  # xlab(expression(paste(Longitude^o, ~'E'))) +
  # ylab(expression(paste(Latitude^o, ~'S'))) +
  geom_point(data=GBR_site.ls, aes(SiteLong, SiteLat), size=1, shape=19, colour="dimgrey") +
  geom_point(data = plot.data.pre[plot.data.pre$N.1 == 0,], aes(SiteLong, SiteLat), size = 2, shape = 19, colour="cornflowerblue")+
  geom_point(data = plot.data.pre[plot.data.pre$N.1 == 1,], aes(SiteLong, SiteLat), size = 3, shape = 19, colour="tomato")+
  theme(text=element_text(size=12,  family="Calibri"), plot.margin = unit(c(.1,.1,.1,.1), "cm"), axis.title = element_blank()) + theme_void()
  

map.1_post <- ggplot() +
  geom_polygon(data=aus2, aes(long, lat, group=group), fill="lightgray", color="darkgray") +
  coord_map(xlim=c(143,156), ylim=c(-22,-10)) +
  # xlab(expression(paste(Longitude^o, ~'E'))) +
  # ylab(expression(paste(Latitude^o, ~'S'))) +
  geom_point(data=GBR_site.ls, aes(SiteLong, SiteLat), size=1, shape=19, colour="dimgrey") +
  geom_point(data = plot.data.post[plot.data.post$N.1 == 0,], aes(SiteLong, SiteLat), size = 2, shape = 19, colour="cornflowerblue")+
  geom_point(data = plot.data.post[plot.data.post$N.1 == 1,], aes(SiteLong, SiteLat), size = 3, shape = 19, colour="tomato")+
  theme(text=element_text(size=12,  family="Calibri"), plot.margin = unit(c(.1,.1,.1,.1), "cm"), axis.title = element_blank())+ theme_void()

map.2_pre <- ggplot() +
  geom_polygon(data=aus2, aes(long, lat, group=group), fill="lightgray", color="darkgray") +
  coord_map(xlim=c(143,156), ylim=c(-22,-10)) +
  xlab(expression(paste(Longitude^o, ~'E'))) +
  ylab(expression(paste(Latitude^o, ~'S'))) +
  geom_point(data=GBR_site.ls, aes(SiteLong, SiteLat), size=1, shape=19, colour="dimgrey") +
  geom_point(data = plot.data.pre[plot.data.pre$N.2 == 0,], aes(SiteLong, SiteLat), size = 2, shape = 19, colour="cornflowerblue")+
  geom_point(data = plot.data.pre[plot.data.pre$N.2 == 1,], aes(SiteLong, SiteLat), size = 3, shape = 19, colour="tomato")+
  theme(text=element_text(size=12,  family="Calibri"), plot.margin = unit(c(.1,.1,.1,.1), "cm"), axis.title = element_blank())+ theme_void()

map.2_post <- ggplot() +
  geom_polygon(data=aus2, aes(long, lat, group=group), fill="lightgray", color="darkgray") +
  coord_map(xlim=c(143,156), ylim=c(-22,-10)) +
  xlab(expression(paste(Longitude^o, ~'E'))) +
  ylab(expression(paste(Latitude^o, ~'S'))) +
  geom_point(data=GBR_site.ls, aes(SiteLong, SiteLat), size=1, shape=19, colour="dimgrey") +
  geom_point(data = plot.data.post[plot.data.post$N.2 == 0,], aes(SiteLong, SiteLat), size = 2, shape = 19, colour="cornflowerblue")+
  geom_point(data = plot.data.post[plot.data.post$N.2 == 1,], aes(SiteLong, SiteLat), size = 3, shape = 19, colour="tomato")+
  theme(text=element_text(size=12,  family="Calibri"), plot.margin = unit(c(.1,.1,.1,.1), "cm"), axis.title = element_blank())+ theme_void()

map.3_pre <- ggplot() +
  geom_polygon(data=aus2, aes(long, lat, group=group), fill="lightgray", color="darkgray") +
  coord_map(xlim=c(143,156), ylim=c(-22,-10)) +
  xlab(expression(paste(Longitude^o, ~'E'))) +
  ylab(expression(paste(Latitude^o, ~'S'))) +
  geom_point(data=GBR_site.ls, aes(SiteLong, SiteLat), size=1, shape=19, colour="dimgrey") +
  geom_point(data = plot.data.pre[plot.data.pre$N.3 == 0,], aes(SiteLong, SiteLat), size = 2, shape = 19, colour="cornflowerblue")+
  geom_point(data = plot.data.pre[plot.data.pre$N.3 == 1,], aes(SiteLong, SiteLat), size = 3, shape = 19, colour="tomato")+
  theme(text=element_text(size=12,  family="Calibri"), plot.margin = unit(c(.1,.1,.1,.1), "cm"), axis.title = element_blank())+ theme_void()

map.3_post <- ggplot() +
  geom_polygon(data=aus2, aes(long, lat, group=group), fill="lightgray", color="darkgray") +
  coord_map(xlim=c(143,156), ylim=c(-22,-10)) +
  xlab(expression(paste(Longitude^o, ~'E'))) +
  ylab(expression(paste(Latitude^o, ~'S'))) +
  geom_point(data=GBR_site.ls, aes(SiteLong, SiteLat), size=1, shape=19, colour="dimgrey") +
  geom_point(data = plot.data.post[plot.data.post$N.3 == 0,], aes(SiteLong, SiteLat), size = 2, shape = 19, colour="cornflowerblue")+
  geom_point(data = plot.data.post[plot.data.post$N.3 == 1,], aes(SiteLong, SiteLat), size = 3, shape = 19, colour="tomato")+
  theme(text=element_text(size=12,  family="Calibri"), plot.margin = unit(c(.1,.1,.1,.1), "cm"), axis.title = element_blank())+ theme_void()

map.4_pre <- ggplot() +
  geom_polygon(data=aus2, aes(long, lat, group=group), fill="lightgray", color="darkgray") +
  coord_map(xlim=c(143,156), ylim=c(-22,-10)) +
  xlab(expression(paste(Longitude^o, ~'E'))) +
  ylab(expression(paste(Latitude^o, ~'S'))) +
  geom_polygon(data = hull.data[hull,], aes(SiteLong, SiteLat), fill = "lightblue", alpha = .8)+
  geom_point(data=GBR_site.ls, aes(SiteLong, SiteLat), size=1, shape=19, colour="dimgrey") +
  geom_point(data = plot.data.pre[plot.data.pre$N.4 == 0,], aes(SiteLong, SiteLat), size = 2, shape = 19, colour="cornflowerblue")+
  geom_point(data = plot.data.pre[plot.data.pre$N.4 == 1,], aes(SiteLong, SiteLat), size = 3, shape = 19, colour="tomato")+
  theme(text=element_text(size=12,  family="Calibri"), plot.margin = unit(c(.1,.1,.1,.1), "cm"), axis.title = element_blank())+ theme_void()

map.4_post <- ggplot() +
  geom_polygon(data=aus2, aes(long, lat, group=group), fill="lightgray", color="darkgray") +
  coord_map(xlim=c(143,156), ylim=c(-22,-10)) +
  xlab(expression(paste(Longitude^o, ~'E'))) +
  ylab(expression(paste(Latitude^o, ~'S'))) +
  geom_polygon(data = hull.data[hull,], aes(SiteLong, SiteLat), fill = "lightblue", alpha = .8)+
  geom_point(data=GBR_site.ls, aes(SiteLong, SiteLat), size=1, shape=19, colour="dimgrey") +
  geom_point(data = plot.data.post[plot.data.post$N.4 == 0,], aes(SiteLong, SiteLat), size = 2, shape = 19, colour="cornflowerblue")+
  geom_point(data = plot.data.post[plot.data.post$N.4 == 1,], aes(SiteLong, SiteLat), size = 3, shape = 19, colour="tomato")+
  theme(text=element_text(size=12,  family="Calibri"), plot.margin = unit(c(.1,.1,.1,.1), "cm"), axis.title = element_blank())+ theme_void()

map.5_pre <- ggplot() +
  geom_polygon(data=aus2, aes(long, lat, group=group), fill="lightgray", color="darkgray") +
  coord_map(xlim=c(143,156), ylim=c(-22,-10)) +
  xlab(expression(paste(Longitude^o, ~'E'))) +
  ylab(expression(paste(Latitude^o, ~'S'))) +
  geom_polygon(data = pre.hull.data[pre.hull,], aes(SiteLong, SiteLat), fill = "lightblue", alpha = .8)+
  geom_point(data=GBR_site.ls, aes(SiteLong, SiteLat), size=1, shape=19, colour="dimgrey") +
  geom_point(data = plot.data.pre[plot.data.pre$N.5 == 0,], aes(SiteLong, SiteLat), size = 2, shape = 19, colour="cornflowerblue")+
  geom_point(data = plot.data.pre[plot.data.pre$N.5 == 1,], aes(SiteLong, SiteLat), size = 3, shape = 19, colour="tomato")+
  theme(text=element_text(size=12,  family="Calibri"), plot.margin = unit(c(.1,.1,.1,.1), "cm"), axis.title = element_blank())+ theme_void()

map.5_post <- ggplot() +
  geom_polygon(data=aus2, aes(long, lat, group=group), fill="lightgray", color="darkgray") +
  coord_map(xlim=c(143,156), ylim=c(-22,-10)) +
  xlab(expression(paste(Longitude^o, ~'E'))) +
  ylab(expression(paste(Latitude^o, ~'S'))) +
  geom_polygon(data = post.hull.data[post.hull,], aes(SiteLong, SiteLat), fill = "lightblue", alpha = .8)+
  geom_point(data=GBR_site.ls, aes(SiteLong, SiteLat), size=1, shape=19, colour="dimgrey") +
  geom_point(data = plot.data.post[plot.data.post$N.5 == 0,], aes(SiteLong, SiteLat), size = 2, shape = 19, colour="cornflowerblue")+
  geom_point(data = plot.data.post[plot.data.post$N.5 == 1,], aes(SiteLong, SiteLat), size = 3, shape = 19, colour="tomato")+
  theme(text=element_text(size=12,  family="Calibri"), plot.margin = unit(c(.1,.1,.1,.1), "cm"), axis.title = element_blank())+ theme_void()

map.pre.post <- plot_grid(plotlist = list(map.1_pre, map.1_post, map.2_pre, map.2_post, map.3_pre, map.3_post, map.4_pre, map.4_post, map.5_pre,
                                          map.5_post), ncol=2, nrow=5)

# Compare species frequency distribution --------

spp_freq_pre <- GBR_fish_site_P %>% filter(Pre.or.post.bleach == "Pre") %>% group_by(SPECIES_NAME) %>% 
            summarize(F.2 = sum(P.2, na.rm = T)/n_distinct(SiteCode[P.2 %in% c(0,1)]),
                      F.3 = sum(P.3, na.rm = T)/n_distinct(SiteCode[P.3 %in% c(0,1)]),
                      F.4 = sum(P.4, na.rm = T)/n_distinct(SiteCode[P.4 %in% c(0,1)]),
                      F.5 = sum(P.5, na.rm = T)/n_distinct(SiteCode[P.5 %in% c(0,1)])) %>%
                      filter(F.5 < 1)

d.2.pre <- with(spp_freq_pre, density(F.2, na.rm = T))
d.3.pre <- with(spp_freq_pre, density(F.3, na.rm = T))
d.4.pre <- with(spp_freq_pre, density(F.4, na.rm = T))
d.5.pre <- with(spp_freq_pre, density(F.5, na.rm = T))

spp_freq_post <- GBR_fish_site_P %>% filter(Pre.or.post.bleach == "Post") %>% group_by(SPECIES_NAME) %>% 
  summarize(F.2 = sum(P.2, na.rm = T)/n_distinct(SiteCode[P.2 %in% c(0,1)]),
            F.3 = sum(P.3, na.rm = T)/n_distinct(SiteCode[P.3 %in% c(0,1)]),
            F.4 = sum(P.4, na.rm = T)/n_distinct(SiteCode[P.4 %in% c(0,1)]),
            F.5 = sum(P.5, na.rm = T)/n_distinct(SiteCode[P.5 %in% c(0,1)])) %>%
            filter(F.5 < 1)

d.2.post <- with(spp_freq_post, density(F.2, na.rm = T))
d.3.post <- with(spp_freq_post, density(F.3, na.rm = T))
d.4.post <- with(spp_freq_post, density(F.4, na.rm = T))
d.5.post <- with(spp_freq_post, density(F.5, na.rm = T))

plot(d.2.pre, type = "l", ylim = c(0,5))
lines(d.2.post, lty = 2)
lines(d.3.pre, col = "green")
lines(d.3.post, col = "green", lty = 2)
lines(d.4.pre, col = "red")
lines(d.4.post, col = "red", lty = 2)
lines(d.5.pre, col = "orange")
lines(d.5.post, col = "orange", lty = 2)

pairs.panels(spp_freq_pre[,-1])
pairs.panels(spp_freq_post[,-1])

spp_freq_pre <- spp_freq_pre %>% left_join(GBR_spp_range.area)
spp_freq_post <- spp_freq_post %>% left_join(GBR_spp_range.area)

density.plot <- ggplot() + 
  geom_density(data = spp_freq_pre, aes(F.2, color = "M.2", lty = "Before"), size = 1)+
  geom_density(data = spp_freq_post, aes(F.2, color = "M.2", lty = "After"), size = 1)+
  geom_density(data = spp_freq_pre, aes(F.3, color = "M.3", lty = "Before"), size = 1)+
  geom_density(data = spp_freq_post, aes(F.3, color = "M.3", lty = "After"), size = 1)+
  geom_density(data = spp_freq_pre, aes(F.4, color = "M.4", lty = "Before"), size = 1)+
  geom_density(data = spp_freq_post, aes(F.4, color = "M.4", lty = "After"), size = 1)+
  geom_density(data = spp_freq_pre, aes(F.5, color = "M.5", lty = "Before"), size = 1)+
  geom_density(data = spp_freq_post, aes(F.5, color = "M.5", lty = "After"), size = 1)+
  xlab("Species frequency")+
  theme(legend.position = 'right') +
  scale_color_manual(values = c('#172A3A','#006166','#508991',"#09BC8A"), 
                     labels = c("M.2", "M.3", "M.4", "M.5"),
                     name = "Method")+
  scale_linetype_manual(values = factor(c(1,3)), labels = c("Before", "After"), name = "Bleaching")

  
biplot <- ggplot() +
  geom_point(data = spp_freq_pre, aes(x = F.2, y = F.4, color = range.area))+
  xlab("Species frequency (M1)") + ylab("Species frequency (M4)")

fig2 <- plot_grid(density.plot, biplot, ncol=1, nrow=2)
fig2


# Build corresponding Site matrix ---------------

# Process SST data from coral trout dataset
dat_plectro2 <- dat_plectro  %>% dplyr::select(SiteCode, PrePost, sst_mean, sst_year, sst_anom, Live_hard_coral) %>%
  group_by(SiteCode, PrePost) %>% summarize_all(mean) %>% data.frame()

dat_plectro2$PrePost <- factor(dat_plectro2$PrePost, levels = c("Before", "After", "Pre", "Post"))
dat_plectro2$PrePost <- ifelse(dat_plectro2$PrePost == "Before", "Pre", "Post")
dat_plectro2$PrePost <- factor(dat_plectro2$PrePost)

all.dat <- GBR_fish_site_P %>% left_join(dat_plectro2, by = c('SiteCode', 'Pre.or.post.bleach' = 'PrePost'))
  
# Logistic GLMM for each species ---------

# See link on Cohen's kappa: https://stats.stackexchange.com/questions/82162/cohens-kappa-in-plain-english

# Calculate spp.frequency (only include spp occurring at 15 sites or more, i.e. 155 species)
spp.frequency <- GBR_fish_site_P %>% group_by(SPECIES_NAME) %>% 
  summarize(Nocc = sum(P.1, na.rm = T), 
            Nocc.pre = sum(P.1[Pre.or.post.bleach == "Pre"], na.rm = T), 
            Nocc.post = sum(P.1[Pre.or.post.bleach == "Post"], na.rm = T)) %>%
  filter(Nocc > 20 & Nocc.pre > 10 & Nocc.post > 10) %>%
  left_join(GBR_spp_range.area)

all.dat.sub <- all.dat %>% filter(SPECIES_NAME %in% spp.frequency$SPECIES_NAME)
all.dat.sub <- data.frame(all.dat.sub)
for (i in 12:15) all.dat.sub[,i] <- scale2(all.dat.sub[,i])

M2.coef <- M3.coef <- M4.coef <- M5.coef <- 
M2.pval <- M3.pval <- M4.pval <- M5.pval <- 
  data.frame(SPECIES_NAME = spp.frequency$SPECIES_NAME, "(Intercept)"=NA, "Live_hard_coral"=NA, "sst_mean"=NA, "sst_anom"=NA, "sst_mean:sst_anom"=NA)

AUC <- accuracy <- kappa <- boyce <- tss <- data.frame(SPECIES_NAME = spp.frequency$SPECIES_NAME, "m2"=NA, "m3"=NA, "m4"=NA, "m5"=NA)

for (i in 1:length(spp.frequency$SPECIES_NAME)) {
    
    sp.dat <- subset(all.dat.sub, SPECIES_NAME == spp.frequency$SPECIES_NAME[i])
    
    m2 <- glmer(P.2 ~ Live_hard_coral + sst_mean*sst_anom + (1|Pre.or.post.bleach), data = sp.dat, family = binomial)
    m3 <- glmer(P.3 ~ Live_hard_coral + sst_mean*sst_anom + (1|Pre.or.post.bleach), data = sp.dat, family = binomial)
    m4 <- glmer(P.4 ~ Live_hard_coral + sst_mean*sst_anom + (1|Pre.or.post.bleach), data = sp.dat, family = binomial)
    m5 <- glmer(P.5 ~ Live_hard_coral + sst_mean*sst_anom + (1|Pre.or.post.bleach), data = sp.dat, family = binomial)
    
    M2.coef[i,-1] <- summary(m2)$coefficients[,1]
    M3.coef[i,-1] <- summary(m3)$coefficients[,1]
    M4.coef[i,-1] <- summary(m4)$coefficients[,1]
    M5.coef[i,-1] <- summary(m5)$coefficients[,1]
    
    M2.pval[i,-1] <- summary(m2)$coefficients[,4]
    M3.pval[i,-1] <- summary(m3)$coefficients[,4]
    M4.pval[i,-1] <- summary(m4)$coefficients[,4]
    M5.pval[i,-1] <- summary(m5)$coefficients[,4]
    
    pred.m2 <- predict(m2, newdata = sp.dat, type = "response")
    pred.m3 <- predict(m3, newdata = sp.dat, type = "response")
    pred.m4 <- predict(m4, newdata = sp.dat, type = "response")
    pred.m5 <- predict(m5, newdata = sp.dat, type = "response")
    
    # Remove NA values from predictions and observations for calculating accuracy metrics
    obs.m2 <- sp.dat$P.2[!is.na(pred.m2)]
    obs.m3 <- sp.dat$P.3[!is.na(pred.m3)]
    obs.m4 <- sp.dat$P.4[!is.na(pred.m4)]
    obs.m5 <- sp.dat$P.5[!is.na(pred.m5)]
    
    pred.m2 <- na.omit(pred.m2)
    pred.m3 <- na.omit(pred.m3)
    pred.m4 <- na.omit(pred.m4)
    pred.m5 <- na.omit(pred.m5)
    
    AUC$m2[i] <- auc(roc(obs.m2, pred.m2, quiet = T))
    AUC$m3[i] <- auc(roc(obs.m3, pred.m3, quiet = T))
    AUC$m4[i] <- auc(roc(obs.m4, pred.m4, quiet = T))
    AUC$m5[i] <- auc(roc(obs.m5, pred.m5, quiet = T))
    
    boyce$m2[i] <- ecospat.boyce(as.numeric(pred.m2), as.numeric(pred.m2[which(obs.m2 == 1)]), 
                                 nclass=0, window.w="default", res=100, PEplot = F)$Spearman.cor
    
    boyce$m3[i] <- ecospat.boyce(as.numeric(pred.m3), as.numeric(pred.m3[which(obs.m3 == 1)]), 
                                 nclass=0, window.w="default", res=100, PEplot = F)$Spearman.cor
    
    boyce$m4[i] <- ecospat.boyce(as.numeric(pred.m4), as.numeric(pred.m4[which(obs.m4 == 1)]), 
                                 nclass=0, window.w="default", res=100, PEplot = F)$Spearman.cor
    
    boyce$m5[i] <- ecospat.boyce(as.numeric(pred.m5), as.numeric(pred.m5[which(obs.m5 == 1)]), 
                                 nclass=0, window.w="default", res=100, PEplot = F)$Spearman.cor
    
    p2 <- as.numeric(pred.m2>0.5)
    accuracy$m2[i] <- mean(p2==obs.m2, na.rm = T)
    c2 <- confusionMatrix(factor(p2), factor(obs.m2))
    kappa$m2[i] <- c2$overall[2]
    tss$m2[i] <- c2$byClass["Sensitivity"] + c2$byClass["Specificity"] - 1
    
    p3 <- as.numeric(pred.m3>0.5)
    accuracy$m3[i] <- mean(p3==obs.m3, na.rm = T)
    c3 <- confusionMatrix(factor(p3), factor(obs.m3))
    kappa$m3[i] <- c3$overall[3]
    tss$m3[i] <- c3$byClass["Sensitivity"] + c3$byClass["Specificity"] - 1
    
    p4 <- as.numeric(pred.m4>0.5)
    accuracy$m4[i] <- mean(p4==obs.m4, na.rm = T)
    c4 <- confusionMatrix(factor(p4), factor(obs.m4))
    kappa$m4[i] <- c4$overall[4]
    tss$m4[i] <- c4$byClass["Sensitivity"] + c4$byClass["Specificity"] - 1
    
    p5 <- as.numeric(pred.m5>0.5)
    accuracy$m5[i] <- mean(p5==obs.m5, na.rm = T)
    c5 <- confusionMatrix(factor(p5), factor(obs.m5))
    kappa$m5[i] <- c5$overall[5]
    tss$m5[i] <- c5$byClass["Sensitivity"] + c5$byClass["Specificity"] - 1
    
    print(i)
}

M2.coef[,-1][M2.pval[,-1] > 0.05] <- NA
M3.coef[,-1][M3.pval[,-1] > 0.05] <- NA
M4.coef[,-1][M4.pval[,-1] > 0.05] <- NA
M5.coef[,-1][M5.pval[,-1] > 0.05] <- NA

# Logistic GLMM for each species: fit on PRE data, predict POST data ---------


kappa.prepost <- tss.prepost <- data.frame(SPECIES_NAME = spp.frequency$SPECIES_NAME, "m2"=NA, "m3"=NA, "m4"=NA, "m5"=NA)

for (i in 1:length(spp.frequency$SPECIES_NAME)) {
  
  sp.dat <- subset(all.dat.sub, SPECIES_NAME == spp.frequency$SPECIES_NAME[i])
  
  m2 <- glm(P.2 ~ Live_hard_coral + sst_mean*sst_anom, data = sp.dat[sp.dat$Pre.or.post.bleach == "Pre",], family = binomial)
  m3 <- glm(P.3 ~ Live_hard_coral + sst_mean*sst_anom, data = sp.dat[sp.dat$Pre.or.post.bleach == "Pre",], family = binomial)
  m4 <- glm(P.4 ~ Live_hard_coral + sst_mean*sst_anom, data = sp.dat[sp.dat$Pre.or.post.bleach == "Pre",], family = binomial)
  m5 <- glm(P.5 ~ Live_hard_coral + sst_mean*sst_anom, data = sp.dat[sp.dat$Pre.or.post.bleach == "Pre",], family = binomial)
  
  pred.m2 <- predict(m2, newdata = sp.dat[sp.dat$Pre.or.post.bleach == "Post",], type = "response")
  pred.m3 <- predict(m3, newdata = sp.dat[sp.dat$Pre.or.post.bleach == "Post",], type = "response")
  pred.m4 <- predict(m4, newdata = sp.dat[sp.dat$Pre.or.post.bleach == "Post",], type = "response")
  pred.m5 <- predict(m5, newdata = sp.dat[sp.dat$Pre.or.post.bleach == "Post",], type = "response")
  
  # Remove NA values from predictions and observations for calculating accuracy metrics
  obs.m2 <- sp.dat$P.2[sp.dat$Pre.or.post.bleach == "Pre"][!is.na(pred.m2)]
  obs.m3 <- sp.dat$P.3[sp.dat$Pre.or.post.bleach == "Pre"][!is.na(pred.m3)]
  obs.m4 <- sp.dat$P.4[sp.dat$Pre.or.post.bleach == "Pre"][!is.na(pred.m4)]
  obs.m5 <- sp.dat$P.5[sp.dat$Pre.or.post.bleach == "Pre"][!is.na(pred.m5)]
  
  pred.m2 <- na.omit(pred.m2)
  pred.m3 <- na.omit(pred.m3)
  pred.m4 <- na.omit(pred.m4)
  pred.m5 <- na.omit(pred.m5)
  
  p2 <- as.numeric(pred.m2>0.5)
  c2 <- confusionMatrix(factor(p2), factor(obs.m2))
  kappa.prepost$m2[i] <- c2$overall[2]
  tss.prepost$m2[i] <- c2$byClass["Sensitivity"] + c2$byClass["Specificity"] - 1
  
  p3 <- as.numeric(pred.m3>0.5)
  c3 <- confusionMatrix(factor(p3), factor(obs.m3))
  kappa.prepost$m3[i] <- c3$overall[3]
  tss.prepost$m3[i] <- c3$byClass["Sensitivity"] + c3$byClass["Specificity"] - 1
  
  p4 <- as.numeric(pred.m4>0.5)
  c4 <- confusionMatrix(factor(p4), factor(obs.m4))
  kappa.prepost$m4[i] <- c4$overall[4]
  tss.prepost$m4[i] <- c4$byClass["Sensitivity"] + c4$byClass["Specificity"] - 1
  
  p5 <- as.numeric(pred.m5>0.5)
  c5 <- confusionMatrix(factor(p5), factor(obs.m5))
  kappa.prepost$m5[i] <- c5$overall[5]
  tss.prepost$m5[i] <- c5$byClass["Sensitivity"] + c5$byClass["Specificity"] - 1
  
  print(i)
}

kappa.prepost[kappa.prepost < 0] <- 0
tss.prepost[tss.prepost < 0] <- 0

# Heat matrices of model coefficients (Fig. 3) ---------------
# Heat matrix for %Coral
coral.all <- data.frame(SPECIES_NAME = M2.coef$SPECIES_NAME,
                        M2 = M2.coef$Live_hard_coral,
                        M3 = M3.coef$Live_hard_coral,
                        M4 = M4.coef$Live_hard_coral,
                        M5 = M5.coef$Live_hard_coral)

coral.all <- subset(coral.all, !(is.na(M3) & is.na(M4) & is.na(M5)))
coral.all[is.na(coral.all)] <- 0

coral.all <- coral.all[order(coral.all$M5, coral.all$M4, coral.all$M3, decreasing = F),]

coral.all <- data.matrix(coral.all[,-1])
coral.all[coral.all > quantile(coral.all, .9)] <- quantile(coral.all, .9)
coral.all[coral.all < quantile(coral.all, .1)] <- quantile(coral.all, .1)

corrplot(t(coral.all), col = rev(brewer.pal(11,"RdBu")),
         method = "color", is.corr= FALSE,
         cl.pos = "n", addgrid.col = "lightgrey", 
         tl.pos = "n")

cor(coral.all)

coral.all.recl <- coral.all
coral.all.recl[coral.all.recl < 0] <- 1
coral.all.recl[coral.all.recl > 0] <- 1

# False positives
# M3 = 10%
length(coral.all.recl[,3][coral.all.recl[,3] == 1 & coral.all.recl[,4] == 0]) * 100/dim(coral.all.recl)[1]
# M2 = 15.8%
length(coral.all.recl[,2][coral.all.recl[,2] == 1 & coral.all.recl[,4] == 0]) * 100/dim(coral.all.recl)[1]
# M1 = 13.3%
length(coral.all.recl[,1][coral.all.recl[,1] == 1 & coral.all.recl[,4] == 0]) * 100/dim(coral.all.recl)[1]

# False negatives
# M3 = 10%
length(coral.all.recl[,3][coral.all.recl[,3] == 0 & coral.all.recl[,4] == 1]) * 100/dim(coral.all.recl)[1]
# M2 = 15.8%
length(coral.all.recl[,2][coral.all.recl[,2] == 0 & coral.all.recl[,4] == 1]) * 100/dim(coral.all.recl)[1]
# M1 = 13.3%
length(coral.all.recl[,1][coral.all.recl[,1] == 0 & coral.all.recl[,4] == 1]) * 100/dim(coral.all.recl)[1]



# Heat matrix for %sst
sst.all <- data.frame(SPECIES_NAME = M2.coef$SPECIES_NAME,
                        M2 = M2.coef$sst_anom,
                        M3 = M3.coef$sst_anom,
                        M4 = M4.coef$sst_anom,
                        M5 = M5.coef$sst_anom)

sst.all <- subset(sst.all, !(is.na(M3) & is.na(M4) & is.na(M5)))
sst.all[is.na(sst.all)] <- 0

sst.all <- sst.all[order(sst.all$M5, sst.all$M4, sst.all$M3, decreasing = F),]

sst.all <- data.matrix(sst.all[,-1])
sst.all[sst.all > quantile(sst.all, .9)] <- quantile(sst.all, .9)
sst.all[sst.all < quantile(sst.all, .1)] <- quantile(sst.all, .1)

corrplot(t(sst.all), col = rev(brewer.pal(11,"RdBu")),
         method = "color", is.corr= FALSE,
         cl.pos = "n", addgrid.col = "lightgrey", 
         tl.pos = "n")

cor(sst.all)

sst.all.recl <- sst.all
sst.all.recl[sst.all.recl < 0] <- 1
sst.all.recl[sst.all.recl > 0] <- 1

# False positives
# M3 = 10%
length(sst.all.recl[,3][sst.all.recl[,3] == 1 & sst.all.recl[,4] == 0]) * 100/dim(sst.all.recl)[1]
# M2 = 15.8%
length(sst.all.recl[,2][sst.all.recl[,2] == 1 & sst.all.recl[,4] == 0]) * 100/dim(sst.all.recl)[1]
# M1 = 13.3%
length(sst.all.recl[,1][sst.all.recl[,1] == 1 & sst.all.recl[,4] == 0]) * 100/dim(sst.all.recl)[1]

# False negatives
# M3 = 10%
length(sst.all.recl[,3][sst.all.recl[,3] == 0 & sst.all.recl[,4] == 1]) * 100/dim(sst.all.recl)[1]
# M2 = 15.8%
length(sst.all.recl[,2][sst.all.recl[,2] == 0 & sst.all.recl[,4] == 1]) * 100/dim(sst.all.recl)[1]
# M1 = 13.3%
length(sst.all.recl[,1][sst.all.recl[,1] == 0 & sst.all.recl[,4] == 1]) * 100/dim(sst.all.recl)[1]

# Identify range shifting species -------

Species_lat_ranges <- GBR_fish_site_P %>% 
  #filter(SPECIES_NAME %in% spp.frequency$SPECIES_NAME) %>%
  group_by(SPECIES_NAME) %>%
  summarize(min.lat.pre = min(SiteLat[P.1 == 1 & Pre.or.post.bleach == "Pre"], na.rm = T),
            max.lat.pre = max(SiteLat[P.1 == 1 & Pre.or.post.bleach == "Pre"], na.rm = T),
            min.lat.post = min(SiteLat[P.1 == 1 & Pre.or.post.bleach == "Post"], na.rm = T),
            max.lat.post = max(SiteLat[P.1 == 1 & Pre.or.post.bleach == "Post"], na.rm = T)) %>%
  filter(is.finite(min.lat.pre))

Species_lat_ranges$lat.ext.pre <- with(Species_lat_ranges, max.lat.pre - min.lat.pre)
Species_lat_ranges$lat.ext.post <- with(Species_lat_ranges, max.lat.post - min.lat.post)
Species_lat_ranges$lat.ext.change <- with(Species_lat_ranges, lat.ext.post - lat.ext.pre)

length(Species_lat_ranges$lat.ext.change[Species_lat_ranges$lat.ext.change < -1]) #25 species with range contraction (9.9%)
length(Species_lat_ranges$lat.ext.change[Species_lat_ranges$lat.ext.change > 1]) #29 species with range extension (11.5%)

Species_lat_ranges$lat.mp.pre <- with(Species_lat_ranges, (min.lat.pre + max.lat.pre)/2)
Species_lat_ranges$lat.mp.post <- with(Species_lat_ranges,(min.lat.post + max.lat.post)/2)
Species_lat_ranges$lat.mp.change <- with(Species_lat_ranges, lat.mp.post - lat.mp.pre)

length(Species_lat_ranges$lat.mp.change[Species_lat_ranges$lat.mp.change < -1]) #22 species with southern midpoint displacement
length(Species_lat_ranges$lat.mp.change[Species_lat_ranges$lat.mp.change > 1]) #8 species with northern midpoint displacement

Species_lat_ranges$range.cont <- ifelse(Species_lat_ranges$lat.ext.change < -1, 1, 0)
Species_lat_ranges$range.ext <- ifelse(Species_lat_ranges$lat.ext.change > 1, 1, 0)
Species_lat_ranges$range.displ <- ifelse(Species_lat_ranges$lat.mp.change > 1 | Species_lat_ranges$lat.mp.change < -1, 1, 0)

Species_lat_ranges$range.change <- ifelse(Species_lat_ranges$range.ext == 1 | Species_lat_ranges$range.cont == 1 | Species_lat_ranges$range.displ == 1, 1, 0)

length(Species_lat_ranges$range.ext[Species_lat_ranges$range.ext ==1])/length(Species_lat_ranges$range.ext)
length(Species_lat_ranges$range.cont[Species_lat_ranges$range.cont ==1])/length(Species_lat_ranges$range.ext)
length(Species_lat_ranges$range.displ[Species_lat_ranges$range.displ ==1])/length(Species_lat_ranges$range.ext)
length(Species_lat_ranges$range.change[Species_lat_ranges$range.change ==1])/length(Species_lat_ranges$range.ext)

# Distribution of TSS and Kappa for all species vs. range-shifting species --------

stats <- rbind(data.frame(method = "M1", SPECIES_NAME = tss$SPECIES_NAME, tss = tss$m2, kappa = kappa$m2),
               data.frame(method = "M2", SPECIES_NAME = tss$SPECIES_NAME, tss = tss$m3, kappa = kappa$m3),
               data.frame(method = "M3", SPECIES_NAME = tss$SPECIES_NAME, tss = tss$m4, kappa = kappa$m4),
               data.frame(method = "M4", SPECIES_NAME = tss$SPECIES_NAME, tss = tss$m5, kappa = kappa$m5))

stats$tss[stats$tss < 0] <- 0.001
stats$kappa[stats$kappa < 0] <- 0.001

g.kappa <- ggplot(stats, aes(factor(method), kappa)) +
  #geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))
  geom_boxplot(notch = T) + xlab("Method") + ylab("Kappa")  + ylim(0,1) + ggtitle("All species")

g.kappa.range.change <- ggplot(stats[stats$SPECIES_NAME %in% Species_lat_ranges$SPECIES_NAME[Species_lat_ranges$range.change ==1],], aes(factor(method), kappa)) +
  #geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))
  geom_boxplot(notch = T) + xlab("Method") + ylab("Kappa")  + ylim(0,1) + ggtitle("Range-shifting species") 

g.tss <- ggplot(stats, aes(factor(method), tss)) +
      #geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))
      geom_boxplot(notch = T) + xlab("Method") + ylab("TSS") + ylim(0,.6)

g.tss.range.change <- ggplot(stats[stats$SPECIES_NAME %in% Species_lat_ranges$SPECIES_NAME[Species_lat_ranges$range.change ==1],], aes(factor(method), tss)) +
  #geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))
  geom_boxplot(notch = T) + xlab("Method") + ylab("TSS") + ylim(0,.6)

plot_grid(g.kappa, g.kappa.range.change, g.tss, g.tss.range.change, nrow = 2, ncol = 2)


g.biplot <- ggplot() +
  geom_point(data = tss, aes(x = m2, y = m5), color = "dimgrey")+
  xlab("TSS (M.2)") + ylab("TSS (M.4)") +
  geom_abline(aes(slope = 1, intercept = 0), color = "dimgrey") +
  geom_density_2d(data = tss, aes(x = m2, y = m5), color = "blue")+
  geom_density_2d(data = tss[tss$SPECIES_NAME %in% Species_lat_ranges$SPECIES_NAME[Species_lat_ranges$range.change ==1],], aes(x = m2, y = m5), color = "green")+
  #stat_density_2d(data = tss, aes(x = m2, y = m5, fill = after_stat(level)), alpha = .5, geom = "polygon")+
  coord_cartesian(xlim = c(0,1), ylim = c(0,1))

ggplot() +
  geom_point(data = tss, aes(x = m2, y = m4, color = log10(spp.frequency$Nocc)))+
  xlab("TSS (M.2)") + ylab("TSS (M.4)")



