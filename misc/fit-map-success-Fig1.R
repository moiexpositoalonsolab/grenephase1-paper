
################################################################################
### Goal
### Fig 1 with survival
################################################################################
library(tidyverse)
library(raster)
library(RColorBrewer)
library(ggplot2)
library(cowplot)

################################################################################
#### Location
# myfolder<-"/content/drive/Shareddrives/MOI-LAB/grenephase1-analyses" # GOOGLE DRVIE
# myfolder<-"~/grenephase1-analyses"
myfolder<-"~/Shareddrives//MOI-LAB/PROJECTS/grenenet/GRENE-net_PHASE1/grenephase1-analyses/" # GOOGLE DRVIE
setwd(myfolder)

################################################################################
## Survival data
sur=read.csv("data/survival.csv")

# ## Survival predicted
# popsize=read.csv("pop_size_estimation/pop_size_estimations.csv")
# popsizeraw=popsize

## All sites we started with
load("grene/data/locations_data.rda")
locations_data %>% head

sitessimplenames<-read.csv("grene/data/sites_simple_names.csv")
head(sitessimplenames)
dim(sitessimplenames)

narrative<-read_xlsx("../GrENE-net_manuscript/Supplemental_tables/Supplemental_Table_site_information.xlsx")

## Get accession locations
load("grene/data/ecotypes_data.rda")
ecotypes_data %>% head

################################################################################
# Get at least survival 1 replicate for

#quickcly look 
sur %>% 
  merge(.,sitessimplenames) %>% 
  dplyr::select(city, country, ends_with('survival')) %>% 
  view

# Look into who survived
maxsur<-
  sur %>% 
  group_by(site) %>% 
  dplyr::summarise(maxsur= any(X1_survival==1)) 

table(maxsur$maxsur)


maxsur<-
  sur %>% 
  group_by(site) %>% 
  dplyr::summarise(maxsur= any(X2_survival==1)) 

table(maxsur$maxsur)


maxsur<-
  sur %>% 
  group_by(site) %>% 
  dplyr::summarise(maxsur= any(X3_survival==1)) 

table(maxsur$maxsur)


maxsur<-
  sur %>% 
  group_by(site) %>% 
  dplyr::summarise(maxsur1= as.numeric(any(X1_survival==1)),
                   maxsur2= as.numeric(any(X2_survival==1)),
                   maxsur3= as.numeric(any(X3_survival==1))
                   ) %>% 
  dplyr::mutate(yearslong=maxsur1+maxsur2+maxsur3)

maxsur<-
maxsur %>% 
  merge(.,locations_data)


################################################################################
# Map
# moiR::ggplot_world_map()
# 
# ggplot_world_map<-
# function (mapcol = "seashell3", backcol = "white", countrycol = mapcol, 
#           projection = "cartesian", xlim = c(-200, +200), ylim = c(-60, 
#                                                                    +80), orientation = c(80, 0, 0)) 
# {
#   library(ggplot2)
#   library(cshapes)
#   # library(gpclib)
#   # gpclibPermit()
#   world <- cshp(date = as.Date("2008-1-1"))
#   world.points <- fortify(world, region = "COWCODE")
#   xticks = round(seq(xlim[1], xlim[2], length.out = 10), digits = 0)
#   yticks = round(seq(ylim[1], ylim[2], length.out = 10), digits = 0)
#   p1 <- ggplot()
#   p1 <- p1 + geom_polygon(data = world.points, aes(long, lat, 
#                                                    group = group), fill = mapcol, col = countrycol)
#   if (projection == "cartesian") {
#     p1 <- p1 + coord_cartesian(xlim = xlim, ylim = ylim)
#   }
#   else if (projection == "mercator") {
#     p1 <- p1 + coord_map("mercator", xlim = xlim, ylim = ylim)
#   }
#   else if (projection == "perspective") {
#     p1 <- p1 + coord_map("perspective", dist = 1, xlim = xlim, 
#                          ylim = ylim)
#   }
#   else if (projection == "ortho") {
#     p1 <- p1 + coord_map("ortho", xlim = xlim, ylim = ylim, 
#                          orientation = orientation)
#   }
#   else {
#     stop("Unknown projection")
#   }
#   p1 <- p1 + xlab("") + ylab("")
#   p1 <- p1 + scale_x_continuous(breaks = xticks, labels = paste(abs(xticks), 
#                                                                 ifelse(xticks < 0, "W", "E"), sep = " "))
#   p1 <- p1 + scale_y_continuous(breaks = yticks, labels = paste(yticks, 
#                                                                 "N", sep = " "))
#   p1 <- p1 + theme(panel.grid.major = element_line(colour = "white", 
#                                                    size = 0, linetype = "dashed"))
#   p1 <- p1 + theme(panel.grid.minor = element_line(colour = "white", 
#                                                    size = 0, linetype = "dashed"))
#   p1 <- p1 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
#                    panel.background = element_rect(fill = backcol, colour = "white"), 
#                    panel.border = element_rect(colour = "black", fill = NA), 
#                    legend.title = element_blank(), legend.background = element_blank(), 
#                    legend.key = element_blank())
#   p1
#   return(p1)
# }


world <- map_data("world")
worldplot <- ggplot() +
  geom_polygon(data = world, aes(x=long, y = lat, group = group),fill = "seashell3",col = "seashell3") + 
  coord_fixed(1.3)
worldplot

xlim = c(-200, +200)
ylim = c(-60,+80)
orientation = c(80, 0, 0)

mapcol = "seashell3"
backcol = "white"


p1 <- worldplot +
  coord_map("ortho", xlim = xlim, ylim = ylim, orientation = orientation) +
  xlab("") + ylab("") +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),  # Remove the background
    panel.border = element_blank(),  # Remove the border
    axis.line = element_blank(),  # Remove the axis line
    plot.background = element_blank(),  # Remove outer plot background
    legend.title = element_blank(),
    legend.background = element_blank(),
    legend.key = element_blank()
  )


################################################################################
survivalcolors<-c("-1"="#8A181A","0"="#FAA51A", "1"="#D9E7A4","2"= "#7AC479","3"= "#025B32")
# Add the survival of sites
finalplot<-
p1+
  # Ecotypes
  geom_point(data=ecotypes_data, aes(y=latitude,x=longitude), shape=3)+
  # Sites
  geom_point(data=narrative,aes(y=LATITUDE,x=LONGITUDE),
             size=3.2, alpha=1, color="black", shape=16
  )+
  geom_point(data=narrative,aes(y=LATITUDE,x=LONGITUDE),
             size=3, alpha=1, color="white", shape=16
  )+
  geom_point(data=narrative,aes(y=LATITUDE,x=LONGITUDE, color=factor(SURVIVALYEARS)),
             size=2.5, alpha=0.8, shape=16
             )+
  scale_color_manual(values = survivalcolors,
                     labels=c(
                       "Logistical fail",
                       "No survival",
                       "Survival 1 year",
                       "Survival 2 year",
                       "Survival 3 year"
                     )
                     )

finalplot

save_plot("figs/fig-fig1-map-surival.pdf",finalplot, base_height = 7, base_width = 7)
save_plot("figs/fig-fig1-map-surival.png",finalplot, base_height = 7, base_width = 7)
