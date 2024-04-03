# extract data from World Ocean Atlas along model trajectories
# nitrogen, salinity, temperature


# setwd("Documents/microARC model/")
library(tidyverse)
library(raster)
library(sp)
require(rgeos)
library(sf)

# load trajectories
trajCoords <- read_csv("trajCoordsAutumn.csv")

xy = trajCoords[c("lonTraj", "latTraj")]
trajSPDF <- trajCoords
coordinates(trajSPDF) <- xy

Trajs <- unique(trajCoords$idTraj)
Months <- unique(trajCoords$monthTraj)

# Documentation of WOA data
# https://www.ncei.noaa.gov/data/oceans/woa/WOA18/DOC/woa18documentation.pdf

# specify general path to the data (here WOA18 data at server of NOAA, USA) 
path_WOA18 <- "https://www.ncei.noaa.gov/data/oceans/woa/WOA18/DATA"
# path where file should be saved
save_path  <- "./WOA data/" 

# specify files for download
urls  <- paste0(path_WOA18,
               "/nitrate/netcdf/all/1.00/woa18_all_n",
               sprintf("%02d", c(1:12)), "_01.nc")
files <- paste0(save_path,
               "monthly_WOA18_nitrate_",
               sprintf("%02d", c(1:12)), ".nc")
# # now downloading files from the URLs provided
# download.file(urls, files)

tst <- brick(paste0(save_path, "monthly_WOA18_nitrate_06.nc"), varname = "n_an")
# layers (z dimension) discribe depths
tst@data@unit
Depths <- getZ(tst)

monthlyWOA_N <- tibble("Traj" = numeric(), "wmTraj" = character(), 
                       "month" = numeric(), "depth" = numeric(), 
                       "mean_N" = numeric())
## for testing
# id <- 70
# m <- 7
# di <- 1

pb = txtProgressBar(min = 0, max = length(Trajs), initial = 0) 
stepi = 0

for (id in Trajs) {
  idOrigin <- unique(trajCoords$wmTraj[trajCoords$idTraj == id])
  stepi = stepi + 1
  for (m in Months) {
    
    WOA_N.nc <- brick(files[m], varname = "n_an")
    WOA_N.nc <- crop(WOA_N.nc, trajSPDF)
    coords_sub <- filter(trajCoords, idTraj == id, monthTraj == m)
    # put coordinates into SpatialPoints-Class
    coords_pts <- SpatialPoints(data.frame(x = coords_sub$lonTraj, 
                                           y = coords_sub$latTraj))
    crs(coords_pts) <- crs(WOA_N.nc)
    
    # create buffer around points
    coords_buffer <- buffer(coords_pts, width = 20000)
    
    for (di in 1:length(Depths)) {
      # ...
      WOA_N_d <- WOA_N.nc[[di]]
      d <- getZ(WOA_N.nc[[di]])
      
      # plot(WOA_N_d)
      # plot(coords_buffer, add = T, col = NA)
      # points(coords_pts, add = T, cex = .1)
      # title(paste0("Trajectory ", as.character(id), " and 20 km radius in month ",
      #              as.character(m)))
      
      # plot(crop(WOA_N_d, coords_buffer))
      # plot(coords_buffer, add = T, col = NA)
      
      ## this method keeps only cells that have their majority within the polygon. 
      ## this results in NA values for small polygons (ie slower parts of the trajectories)
      # # use coords_buffer as a clipping mask on WOA_N_d
      # N_coords <- mask(WOA_N_d, coords_buffer)
      
      # alternative: extract all cells intersecting the polygons
      cls <- cellFromPolygon(WOA_N_d, coords_buffer, weights = T)[[1]][, "cell"]
      N_coords <- WOA_N_d
      N_coords[][-cls] <- NA
      # plot(N_coords)
      # plot(trim(N_coords))
      # plot(coords_buffer, add = T)
      
      # plot(N_coords)
      # title(paste0("Nitrogen in trajectory ", as.character(id), " and it's 20 km radius in month ",
      #              as.character(m), " at ", as.character(d), "m depth"))
      # 
      # get mean of N_coords
      meanTrajN <- mean(getValues(N_coords), na.rm = T)
      # meanMonthN <- mean(getValues(N_month), na.rm = T)
      
      # save restults
      monthlyWOA_N <- bind_rows(monthlyWOA_N,
                       tibble("Traj" = id, "wmTraj" = idOrigin,
                              "month" = m, "depth" = d, 
                              "mean_N" = meanTrajN))
    }
  }
  setTxtProgressBar(pb,stepi)
}



# tabelle:
# traj | wmOrigin | month | depth | mean var

# save output
write_csv(monthlyWOA_N, path = "monthlyTrajNAutumn.csv")

# load ouput
monthlyWOA_N <- read_csv("monthlyTrajN0.csv", col_types = c("ncnnn"))
# plot N concentration vs time, for each traj and depth
# HIER WEITER
# oben noch mal kontrollieren, plots fertig und speichern, dann vergleich mit modellergebnissen? 



# plot by id (linetype) and depth (colour)
monthlyWOA_N %>%
  mutate(id = as.factor(Traj), 
         depth = as.factor(depth)) %>%
  ggplot(aes(x = month, y = mean_N, color = depth, linetype = id)) +
   geom_line()+
   scale_x_continuous(limits = c(1, 10), breaks = (c(1:10)), labels = c(1:10)) +
   labs(title = "Nitrogen (WOA18) concentrations along model \ntrajectories (20 km radius)",
        y = "Nitrogen concentration [µmol kg^-1]",
        colour = "Depth [m]",
         linetype = "Traj ID") +
    theme(legend.key.size = unit(3, "points"),
         text = element_text(size = 9))

ggsave(filename = "monthlyTrajN.png", width = 11, height = 7, units = "cm")

# plot by watermass (linetype) and depth (colour)

summarised_N <- 
  monthlyWOA_N %>%
  group_by(wmTraj, month, depth) %>% #, wmTraj, depth
    summarise(mean = mean(mean_N, na.rm = T)) %>%
    ungroup() 
 
summarised_N %>%
  ggplot(aes(x = month, y = mean, color = as.factor(depth), linetype = as.factor(wmTraj))) +
    geom_line()+
    scale_x_continuous(limits = c(1, 10), breaks = (c(1:10)), labels = c(1:10)) +
    labs(title = "Nitrogen (WOA18) concentrations along model \ntrajectories (20 km radius)",
         subtitle = "by trajectory origin",
         y = "Nitrogen concentration [µmol kg^-1]",
         colour = "Depth [m]",
         linetype = "origin") +
  scale_color_viridis_d() +
    theme(legend.key.size = unit(3, "points"),
          text = element_text(size = 9))
 
summarised_N %>%
  #filter(wmTraj == "Arctic") %>%
  ggplot(aes(x = month, y = depth, z = mean)) +
    scale_y_reverse() +
    scale_x_continuous(breaks = c(1:10)) +
    # geom_raster(aes(fill = mean))
    geom_tile(aes(fill = mean)) +
    scale_fill_viridis_c()  +
    facet_wrap(~wmTraj) +
    labs(title = "Nitrogen concentration along model \ntrajectories (20 km radius)",
         subtitle = "averaged by trajectory origin", 
         y = "Depth [m]", x = "month", fill = "Nitrogen \nconcentration \n[µmol kg^-1]")
 
summarised_N %>%
  #filter(wmTraj == "Arctic") %>%
  ggplot(aes(x = month, y = depth, color = mean)) +
    scale_y_reverse() +
    scale_x_continuous(breaks = c(1:10)) +
    scale_color_viridis_c(aesthetics = c("color", "fill")) +
    geom_point()+ 
    facet_wrap(~wmTraj) +
    labs(title = "Nitrogen concentration along model \ntrajectories (20 km radius)",
         subtitle = "averaged by trajectory origin", 
         y = "Depth [m]", x = "month", 
         color = "Nitrogen \nconcentration \n[µmol kg^-1]")
 
summarised_N %>%
  ggplot(aes(x = month, y = depth, z = mean)) +
  scale_y_reverse() +
  geom_contour_filled() +
  geom_contour(colour = "white", size = .3)+
  metR::geom_text_contour(colour = "white")+
  scale_x_continuous(breaks = c(1:10)) +
  scale_color_viridis_d(aesthetics = c("color", "fill")) +
  geom_point(size = .1)+ 
  facet_wrap(~wmTraj) +
  labs(title = "Nitrogen concentration along model \ntrajectories (20 km radius)",
       subtitle = "averaged by trajectory origin", 
       y = "Depth [m]", x = "month", 
       fill = "Nitrogen \nconcentration \n[µmol kg^-1]")
# this did not work before updates, therefore I also tried interpolation to create
# a section plot usin geom_contour_filled(). It looks better when the interpolation 
# is done in geom_contour_filled() (no extra interpolation step). 


# interpolation to fill data gaps
# for Arctic and Atlantic separately 
origdata_arc <- summarised_N %>%
  filter(wmTraj == "Arctic")%>%
  dplyr::select(-wmTraj)
   
grid_arc <- with(origdata_arc, interp::interp(month, depth, mean), extrap =T)
griddf_arc <- subset(data.frame(x = rep(grid_arc$x, nrow(grid_arc$z)),
                            y = rep(grid_arc$y, each = ncol(grid_arc$z)),
                            z = as.numeric(grid_arc$z)),
                 !is.na(z))

origdata_atl <- summarised_N %>%
  filter(wmTraj == "Atlantic")%>%
  dplyr::select(-wmTraj)

grid_atl <- with(origdata_atl, interp::interp(month, depth, mean), extrap =T)
griddf_atl <- subset(data.frame(x = rep(grid_atl$x, nrow(grid_atl$z)),
                                y = rep(grid_atl$y, each = ncol(grid_atl$z)),
                                z = as.numeric(grid_atl$z)),
                     !is.na(z))
 
griddf <- bind_rows("Arctic" = griddf_arc, "Atlantic" = griddf_atl, .id = "origin")
# create section plot

ggplot(griddf, aes(x, y, z = z)) +
  geom_contour_filled() +
  geom_point(data = summarised_N, aes(x = month, y = depth, z = NULL), size = .1) +
  scale_y_reverse() +
  scale_x_continuous(breaks = c(1:10)) +
  scale_color_viridis_d(aesthetics = c("color", "fill")) +
  facet_wrap(~ origin) +
  labs(title = "Nitrogen concentration along model \ntrajectories (20 km radius)",
       subtitle = "averaged by trajectory origin", 
       y = "Depth [m]", x = "month", 
       fill = "Nitrogen \nconcentration \n[µmol kg^-1]")

# noch mal das gleiche als tile plot
ggplot(griddf, aes(x, y, z = z)) +
  geom_tile(aes(fill = z)) +
  geom_point(data = summarised_N, aes(x = month, y = depth, z = NULL), size = .1) +
  geom_contour(colour = "white") +
  scale_y_reverse() +
  scale_x_continuous(breaks = c(1:10)) +
  scale_color_viridis_c(aesthetics = c("color", "fill")) +
  facet_wrap(~ origin) +
  labs(title = "Nitrogen concentration along model \ntrajectories (20 km radius)",
       subtitle = "averaged by trajectory origin", 
       y = "Depth [m]", x = "month", 
       fill = "Nitrogen \nconcentration \n[µmol kg^-1]")

