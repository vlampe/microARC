# extract chl a values from satellite observations alsong model trajectories
setwd("Documents/microARC model/")
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


# load chl data
# data source:
# https://coastwatch.pfeg.noaa.gov/erddap/griddap/nesdisVHNSQchlaMonthly.nc?chlor_a%5B(2018-01-01T12:00:00Z):1:(2018-12-31T12:00:00Z)%5D%5B(0.0):1:(0.0)%5D%5B(89.75625):1:(50)%5D%5B(-179.9812):1:(179.9813)%5D
# chl_nc <- brick("nesdisVHNSQchlaMonthly_897a_a624_cd26.nc", values = T) 
# chl_nc <- crop(chl_nc, trajSPDF)


## andere Daten benutzen:
# https://centaur.reading.ac.uk/86253/
# https://oceancolour.org/thredds/ncss/grid/CCI_ALL-v5.0-MONTHLY/dataset.html
# the european chl data product
# download link https://oceancolour.org/thredds/ncss/CCI_ALL-v5.0-MONTHLY?var=chlor_a&north=89.9791&west=-10&east=40&south=65&disableProjSubset=on&horizStride=1&time_start=2018-01-01T00%3A00%3A00Z&time_end=2018-12-31T00%3A00%3A00Z&timeStride=1&addLatLon=true
chl_nc <- brick("CCI_ALL-v5.0-MONTHLY.nc", values = T)
chl_nc <- crop(chl_nc, trajSPDF)

# check extent and values
plot(chl_nc[[5]], zlim = c(0,10))


Trajs <- unique(trajCoords$idTraj)
Months <- unique(trajCoords$monthTraj)

# # for testing
# id <- 70
# id <- 1820
# id <- 883
# m <- 7

monthlyChl <- tibble("Traj" = numeric(), "month" = numeric(), 
                     "mean_chl" = numeric(), "month_chl" = numeric())

for(id in Trajs) {
  
  for(m in Months) {
    # subset input data
    chl_month <- chl_nc[[m]]
    coords_sub <- filter(trajCoords, idTraj == id, monthTraj == m)
     # put coordinates into SpatialPoints-Class
    coords_pts <- SpatialPoints(data.frame(x = coords_sub$lonTraj, y = coords_sub$latTraj))
    crs(coords_pts) <- crs(chl_month)
    
    # create buffer around points
    coords_buffer <- buffer(coords_pts, width = 20000) # width ist distanz in m 
    
    # plot
    ex_coords <- extent(coords_pts)
    # plot(chl_month, xlim = c(ex_coords@xmin, ex_coords@xmax), 
    #      ylim = c(ex_coords@ymin, ex_coords@ymax))
    ## Loop funktioniert, aber wäre cool wenn gute grafiken mit produziert werden
    # extent anpassen und titel hinzufügen
    plot(chl_month)
    plot(coords_buffer, add = T, col = NA)
    points(coords_pts, add = T, cex = .1)
    title(paste0("Trajectory ", as.character(id), " and 20 km radius in month ",
                 as.character(m)))
    
    # use coords_buffer as a clipping mask on chl_month
    chl_coords <- mask(chl_month, coords_buffer)
    plot(chl_coords)
    title(paste0("Chlorophyll in trajectory ", as.character(id), " and it's 20 km radius in month ",
                 as.character(m)))
    
    # get mean of chl_coords
    meanTrajChl <- mean(getValues(chl_coords), na.rm = T)
    meanMonthChl <- mean(getValues(chl_month), na.rm = T)
    
    monthlyChl <- bind_rows(monthlyChl,
                     tibble("Traj" = id, "month" = m, "mean_chl" = meanTrajChl,
                     "month_chl" = meanMonthChl))
    rm(meanTrajChl, meanMonthChl)
  }
}

# save output
write_csv(monthlyChl, path = "monthlyTrajChlAutumn.csv")

# plot chl concentration vs time, for each traj

monthlyChl %>%
  mutate(id = as.factor(Traj)) %>%
  group_by(id) %>%
  ggplot(aes(x = month, y = mean_chl, color = id)) +
    geom_line()+
    geom_line(aes(y = month_chl), color = "black")+ # mean chl of whole extent of all trajectories
    scale_x_continuous(limits = c(1, 10), breaks = (c(1:10)), labels = c(1:10)) +
    labs(title = "Chlorophyll (remote sensing) concentrations along model \ntrajectories (20 km radius)",
         y = "Chlorophyll concentration [mg m^-3]",
         colour = "Trajectory ID",
         caption = "black line indicates mean Chl concentration within -34.2,41.3, 66.2, 85.0") +
    theme(legend.key.size = unit(3, "points"),
          text = element_text(size = 9))

ggsave(filename = "monthlyTrajChl.png", width = 11, height = 7, units = "cm")


# plot with color showing watermass
wmTrajs <- unique(dplyr::select(trajCoords, idTraj, wmTraj))
wmTrajs$wmTraj <- factor(wmTrajs$wmTraj, levels = c("Atlantic", "Arctic"))
monthlyChl %>%
  right_join(wmTrajs, by = c("Traj" = "idTraj")) %>%
  mutate(id = as.factor(Traj),
         wmTraj = as.factor(wmTraj)) %>%
  group_by(id, wmTraj) %>%
  ggplot(aes(x = month, y = mean_chl, color = wmTraj, group = id)) +
  geom_line()+
  geom_line(aes(y = month_chl), color = "black")+ # mean chl of whole extent of all trajectories
  scale_x_continuous(limits = c(1, 10), breaks = (c(1:10)), labels = c(1:10)) +
  labs(title = "Chlorophyll (remote sensing) concentrations along model \ntrajectories (20 km radius)",
       subtitle = "by trajectory origin",
       y = "Chlorophyll concentration [mg m^-3]", 
       colour = "origin",
       caption = "black line indicates mean Chl concentration within -34.2,41.3, 66.2, 85.0") +
  theme(#legend.key.size = unit(3, "points"),
        text = element_text(size = 9))
ggsave(filename = "monthlyTrajChl_origin.png", width = 11, height = 7, units = "cm")


# 
# 
# # create spatial lines from coords_sub, and create buffer polygon around it
# l <- Line(cbind(coords_sub$lonTraj, coords_sub$latTraj))
# ls <- Lines(list(l), 1)
# sls <- SpatialLines(list(ls))
# proj4string(sls) <-  crs(chl_month)
# 
# #plot(chl_month)
# #plot(sls, add = T)
# 
# buffer_traj <- gBuffer(sls, capStyle = "FLAT", joinStyle = "MITRE", width = .1)
# 
# 
# 
# #plot(buffer_traj, add = T)
# buffer_tst3 <- buffer(sls, width = 1) ## wie weit ist das? in pixel? km ist unpassend weil karte nicht projeziert ist
# # buffer scheint nicht mit spatial lines klar zu kommen. 
# plot(chl_month)
# plot(buffer_tst3, add = T, col = "red")
# plot(sls, add = T)
# 
# 
# chl_tst <- chl_nc[[5]]
# traj_tst <- filter(trajCoords, idTraj == 70, monthTraj == 5)
# traj_tst <- filter(trajCoords, idTraj == 70)
# 
# # make a polygon with a buffer around the trajectory
# 
# l <- Line(cbind(traj_tst$lonTraj, traj_tst$latTraj))
# ls <- Lines(list(l), 1)
# sls <- SpatialLines(list(ls))
# plot(sls)  
# proj4string(sls) <-  crs(chl_tst)
# 
# plot(chl_tst)
# plot(sls, add = T)
# 
# ex <- extent(-30, 60, 65, 85)
# chl_ex <- crop(chl_tst, ex)
# chl_ex2 <- crop(chl_ex, trajSPDF)
# 
# plot(chl_ex2)
# 
# plot(chl_ex)
# plot(sls, add = T)
# 
# # st_buffer can draw polygon arounf SpatialLines
# coord_sf <- st_sf(sls)
# buffer_tst <- st_buffer(sls, dist = 10000) # km buffer around trajectory
# # gBuffer
# buffer_tst2 <- gBuffer(sls)
# 
# 
# 
# plot(chl_ex)
# plot(sls, add = T)
# plot(buffer_tst2, add = T)
# 
# chl_ex_crop <- crop(chl_ex, buffer_tst2)
# plot(chl_ex_crop)
# plot(buffer_tst2, add = T)
# 
# chl_ex_mask <- mask(chl_ex, buffer_tst2)
# plot(chl_ex_mask)
# 
# mean(getValues(chl_ex), na.rm = T)
# mean(getValues(chl_ex_mask), na.rm = T)
# 
# 
# # plot(chl_nc[[5]])
# # 
# # 
# # tst <- filter(trajCoords, idTraj == 70)
# # lines(x = tst$lonTraj, y = tst$latTraj)
# # ## maybe convert coordinates to SpatialPolygon or SpatialLine first
# # # then use extract to crop raster to trajectory path
# # 
# # 
# # # put coordinates into SpatialPoints-Class
# # pts <- SpatialPoints(data.frame(x = tst$lonTraj, y = tst$latTraj))
# # points(pts)
# # 
# # # # lns <- as.SpatialLines(tst)
# # # lines_tst <- cbind(tst$lonTraj, tst$latTraj)
# # # lines_tst2 <- Line(lines_tst)
# # # lines_tst3 <- SpatialLines(lines_tst2)
# # 
# #  #tstdf <- SpatialLinesDataFrame(tst)
# # # Use the points directly with the buffer argument
# # # extract values within buffer of pts from chl_nc
# # tst2 <- raster::extract(chl_nc[[5]], pts, buffer = 20000)
# # plot(tst2[[1]])
# # # extract with cell numbers
# # tst4 <- raster::extract(chl_nc[[5]], pts, buffer = 20000, cellnumbers = T)
# # # get coordinates from cell numbers
# # tst5 <- lapply(tst4, function(x) xyFromCell(chl_nc[[5]], x[,1]))
# # 
# # tst6 <- rasterize(pts, chl_nc[[5]])
# # plot(tst6)
# # 
# # 
# # 
# # # Let's compute the average value for each buffer
# # tst3 <- sapply(tst2, mean)
# # plot(tst3)
# 
# nc_croppedto_pts <- raster::mask(chl_nc[[5]], pts, inverse = T)
# plot(nc_croppedto_pts)
# lines(x = tst$lonTraj, y = tst$latTraj)
# # 

