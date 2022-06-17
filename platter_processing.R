### pollen platter analysis ###

#
library(lubridate)
library(dplyr)
library(raster)
library(terra)
library(NISTunits)


#sampler number
sampler_n <- 14

#time start
time_start <- mdy_hms("5/9/2022 15:49:01")

#time deployed
time_deploy <- mdy_hms("5/10/2022 13:19:00")

#time retreived
time_retreived <- mdy_hms("5/12/2022 14:47:00")

#angle at retreival 
angle_retreival <- 17.8 + 360

#programmed step time
step_time_min <- 65.717 
step_time_sec = step_time_min * 60
step_angle <- (12/516) * 365

platter_df <- data.frame(time_period = 1:48, timestep_start = rep(NA, 48), timestep_end = rep(NA,48),
                         angle_start = rep(NA, 48), angle_end = rep(NA, 48))


platter_df <- platter_df %>% 
  mutate(timestep_start = time_deploy + step_time_sec * time_period - step_time_sec,
         timestep_end = step_time_min - (as.numeric(difftime(time_deploy, time_start, units = "mins")) %% step_time_min) + #time left before turn on deploy
                        time_deploy + 
                        step_time_sec * time_period - step_time_min,
         angle_start = step_angle * time_period - step_angle + (90 - step_angle), #the blue line starts on the left hand side of the open slit, so the first slot is just under 90 degrees
         angle_end = step_angle * time_period + (90 - step_angle))


#scanned file
file_scanned <- "pp_scan_sampler_14_d220510"
#rotated_image_r <- raster::stack("C:/Users/dsk273/Desktop/pp_scan_sampler_14_d220510_c.tif")
rotated_image <- terra::rast("C:/Users/dsk273/Desktop/pp_scan_sampler_14_d220510_c.tif")
crs(rotated_image) <- NA
terra::ext(rotated_image) <- c(0, ncol(rotated_image), 0, nrow(rotated_image))
plot(rotated_image)
#plotRGB(rotated_image)

platter_centroid_x <- 4366
platter_centroid_y <- 4320
platter_dist_inner <- 1640
platter_dist_outer <- 3600

### divide image up into chunks and save each ##############################################################################
#outer points of the focal slot
outer_points_x <- platter_centroid_x + platter_dist_outer * cos(NISTdegTOradian(platter_df$angle_start))
outer_points_y <- platter_centroid_y + platter_dist_outer * sin(NISTdegTOradian(platter_df$angle_start))

#inner points of the focal slot
inner_points_x <- platter_centroid_x + platter_dist_inner * cos(NISTdegTOradian(platter_df$angle_start))
inner_points_y <- platter_centroid_y + platter_dist_inner * sin(NISTdegTOradian(platter_df$angle_start))

#create a spatial vector out of the points and visual check
pp_aoi <- rbind(c(inner_points_x[1], inner_points_y[1]), 
                c(inner_points_x[2], inner_points_y[2]), 
                c(outer_points_x[2], outer_points_y[2]), 
                c(outer_points_x[1], outer_points_y[1]))
pp_aoi_v <- vect(pp_aoi, type = "polygons")

#plot(pp_aoi_v, col='yellow')
plot(rotated_image) #plot(1:1)
plot(pp_aoi_v, col='yellow', add = TRUE)

#extract image
#test <- terra::extract(rotated_image, pp_aoi_v)
test <- terra::crop(rotated_image, pp_aoi_v)
plot(test)
test2 <- terra::mask(test, pp_aoi_v, updatevalue = 1)
plot(test2)

terra::writeRaster(test2, "C:/Users/dsk273/Desktop/pp_scan_sampler_14_d220510_c_chunk1_b.tif")

#focal species pollen color

#run FIJI/imageJ macro for: remove background

#remove non-pollen colors and replace with a mask

#calculate area of mask

#remove background

#thresholding

#calculate area of non-mask

#segmentation

#calculate number of grains in image segment

#divide grains by non-masked area

#relative grains within platter sample
