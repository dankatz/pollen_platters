### pollen platter analysis ###

#
library(lubridate)
library(dplyr)
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

platter_df <- data.frame(period = 1:48, timestep_start = rep(NA, 48), timestep_end = rep(NA,48),
                         angle_start = rep(NA, 48), angle_end = rep(NA, 48))


platter_df <- platter_df %>% 
  mutate(timestep_start = time_deploy + step_time_sec * period - step_time_sec,
         timestep_end = step_time_min - (as.numeric(difftime(time_deploy, time_start, units = "mins")) %% step_time_min) + #time left before turn on deploy
                        time_deploy + 
                        step_time_sec * period - step_time_min,
         angle_start = step_angle * period - step_angle,
         angle_end = step_angle * period)


#scanned file
file_scanned <- "pp_scan_sampler_14_d220510"
rotated_image <- rast("C:/Users/dsk273/Desktop/pp_scan_sampler_14_d220510_c.tif")

terra::plot(rotated_image)
platter_centroid_x <- 4366
platter_centroid_y <- 4320
platter_dist_inner <- 1640
platter_dist_outer <- 3600

#divide image up into chunks
platter_outer_left_x <- platter_centroid_x + platter_dist_outer
platter_outer_left_y <- platter_centroid_y

platter_outer_right_x <- platter_centroid_x

S <- 2 * platter_dist_outer * sin(platter_df$angle_start[2]/2)


points <- 100
radius <- 10
center_x <- 5
center_y <- 5

drawCirclePoints <- function(points, radius, platter_centroid_x, platter_centroid_y) {
  slice <- 2 * pi / points
  angle <- slice * seq(0, points, by = 1)
  
  newX <- platter_centroid_x + radius * cos(angle)
  newY <- platter_centroid_y + radius * sin(angle)
  
  plot(newX, newY)
print(newX)
}

drawCirclePoints(points, radius = platter_dist_outer, platter_centroid_x, platter_centroid_y)



testx <- platter_centroid_x + platter_dist_outer * cos(NISTdegTOradian(platter_df$angle_start))
testy <- platter_centroid_y + platter_dist_outer * sin(NISTdegTOradian(platter_df$angle_start))

terra::plot(rotated_image)
plot(testx, testy)


# create a point shapefile from the dataframe, then convert to polygon
polygon <- df %>%
  st_as_sf(coords = c("x", "y"), crs = crs(rotated_image)) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POINT") %>% 
  st_buffer(dist=200)




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
