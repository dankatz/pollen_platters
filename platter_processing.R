### pollen platter analysis ###

#
library(lubridate)
library(dplyr)
library(readr)
library(raster)
library(terra)
library(NISTunits) #install.packages("NISTunits")
library(here)
library(purrr)
library(ggplot2)

setwd("C:/Users/danka/Box")
here::i_am("katz_photo.jpg")


#sampler number
sampler_n <- "9"

#time start
time_start <- mdy_hms("6/21/2022 17:15:00")

#time deployed
time_deploy <- mdy_hms("6/24/2022 14:59:00")

#time retreived
time_retreived <- mdy_hms("6/26/2022 13:57:00")

#angle at retreival 
angle_retreival <- 8 + 360

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
         angle_end = step_angle * time_period + (90 - step_angle),
         seg_obj_n = NA,
         seg_obj_size = NA)


#scanned file
# file_scanned <- "pp_scan_sampler9_d220624_r.tif"
# file_path <- "C:/Users/dsk273/Box/Cornell/mentoring/student projects/summer 2022/Kent pollen catcher/platter_scans/processed_scan/"
#here("Cornell", "mentoring", "student projects", "summer 2022", "Kent pollen catcher", "platter_scans", "processed_scan", "pp_scan_sampler9_d220624_r.tif")
#rotated_image_r <- raster::stack("C:/Users/dsk273/Desktop/pp_scan_sampler_14_d220510_c.tif")
rotated_image <- terra::rast(here("Cornell", "mentoring", "student projects", "summer 2022", "Kent pollen catcher", "platter_scans", "processed_scan", 
                                  "pp_scan_sampler9_d220624_r.tif"))
#rotated_image <- terra::rast(paste0(file_path, file_scanned))
crs(rotated_image) <- NA
terra::ext(rotated_image) <- c(0, ncol(rotated_image), 0, nrow(rotated_image))
plot(rotated_image)
#plotRGB(rotated_image)

platter_centroid_x <- 4047
platter_centroid_y <- 4044
platter_dist_inner <- 1640
platter_dist_outer <- 3600

### divide image up into chunks and save each ##############################################################################
platter_df <- platter_df %>% 
  mutate(
    outer_points_x_start = platter_centroid_x + platter_dist_outer * cos(NISTdegTOradian(platter_df$angle_start)), #outer points of the focal slot
    outer_points_y_start = platter_centroid_y + platter_dist_outer * sin(NISTdegTOradian(platter_df$angle_start)),
    inner_points_x_start = platter_centroid_x + platter_dist_inner * cos(NISTdegTOradian(platter_df$angle_start)),#inner points of the focal slot
    inner_points_y_start = platter_centroid_y + platter_dist_inner * sin(NISTdegTOradian(platter_df$angle_start)),
    outer_points_x_end   = platter_centroid_x + platter_dist_outer * cos(NISTdegTOradian(platter_df$angle_end)), #outer points of the focal slot
    outer_points_y_end   = platter_centroid_y + platter_dist_outer * sin(NISTdegTOradian(platter_df$angle_end)),
    inner_points_x_end   = platter_centroid_x + platter_dist_inner * cos(NISTdegTOradian(platter_df$angle_end)),#inner points of the focal slot
    inner_points_y_end  = platter_centroid_y + platter_dist_inner * sin(NISTdegTOradian(platter_df$angle_end)))

#loop through each individual time chunk and save a file for it
for(i in 1:4){ #max 42
  #create a spatial vector out of the points and visual check
  pp_aoi <- rbind(c(platter_df$outer_points_x_start[i], platter_df$outer_points_y_start[i]), 
                  c(platter_df$inner_points_x_start[i], platter_df$inner_points_y_start[i]), 
                  c(platter_df$inner_points_x_end[i], platter_df$inner_points_y_end[i]), 
                  c(platter_df$outer_points_x_end[i], platter_df$outer_points_y_end[i]))
  pp_aoi_v <- vect(pp_aoi, type = "polygons")
   plot(rotated_image) #plot(1:1)
  plot(pp_aoi_v, col='yellow', add = TRUE)
  
  #extract image
  #test <- terra::extract(rotated_image, pp_aoi_v)
  chunk_x_uncropped <- terra::crop(rotated_image, pp_aoi_v)
  #plot(chunk_x_uncropped)
  chunk_x <- terra::mask(chunk_x_uncropped, pp_aoi_v, updatevalue = 1)
  #plot(chunk_x)

  file_save_name <- paste0(gsub(x = file_scanned, pattern = ".tif", replacement = ""), "_chunk_", i, ".tif")
  terra::writeRaster(chunk_x, 
                     
                     here("Cornell", "mentoring", "student projects", "summer 2022", "Kent pollen catcher", "platter_scans", "platter_chunks", 
                          file_save_name),
                     overwrite = TRUE)
}




### read in results from FIJI/Labkit ##############################################################################################
#run the macro in FIJI:
#C:\Users\danka\Box\Cornell\mentoring\student projects\summer 2022\Kent pollen catcher\Labkit_classifications
#Macro_labkit_plantain2.ijm.ijm

result_csvs <- dir(here("Cornell", "mentoring", "student projects", "summer 2022", "Kent pollen catcher", "Labkit_classifications", "classification_chunk_results"), full.names = TRUE) %>%  
                map_dfr(.x = ., .f = read_csv)
result_csvs <- result_csvs %>% 
              mutate(time_period = gsub(pattern = ".*_", replacement = "", x = Slice),
                     time_period = as.numeric(gsub(pattern = ".tif", replacement = "", x = time_period))) %>% 
              arrange(time_period)

result_csvs %>% 
  ggplot(aes(x = time_period, y = Count)) + geom_point() + geom_line() + theme_bw()

pol_dep <- platter_df %>% dplyr::select(time_period, timestep_start, timestep_end) %>% 
            left_join(., result_csvs)

pol_dep %>% 
  ggplot(aes(x = timestep_start, y = Count)) + geom_point() + geom_line() + theme_bw()


