### pollen platter analysis ###

library(lubridate)
library(dplyr)
library(readr)
library(tidyr)
library(raster)
library(terra)
library(NISTunits) #install.packages("NISTunits")
library(here)
library(purrr)
library(ggplot2)
library(googlesheets4)
library(stringr)
library(sf)
beepr::beep(5)

setwd("C:/Users/danka/Box")
here::i_am("katz_photo.jpg")




### connect to deployment spreadsheet #######################################################################################
deployment_sheet <- read_csv(here("texas",  "pollen_platter", "pollen_platters_fs20_21_220809.csv"))

platters_to_process_list_rows <- c(1:110)
for(j in 1:length(platters_to_process_list_rows)){
  print(j)
  platter_row <- platters_to_process_list_rows[j] #platter_row <- platters_to_process_list_rows[1]
  scanned_file_name_focal <- deployment_sheet$scanned_file[platter_row]
  
  #sampler number
  sampler_n <- deployment_sheet$sampler[platter_row]  #"9"
  
  # #time start
  # time_start <- mdy_hm(deployment_sheet$sampler_start_date_time[platter_row])
  
  #time deployed
  time_deploy <- mdy_hm(deployment_sheet$date_deploy_auto[platter_row]) #mdy_hms("6/24/2022 14:59:00")
  
  #time retreived
  time_retreived <- mdy_hm(deployment_sheet$date_retreive_auto[platter_row]) #mdy_hms("6/26/2022 13:57:00")
  
  #angle at retreival 
  #angle_retreival <- deployment_sheet$retreival_angle[platter_row] + 90 #the measured angle is from the blue line to the left side of the opening slit. The angles used in this script start at 90 #8 + 360
  
  #platter centroid (manually measured)
  platter_centroid_x <- deployment_sheet$image_centroid_x[platter_row] # 4047 #
  platter_centroid_y <- deployment_sheet$image_centroid_y[platter_row] # 4044
  
  #programmed step time
  step_time_min <- 112.4394 
  step_time_sec = step_time_min * 60
  step_angle <- (4/516) * 365
  n_slots <- (360/step_angle) #* step_time_min)/60)/24
  
  platter_df <- data.frame(time_period = 1:n_slots, timestep_start = rep(NA, n_slots), timestep_end = rep(NA,n_slots),
                           angle_start = rep(NA, n_slots), angle_end = rep(NA, n_slots))
  
  
  platter_df <- platter_df %>% 
    mutate(timestep_start = time_deploy + step_time_sec * time_period - step_time_sec,
           timestep_end = #step_time_min - (as.numeric(difftime(time_deploy, time_start, units = "mins")) %% step_time_min) + #time left before turn on deploy
                          time_deploy + step_time_sec * time_period - step_time_min,
           angle_start = -step_angle * time_period - step_angle + (90 + step_angle*2), #the blue line starts on the left hand side of the open slit, so the first slot is just over 90 degrees
           angle_end = -step_angle * time_period + (90 + step_angle *2))
  # seg_obj_n = NA,
  # seg_obj_size = NA)
  
  
  #scanned file
  # file_scanned <- "pp_scan_sampler9_d220624_r.tif"
  # file_path <- "C:/Users/dsk273/Box/Cornell/mentoring/student projects/summer 2022/Kent pollen catcher/platter_scans/processed_scan/"
  #here("Cornell", "mentoring", "student projects", "summer 2022", "Kent pollen catcher", "platter_scans", "processed_scan", "pp_scan_sampler9_d220624_r.tif")
  #rotated_image_r <- raster::stack("C:/Users/dsk273/Desktop/pp_scan_sampler_14_d220510_c.tif")
  rotated_image <- terra::rast(here("texas", "pollen platter", "TX_platter_scans", "processed_scans", 
                                    paste0(scanned_file_name_focal, "_c.tif")
                                    #"pp_scan_sampler9_d220624_r.tif"
  )
  )
  #rotated_image <- terra::rast(paste0(file_path, file_scanned))
  crs(rotated_image) <- NA
  terra::ext(rotated_image) <- c(0, ncol(rotated_image), 0, nrow(rotated_image))
  plot(rotated_image)
  #plotRGB(rotated_image)
  
  
  
  ### divide image up into chunks and save each ##############################################################################
  #the width of the sampled area on the platter
  platter_dist_inner <- 1640
  platter_dist_outer <- 3600
  
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
  
  plot(rotated_image) #plot(1:1)
  #loop through each individual time chunk and save a file for it
  for(i in 1:n_slots){ #max 42 non-overlapping strips
    #create a spatial vector out of the points and visual check
    pp_aoi <- rbind(c(platter_df$outer_points_x_start[i], platter_df$outer_points_y_start[i]), 
                    c(platter_df$inner_points_x_start[i], platter_df$inner_points_y_start[i]), 
                    c(platter_df$inner_points_x_end[i], platter_df$inner_points_y_end[i]), 
                    c(platter_df$outer_points_x_end[i], platter_df$outer_points_y_end[i]))
    pp_aoi_v <- vect(pp_aoi, type = "polygons")
    # plot(rotated_image) #plot(1:1)
    plot(pp_aoi_v, col='yellow', add = TRUE)
    
    #extract image
    #test <- terra::extract(rotated_image, pp_aoi_v)
    chunk_x_uncropped <- terra::crop(rotated_image, pp_aoi_v)
    #plot(chunk_x_uncropped)
    chunk_x <- terra::mask(chunk_x_uncropped, pp_aoi_v, updatevalue = 1)
    #plot(chunk_x)
    
    file_save_name <- paste0(gsub(x = scanned_file_name_focal, pattern = ".tif", replacement = ""), "_chunk_", i, ".tif")
    terra::writeRaster(chunk_x, 
                       
                       here("texas", "pollen platter", "TX_platter_scans", "processed_scans", "platter_chunks", 
                            file_save_name),
                       overwrite = TRUE)
  } #end loop for 
}#end loop for a sample


### run classification macro in FIJI/Labkit ############################################################################################

#For some reason the segmentation command isn't picking up anything when I run it through the macro menu. After a lot of failed troubleshooting,
#I just used the probability map and ran it on the directory where all chunks were.

#processing the files now to extract the number of pixels where the Juas probability > .5 for each chunk
#test on a single chunk
# test <- raster::raster("C:/Users/dsk273/Desktop/classified_chunk_images/pp_scanner_sampler_1_d201230_chunk_11_pm.tif")
# hist(test[test[]<.9])
# 
# test[test[] > 0.5] <- NA
# raster::plot(test)
# length(test[test[]<0.5])

files_to_scan <- dir("C:/Users/dsk273/Desktop/classified_chunk_images/")
setwd("C:/Users/dsk273/Desktop/classified_chunk_images/")

chunk_pm_results <- data.frame(files_to_scan, pol_pix_n = rep(NA, length(files_to_scan)))

#running through each probability map. This could be faster with purrr, but a for loop is still adequate here
for(i in 1:nrow(chunk_pm_results)){
  pol_rast <- raster::raster(chunk_pm_results$files_to_scan[i])
  chunk_pm_results$pol_pix_n[i] <- length(pol_rast[pol_rast[]<0.5])
}

write_csv(chunk_pm_results, here("texas", "pollen_platter", "TX_platter_analysis", "Labkit_classifications", "TX_platters_Juas_pixels.csv"))


  
### expand deployment df ##############################################################################################
#deployment_sheet <- read_csv(here("texas", "pheno", "pollen_platter_analysis", "pollen_platters_fs20_21_220807.csv"))
deployment_sheet <- read_csv(here("texas",  "pollen_platter", "pollen_platters_fs20_21_220812.csv"))

#programmed step time
step_time_min <- 112.4394 
step_time_sec = step_time_min * 60
step_angle <- (4/516) * 365
n_slots <- round(360/step_angle , 0) #* step_time_min)/60)/24

rotation_time_sec <- (((step_time_min * n_slots)/60)/24) * 24 * 60 * 60

#file with results of pixel classification in each chunk (see above)
chunk_pm_results <- read_csv(here("texas", "pollen_platter", "TX_platter_analysis", "Labkit_classifications", "TX_platters_Juas_pixels.csv"))

deploy_join <- deployment_sheet %>% 
      dplyr::select(sampler, treatment, site, POINT_X, POINT_Y, long, lat, date_deploy_auto, date_retreive_auto, retreival_angle, scanned_file) 
pd <- expand_grid(deploy_join, data.frame(time_chunk = 1:n_slots)) %>% 
  mutate(files_to_scan = paste0(scanned_file, "_chunk_", time_chunk, "_pm.tif")) %>% 
  left_join(., chunk_pm_results) %>%  #add in the results from the pollen classification
  mutate(date_deploy_auto = mdy_hm(date_deploy_auto),
         date_retreive_auto = mdy_hm(date_retreive_auto),
    chunk_time_start = date_deploy_auto + step_time_sec * time_chunk - step_time_sec,
    chunk_time_end = date_deploy_auto + step_time_sec * time_chunk,
    chunk_hr_med = chunk_time_start + step_time_min/2,
    chunk_hr = hour(chunk_hr_med)) %>% 
  mutate(retreival_angle_c = case_when(retreival_angle < 0 ~ 360 + retreival_angle, TRUE ~ retreival_angle),
         deploy_duration = date_retreive_auto - date_deploy_auto,
         deploy_duration_num = as.numeric(deploy_duration),
         time_past_full_rotation = case_when(deploy_duration_num > as.numeric(rotation_time_sec) ~  (deploy_duration_num - as.numeric(rotation_time_sec)),
                                             TRUE ~ 0), #in seconds
         retreival_angle_predicted = case_when(time_past_full_rotation < 1 ~ as.numeric(deploy_duration/step_time_sec) * step_angle,
                                               time_past_full_rotation > 1 ~ as.numeric(deploy_duration/step_time_sec) * step_angle -360),
         retreival_angle_dif = retreival_angle_c - retreival_angle_predicted,
         
         angle_chunk_programmed = step_angle * time_chunk + step_angle/2,
  
         #for deployments duration that was under a single rotation time
         # which chunks shouldn't have been reached
         chunk_problem = case_when( 
                                time_past_full_rotation > (rotation_time_sec *2) ~ "platter left out long term",
                                time_past_full_rotation < 1 & 
                                angle_chunk_programmed > retreival_angle_predicted  ~ "chunk shouldn't have been reached", 
                                  time_past_full_rotation < 1 & 
                                angle_chunk_programmed >  retreival_angle_c  ~ "platter stopped before reaching chunk",
                                retreival_angle_c > retreival_angle_predicted + 15 ~ "platter moved too fast",
         
         #for deployment duration that exceeded the programmed single rotation
         #when a platter wasn't retrieved in time and potentially overwrote data
                                time_past_full_rotation > 1 & 
                                angle_chunk_programmed < (as.numeric(time_past_full_rotation/step_time_sec) * step_angle) &
                                retreival_angle_c > angle_chunk_programmed ~ "chunk overwritten",
                     
                                time_past_full_rotation > 1 & 
                               angle_chunk_programmed < (as.numeric(time_past_full_rotation/step_time_sec) * step_angle) &
                               retreival_angle_c < angle_chunk_programmed ~ "okay",
       
         #when a platter should have overwritten but stopped first
           
         TRUE ~ "okay"
         )
          # chunk_overwritten = case_when(angle_chunk_programmed < angles_overshot ~ "overwritten", 
          #                               TRUE ~ "okay")
        # all_types_bad_data = case_w
         )

#normalizing by platter max
max_per_platter <- pd %>%  group_by(scanned_file) %>% 
  summarize(p_max = max(pol_pix_n))

pd <- left_join(pd, max_per_platter)

pd <- pd %>% group_by(scanned_file) %>% 
  mutate( 
    pol_pix_ma24 = slider::slide_dbl(pol_pix_n, mean, 
                                 .before = round(12/(step_time_min/60), 0),
                                 .after =  round(12/(step_time_min/60), 0)),
    pol_pix_ma24_rel = pol_pix_n/pol_pix_ma,
    pol_pix_ma48 = slider::slide_dbl(pol_pix_n, mean, 
                                     .before = round(24/(step_time_min/60), 0),
                                     .after =  round(24/(step_time_min/60), 0)),
    pol_pix_ma48_rel = pol_pix_n/pol_pix_ma
    #roll_day = slider::slide_index_dbl(pol_pix_n/p_max, chunk_hr_med, mean, .before = days(.5))
  )



### connecting with other data: individual and site pheno curves ####################################################
deployment_sheet <- read_csv(here("texas",  "pollen_platter", "pollen_platters_fs20_21_220809.csv"))

#read in the site polygons 
site_poly <- st_read("C:/Users/danka/Box/texas/pheno/site_names_fs2020_21.shp")





# #extract the coordinates (lat, long, from the projected coordinates)
# #I just added it back to the original csv for convenience
deployment_sheet_coords_sf <- group_by(deployment_sheet, GlobalID) %>%
  dplyr::select(POINT_X, POINT_Y) %>%
  sf::st_as_sf(., coords = c("POINT_X", "POINT_Y"), crs = 32614)  #the coordinates are in 32614 (WGS UTM 14N)

deployment_sheet_coords_sf %>%  
  sf::st_transform(., crs = 4326) %>% #projecting back to 4326
  sf::st_coordinates(.) %>% as.data.frame() #extracting the coordinates
write.table(deployment_sheet_coords, "clipboard", sep="\t", row.names=FALSE)
# #deployment_sheet <- deployment_sheet %>% mutate(long = deployment_sheet_coords$X, lat = deployment_sheet_coords$Y)


deployment_sheet_site_name_by_coords <- st_join(deployment_sheet_coords_sf, site_poly) %>% 
  dplyr::select(GlobalID, site_name)
st_geometry(deployment_sheet_site_name_by_coords) <- NULL

deployment_sheet <- left_join(deployment_sheet, deployment_sheet_site_name_by_coords)  
#test <- deployment_sheet %>% dplyr::select(site, site_name)


### add in met data from daymet #######################################################################################
library(daymetr)
#start with the pixels of each NAB station
 site_sf_coords <- deployment_sheet %>% dplyr::select(site_name, lat, long) %>% 
                    mutate(lat = round(lat, 2),
                           long = round(long, 2)) %>% distinct()
#test <- download_daymet(lon = deployment_sheet$long[1], lat = deployment_sheet$lat[1], start =2015, end = 2017, simplify = TRUE)
setwd(here("texas",  "pheno", "met_data"))
write_csv(site_sf_coords, "platters_site_sf_coords.csv")
weather_at_platters <- download_daymet_batch(file_location = "platters_site_sf_coords.csv", start =2019, end = 2021, simplify = TRUE)
write_csv(weather_at_platters, "weather_at_platter_sites_220809.csv")
unique(weather_at_platters$measurement)

weather_at_platters <- read_csv(here("texas",  "pheno", "met_data", "weather_at_platter_sites_220809.csv") )%>% 
  mutate(date = as.Date(paste(year, yday, sep = "-"), "%Y-%j")) %>%
  mutate(measurement = gsub(pattern = ".", replacement = "", x = measurement, fixed = TRUE)) %>%
  dplyr::select(site_name = site, date, measurement, value) %>%
  pivot_wider(id_cols = c(site_name, date), names_from = measurement, values_from = value, names_prefix = "met_")
#head(weather_at_stations)




### read in temperature and humidity from data loggers and local stations ##################################################################################

datalogger_meta <- read_csv(here("texas",  "pheno", "met_data", "site_data_logger", "data_logger_metadata_220812_f.csv")) %>% 
  mutate(time_date_deploy = mdy_hm(time_date_deploy),
         time_date_retreive = mdy_hm(time_date_retreive))

#could make a function and extract ID from file name, but with just a few, I'm doing it manually for the moment
comal_data_logger_1 <- read_csv(here("texas",  "pheno", "met_data", "site_data_logger", "Comal_tempRH_210107_210212_inside_20987684.csv"),
                                skip = 2, col_names = c("ID", "date_time_", "temp_c", "rh","dewpoint"), col_select = 1:5) %>% 
  mutate(data_logger_ID = 20987684)

comal_data_logger_2 <- read_csv(here("texas",  "pheno", "met_data", "site_data_logger", "Comal_tempRH_210107_210212_outside_20987683.csv"),
                                skip = 2, col_names = c("ID", "date_time_", "temp_c", "rh","dewpoint"), col_select = 1:5) %>% 
  mutate(data_logger_ID = 20987683)

hays_data_logger_1 <- read_csv(here("texas",  "pheno", "met_data", "site_data_logger", "Hays_tempRH_210105_210209_inside_20987742.csv"),
                               skip = 2, col_names = c("ID", "date_time_", "temp_c", "rh","dewpoint"), col_select = 1:5) %>% 
  mutate(data_logger_ID = 20987742)

hays_data_logger_2 <- read_csv(here("texas",  "pheno", "met_data", "site_data_logger", "Hays_tempRH_210105_210209_outside_20987686.csv"),
                               skip = 2, col_names = c("ID", "date_time_", "temp_c", "rh","dewpoint"), col_select = 1:5) %>% 
  mutate(data_logger_ID = 20987686)

wade_data_logger_1 <- read_csv(here("texas",  "pheno", "met_data", "site_data_logger", "Wade_tempRH_210112_210210_inside_20987743.csv"),
                               skip = 2, col_names = c("ID", "date_time_", "temp_c", "rh","dewpoint"), col_select = 1:5) %>% 
  mutate(data_logger_ID = 20987743)

wade_data_logger_2 <- read_csv(here("texas",  "pheno", "met_data", "site_data_logger", "Wade_tempRH_210112_210210_outside_20987685.csv"),
                               skip = 2, col_names = c("ID", "date_time_", "temp_c", "rh","dewpoint"), col_select = 1:5) %>% 
  mutate(data_logger_ID = 20987685)

hobos <- bind_rows(comal_data_logger_1, comal_data_logger_2, hays_data_logger_1, hays_data_logger_2, wade_data_logger_1, wade_data_logger_2)

hobos <- left_join(hobos, datalogger_meta) %>% 
  mutate(time_date_meas = mdy_hms(date_time_),
         temp_f = temp_c * (9/5) + 32) %>% 
  filter(time_date_meas > time_date_deploy & time_date_meas < time_date_retreive)

# ggplot(hobos, aes(x= time_date_meas, y = temp_f, color = treatment)) + geom_line() + theme_bw() + facet_wrap(~site)

### wind logger
windlogger <- read_csv(here("texas",  "pheno", "met_data", "site_data_logger", "Wade_windlogger_jan2021.csv")) %>% 
  mutate(time_date_meas = mdy_hm(dt_mdy_hm),
         data_logger_ID = 999) %>% dplyr::select(-dt_mdy_hm)
windlogger <- left_join(windlogger, datalogger_meta)  %>% 
  filter(time_date_meas > time_date_deploy & time_date_meas < time_date_retreive)

# ggplot(windlogger, aes(x = time_date_meas, y = Speed)) + geom_line() + theme_bw()

### ambient weather stations (courtesy of Joe Wilson at Comal)
hays_ambient <- read_csv(here("texas",  "pheno", "met_data", "site_data_logger", "Hays_ambient-weather-20200901-20210211.csv"),
                         name_repair = "universal") %>%   mutate(data_logger_ID = 9991)
wade_ambient <- read_csv(here("texas",  "pheno", "met_data", "site_data_logger", "Wade_ambient-weather-20200901-20210211_30.8493_neg_98.0088.csv"),
                         name_repair = "universal") %>%   mutate(data_logger_ID = 9992)
comal_ambient <- read_csv(here("texas",  "pheno", "met_data", "site_data_logger", "Comal_ambient-weather-20200901-20210211.csv"),
                          name_repair = "universal") %>%   mutate(data_logger_ID = 9993)

ambient <- bind_rows(hays_ambient, wade_ambient, comal_ambient)
ambient2 <- ambient %>% 
  mutate(time_date_meas = substr(Date, 1, 16),
         time_date_meas = ymd_hm(time_date_meas)) %>% 
  rename(temp_f = Outdoor.Temperature,
         rh = Outdoor.Humidity) 
ambient3 <- ambient2 %>% 
  left_join(., datalogger_meta) %>% 
  filter(time_date_meas > time_date_deploy & time_date_meas < time_date_retreive)
# 
# ggplot(ambient3, aes(x = time_date_meas,  y = temp_f, color = treatment)) + geom_line() + theme_bw() + facet_wrap(~site) 
# 
# ggplot(ambient3, aes(x = time_date_meas,  y = temp_f)) + geom_point(color = "red", size = 0.4) + theme_bw() + facet_wrap(~site) +
#   geom_point(data = hobos2, aes(x = time_date_meas, y = temp_f), color = "blue", size = 0.2)



### add the met data to the pollen platter data #######################################################################

### Hobo data
# for each row, extract the chunk start time and the stop time
hobo_extract_fun <- function(df){
  f_sampler <- df$sampler
  f_start <- df$chunk_time_start
  f_end <- df$chunk_time_end
  f_site <- df$site
  f_tmt <- "outside" # temporary placeholder, saying all sites are "outside" 
  
  met <-
    hobos %>% 
    filter(site == f_site) %>% 
    filter(treatment == f_tmt) %>% 
    filter(time_date_meas > f_start & time_date_meas < f_end) %>% 
    summarize(temp_c_hobo = mean(temp_c, na.rm = TRUE),
              rh_hobo = mean(rh, na.rm = TRUE)) %>% 
    mutate(chunk_time_start = f_start,
           chunk_time_end = f_end,
           site = f_site,
           sampler = f_sampler) %>% 
    ungroup()
  return(met)
}

#test <- pd[8812:8912, ] %>% ungroup()
hobo_chunks <- pd %>% 
  split(1:nrow(.)) %>% #needs to be split to run on a df: https://stackoverflow.com/questions/55018739/run-purrrmap-dfr-on-dataframe-rows
  map_dfr(.x = ., .f = hobo_extract_fun)
beepr::beep(3)

pd <- left_join(pd, hobo_chunks)


### Ambient data
# for each dataframe row, extract the chunk start time and the stop time
ambient_extract_fun <- function(df){
  f_sampler <- df$sampler
  f_start <- df$chunk_time_start
  f_end <- df$chunk_time_end
  f_site <- df$site
  f_tmt <- "outside" # temporary placeholder, saying all sites are "outside" 
  
# average all env vars in within the start and stop time
  met <-
    ambient3 %>%  #names(ambient3)
    filter(site == f_site) %>% 
    filter(treatment == f_tmt) %>% 
    filter(time_date_meas > f_start & time_date_meas < f_end) %>% 
    summarize(temp_f_amb = mean(temp_f, na.rm = TRUE),
              rh_f_amb = mean(rh, na.rm = TRUE),
              feelslike_f_amb = mean(Feels.Like, na.rm = TRUE),
              Dew.Point_amb = mean(Dew.Point, na.rm = TRUE),
              windspeed_amb = mean(Wind.Speed, na.rm = TRUE),
              gust_amb = mean(Wind.Gust, na.rm = TRUE),
              wind_direction_amb = mean(Wind.Direction, na.rm = TRUE),
              rain_hourly_amb = mean(Hourly.Rain, na.rm = TRUE),
              pressure_rel_amb = mean(Relative.Pressure, na.rm = TRUE),
              pressure_abs_amb = mean(Absolute.Pressure, na.rm = TRUE),
              sol_rad_amb = mean(Solar.Radiation, na.rm = TRUE)) %>% 
    mutate(chunk_time_start = f_start,
           chunk_time_end = f_end,
           site = f_site,
           sampler = f_sampler) %>% 
    ungroup()
  return(met)
}

ambient_chunks <- pd %>%  #pd[8812:8912, ] 
  split(1:nrow(.)) %>% #needs to be split to run on a df: https://stackoverflow.com/questions/55018739/run-purrrmap-dfr-on-dataframe-rows
  map_dfr(.x = ., .f = ambient_extract_fun)
beepr::beep(3)

pd <- left_join(pd, ambient_chunks)


# write_csv(pd, here("texas",  "pollen_platter", "pollen_platters_fs20_21_processed_weather_220812.csv"))

### preliminary data vis ##############################################################################################
#pd <- read_csv(here("texas",  "pollen_platter", "pollen_platters_fs20_21_processed_weather_220812.csv"))
#pd2 <- filter(pd, chunk_hr_med > mdy_hm("1/1/2021 9:00")) 
get.es <- function(temp){
  es <- 6.11 * exp((2.5e6 / 461) * (1 / 273 - 1 / (273 + temp)))
  return(es)
}

get.vpd <- function(rh, temp){
  ## calculate saturation vapor pressure
  es <- get.es(temp)
  ## calculate vapor pressure deficit
  vpd <- ((100 - rh) / 100) * es
  return(vpd)
}

pd2 <- filter(pd, chunk_problem == "okay") %>% 
  mutate(temp_c_amb = (temp_f_amb - 32)*(5/9),
         es_amb = 0.6108 * exp(17.27 * temp_c_amb / (temp_c_amb + 237.3)),
         vpd_amb = rh_f_amb/(100 * es_amb) - es_amb) %>% 
  #https://physics.stackexchange.com/questions/4343/how-can-i-calculate-vapor-pressure-deficit-from-temperature-and-relative-humidit
  mutate(vpd2_amb = get.vpd(rh_f_amb, temp_c_amb))
  

hist(pd2$vpd2_amb)
hist(pd2$vpd_amb)




  
ggplot(pd2, aes(x = chunk_hr_med, y = pol_pix_n, col = treatment)) + geom_line() + facet_wrap(~scanned_file, scales = "free") + theme_bw()



# hour of observations
pd2 %>% 
  group_by(chunk_hr) %>% 
  summarize(p_per_hour_ma_rel_mean = mean(pol_pix_ma_rel),
            p_per_hour_ma_rel_sd = sd(pol_pix_ma_rel)) %>% 
ggplot(aes(x =  chunk_hr, y = p_per_hour_ma_rel_mean)) + geom_point() + theme_bw() + geom_smooth()


pd2 %>% 
  ggplot(aes(x = as.factor(chunk_hr), y = pol_pix_ma_rel)) + geom_boxplot() + facet_wrap(~treatment, ncol = 1)
  

### looking at a single platter
pd2 %>% filter(scanned_file == "pp_scanner_sampler_14_d210114") %>% 
  ggplot(aes(x = chunk_hr_med, y = pol_pix_n, group = scanned_file)) + geom_line(size = 1)  + 
  theme_bw() + #facet_wrap(~scanned_file, scales = "free")+
  scale_color_viridis_c() + xlab("date") + ylab("pollen abundance (pixels)")

### how effective was the porch screen at reducing pollen? (open air vs porch control)
pd2 %>% 
  mutate(date_deploy = as_date(date_deploy_auto),
         date_deploy_loc = paste(date_deploy, POINT_X, POINT_Y)) %>% 
  filter(date_deploy_loc == "2021-01-12 589852.0617 3410760.884" | 
         date_deploy_loc == "2020-12-29 589852.0704 3410760.872" |
         date_deploy_loc == "2021-01-22 580013.9055 3299866.384") %>% 
  filter(treatment == "porch control" | treatment == "open air") %>% 
  group_by(treatment, date_deploy_loc) %>% 
  summarize(n = n(),
            p_mean = mean(pol_pix_n, na.rm = TRUE),
            p_se = sd(pol_pix_n, na.rm = TRUE)/sqrt(n)) %>% 
  ggplot(aes(x = date_deploy_loc, y = p_mean, fill = treatment, ymin = p_mean - p_se, ymax = p_mean + p_se)) + geom_col(position = position_dodge())  + 
  geom_errorbar(position = position_dodge(1), width = 0.2)+
  theme_bw(base_size = 16)  + xlab("deployments") + ylab("pollen concentration (number of pixels)") + scale_x_discrete(labels = c("", "", ""))+
  scale_fill_discrete( labels = c("without screen", "with screen"))


### A long time series
pd2 %>% 
  mutate(date_deploy = as_date(date_deploy_auto),
         date_site = paste(date_deploy, site),
         location = paste(POINT_X, POINT_Y)) %>% 
  filter(treatment != "contamination control" & treatment != "porch control") %>% 
  filter(location == "589853.9347 3410801.385") %>% 
  filter(site == "Wade") %>%  #| site == "Comal" | site == "Wade") %>%  #pol_pix_n/p_max
  ggplot(aes(x = chunk_hr_med, y = pol_pix_n, group = scanned_file)) + geom_line(size = 1)  + 
  scale_x_datetime(date_breaks = "1 day", date_labels =  "%d %b")  + #scale_color_viridis_c() +
  theme_bw(base_size = 16) + xlab("date") + ylab("pollen concentration (number of pixels)") #+ facet_wrap(~location, scales = "free") 

pd2 %>% 
  mutate(date_deploy = as_date(date_deploy_auto),
         date_site = paste(date_deploy, site),
         location = paste(POINT_X, POINT_Y)) %>% 
  filter(treatment != "contamination control" & treatment != "porch control") %>% 
  filter(location == "589853.9347 3410801.385") %>% 
  ungroup() %>% 
  arrange(chunk_hr_med) %>% 
  mutate(pol_pix_cum = cumsum(pol_pix_n)) %>% 
  filter(site == "Wade") %>%  #| site == "Comal" | site == "Wade") %>%  #pol_pix_n/p_max
  ggplot(aes(x = chunk_hr_med, y = pol_pix_cum)) + geom_line(size = 1)  + 
  scale_x_datetime(date_breaks = "4 day", date_labels =  "%d %b")  + scale_color_viridis_c() +
  theme_bw() + xlab("date") + ylab("cumulative pollen (number of pixels)") #+ facet_wrap(~location, scales = "free") 


#env data at a long time series  
pd2 %>% 
  mutate(date_deploy = as_date(date_deploy_auto),
         date_site = paste(date_deploy, site),
         location = paste(POINT_X, POINT_Y)) %>% 
  filter(treatment != "contamination control" & treatment != "porch control") %>% 
  filter(location == "589853.9347 3410801.385") %>% 
  filter(site == "Wade") %>%  #| site == "Comal" | site == "Wade") %>%  #pol_pix_n/p_max
  ggplot(aes(x = chunk_hr_med, y = vpd2_amb, group = scanned_file, col = pol_pix_ma48_rel)) + geom_line(size = 1)  + 
  scale_x_datetime(date_breaks = "1 day", date_labels =  "%d %b")  + scale_color_viridis_c(name = "pollen (pixels)") +
  theme_bw(base_size = 16) + xlab("date") + ylab("relative humidity (%)")#ylab("pressure (inHg)") #+ facet_wrap(~location, scales = "free") 




pd2 %>% 
  mutate(date_deploy = as_date(date_deploy_auto),
         date_site = paste(date_deploy, site),
         location = paste(POINT_X, POINT_Y)) %>% 
  filter(location == "589853.9347 3410801.385") %>% 
  ggplot(aes(x = pressure_rel_amb, y = pol_pix_n)) + geom_point(size = 1)  + 
  #scale_x_datetime(date_breaks = "1 day", date_labels =  "%d %b")  + scale_color_viridis_c(name = "pollen (pixels)") +
  theme_bw() + xlab("pressure (inHg)") + ylab("pollen release (pixels") + geom_smooth(method ="lm")



###### 
pd2 %>% 
  mutate(date_deploy = as_date(date_deploy_auto),
         date_site = paste(date_deploy, site),
         location = paste(POINT_X, POINT_Y)) %>% 
  filter(treatment != "contamination control" & treatment != "porch control") %>% 
 # filter(location == "589853.9347 3410801.385") %>% 
  filter(site == "Hays" | site == "Comal" | site == "Wade") %>%  #pol_pix_n/p_max
  ggplot(aes(x = chunk_hr_med, y = pol_pix_n, col = vpd2_amb,group = scanned_file)) + geom_line(size = 1)  + 
  #scale_x_datetime(date_breaks = "1 day", date_labels =  "%d %b")  + 
  scale_color_viridis_c(name = "vpd (mb)") +
  theme_bw() + xlab("date") + ylab("pollen concentration (number of pixels)") + facet_wrap(~site, scales = "free_y", ncol = 1) 


#how well correlated are environmental variables between sites
pd2 %>% 
  ggplot(aes(x = date))


### looking at time series and raw data
names(pd2)
pd2 %>% 
  filter(treatment != "contamination control") %>% 
  #filter(treatment != "porch control") %>% 
  #filter(treatment == "open air") %>% 
  filter(retreival_angle_dif > - 6 & retreival_angle_dif < 6) %>% 
  filter(site == "Hays" | site == "Comal" | site == "Wade") %>%  #pol_pix_n/p_max
  ggplot(aes(x = chunk_hr_med, y = pol_pix_ma_rel, col = vpd_amb, group = scanned_file)) + geom_line(size = 1)  + 
  theme_bw() + facet_grid(treatment~site, scales = "free")+
  scale_color_viridis_c()


### looking at platters that should have good temporal accuracy (angle was within 4 hrs)
hist(pd2$retreival_angle_dif, breaks = 200)
names(pd2)
pd2 %>% 
  filter(treatment != "contamination control") %>% 
  #filter(treatment != "porch control") %>% 
  #filter(treatment == "open air") %>% 
  filter(retreival_angle_dif > 0 & retreival_angle_dif < 15) %>% 
  filter(site == "Hays" | site == "Comal" | site == "Wade") %>%  #pol_pix_n/p_max
  ggplot(aes(x = chunk_hr_med, y = pol_pix_n, col = temp_f_amb, group = scanned_file)) + geom_line(size = 1)  + 
  theme_bw() + facet_wrap(~scanned_file, scales = "free")+ scale_color_viridis_c()



### looking at a few that seem to have especially clear signals
some_clear_platters <- c("pp_scanner_sampler_22_d210112", "pp_scanner_sampler_12_d210114",
                         "pp_scanner_sampler_14_d210114", "pp_scanner_sampler_14_d210105",
                         "pp_scanner_sampler_11_d210129", "pp_scanner_sampler_11_d210121",
                         "pp_scanner_sampler_16_d210112", "pp_scanner_sampler_19_d201230")
pd2 %>% filter(scanned_file %in% some_clear_platters) %>% 
  filter(retreival_angle_dif > - 6 & retreival_angle_dif < 6) %>% 
  ggplot(aes(x = chunk_hr_med, y = pol_pix_n, col = temp_f_amb, group = scanned_file)) + geom_line(size = 1)  + 
  theme_bw() + facet_wrap(~scanned_file, scales = "free")+
  scale_color_viridis_c()
  
  
#direct comparison of chunk vs the 24 hour period
pd2 %>% 
  filter(treatment != "contamination control") %>%
  #filter(retreival_angle_dif > - 10 & retreival_angle_dif < 10) %>% 
  #filter(treatment != "porch control") %>% 
  #filter(treatment == "open air") %>% 
  filter(site == "Hays" | site == "Comal" | site == "Wade") %>%  #pol_pix_n/p_max
ggplot(aes(x = vpd_amb, y = pol_pix_ma48_rel, col = scanned_file)) + geom_point() + theme_bw() + facet_grid(treatment~site, scales = "free")+
  scale_color_viridis_d() + geom_smooth(method = "lm")


ggplot(pd2, aes(x = rh_f_amb, y = vpd_amb, col = temp_c_amb)) + geom_point()



### old stuff ##############################################################################################

pol_dep %>% 
  filter(processed_file == "pp_scan_samp_11_d220506_r") %>% 
  ggplot(aes(x = time_window_start, y = Count)) + geom_point(aes(color = time_window_hour)) +  theme_bw() + facet_wrap(~processed_file, scales = "free_y") +
  scale_color_viridis_c( ) + 
  geom_line(aes(x = time_window_start, y = zoo::rollmean(Count, 5, align = "center", fill = NA))) 


pol_dep %>% 
  filter(processed_file == "pp_scan_samp_23_d220517_r") %>% 
  # filter(deployment_time > mdy_hm("6/1/22 9:00")) %>% 
  # filter(deployment_time < mdy_hm("7/17/22 9:00")) %>% 
  filter(species == "Quru") %>% 
  #filter(time_period < 12) %>% 
  ggplot(aes(x = time_into_deploy_falsedate, y = Count)) + geom_point(aes(color = time_window_hour)) +  theme_bw() + facet_wrap(~processed_file, scales = "free_y") +
  scale_color_viridis_c( ) + 
  geom_line(aes(x = time_into_deploy_falsedate, y = zoo::rollmean(Count, 5, align = "center", fill = NA))) 
#coord_cartesian(ylim= c(0, 1000)) 
# geom_vline(xintercept = mdy_hm("6/1/22 18:00"), color = "blue") +
# geom_vline(xintercept = mdy_hm("6/2/22 6:00"), color = "yellow") +
# geom_vline(xintercept = mdy_hm("6/2/22 18:00"), color = "blue") +
# geom_vline(xintercept = mdy_hm("6/3/22 6:00"), color = "yellow") +
# geom_vline(xintercept = mdy_hm("6/3/22 18:00"), color = "blue") 


#combine env vars with pollen platters
pol_dep$time_window_med[100]
temp_twigs$date_time

x <- pol_dep$time_window_med[100]
env_index <- which(abs(temp_twigs$date_time-x) == min(abs(temp_twigs$date_time - x)))
temp_twigs$date_time[env_index]

# for(i in 1:10){
#   
# }##
