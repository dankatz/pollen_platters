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

setwd("C:/Users/danka/Box")
here::i_am("katz_photo.jpg")


# #
# test <- read.delim("clipboard")
# test4 <- test %>% as_tibble(.) %>% 
# mutate(
# site = gsub("(.*),.*", "\\1", sampler.text),
# test2 = gsub(".*,","",sampler.text),
# test3 = mdy(test2))
# test4
# #writeClipboard(test4$test3)
# write.table(test4, "clipboard", sep="\t", row.names=FALSE, col.names=FALSE)


### connect to deployment spreadsheet #######################################################################################
deployment_sheet <- read_csv(here("texas", "pheno", "pollen_platter_analysis", "pollen_platters_fs20_21_220808.csv"))

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

### read in temperature and humidity from logger ##################################################################################
# temp_twigs <- read_csv(here("Cornell", "Local projects", "twig pheno 2022", "intertwined_with_twigs_May2022_20987743_c.csv")) %>% 
#   mutate(date_time = mdy_hm(date_time))
# temp_outside <- read_csv(here("Cornell", "Local projects", "twig pheno 2022", "outside_high_tunnel_May2022_20987685_c.csv")) %>% 
#   mutate(date_time = mdy_hm(date_time))
# temp_inside <- read_csv(here("Cornell", "Local projects", "twig pheno 2022", "inside_high_tunnel_May22_20987683_c.csv")) %>% 
#   mutate(date_time = mdy_hm(date_time))

# 
# ggplot(temp_outside, aes(x= date_time, y = Temp_c)) + geom_line(col = "blue") + theme_bw() +
#   geom_line(data = temp_inside, aes(x = date_time, y = Temp_c), col = "red") +
#   geom_line(data = temp_twigs, aes(x = date_time, y = Temp_c), col = "green")
# 
# ggplot(temp_outside, aes(x= date_time, y = RH_p)) + geom_line(col = "blue") + theme_bw() +
#   geom_line(data = temp_inside, aes(x = date_time, y = RH_p), col = "red") +
#   geom_line(data = temp_twigs, aes(x = date_time, y = RH_p), col = "green")


### expand deployment df ##############################################################################################
deployment_sheet <- read_csv(here("texas", "pheno", "pollen_platter_analysis", "pollen_platters_fs20_21_220807.csv"))

#programmed step time
step_time_min <- 112.4394 
step_time_sec = step_time_min * 60
step_angle <- (4/516) * 365
n_slots <- round(360/step_angle , 0) #* step_time_min)/60)/24

rotation_time_sec <- (((step_time_min * n_slots)/60)/24) * 24 * 60 * 60

#file with results of pixel classification in each chunk (see above)
chunk_pm_results <- read_csv(here("texas", "pollen_platter", "TX_platter_analysis", "Labkit_classifications", "TX_platters_Juas_pixels.csv"))

deploy_join <- deployment_sheet %>% 
      dplyr::select(sampler, site, POINT_X, POINT_Y, date_deploy_auto, date_retreive_auto, retreival_angle, scanned_file) 
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
    pol_pix_ma = slider::slide_dbl(pol_pix_n, mean, 
                                 .before = round(12/(step_time_min/60), 0),
                                 .after =  round(12/(step_time_min/60), 0)),
    pol_pix_ma_rel = pol_pix_n/pol_pix_ma
    #roll_day = slider::slide_index_dbl(pol_pix_n/p_max, chunk_hr_med, mean, .before = days(.5))
  )

### preliminary data vis ##############################################################################################
#pd2 <- filter(pd, chunk_hr_med > mdy_hm("1/1/2021 9:00")) 
pd2 <- filter(pd, chunk_problem == "okay") 
  
ggplot(pd2, aes(x = chunk_hr_med, y = pol_pix_n/p_max, col = chunk_problem)) + geom_line() + facet_wrap(~scanned_file, scales = "free") + theme_bw()


# hour of observations
pd2 %>% 
  group_by(chunk_hr) %>% 
  summarize(p_per_hour_ma_rel_mean = mean(pol_pix_ma_rel),
            p_per_hour_ma_rel_sd = sd(pol_pix_ma_rel)) %>% 
ggplot(aes(x =  chunk_hr, y = p_per_hour_ma_rel_mean)) + geom_point() + theme_bw() + geom_smooth()


pd2 %>% 
  ggplot(aes(x = as.factor(chunk_hr), y = pol_pix_ma_rel)) + geom_boxplot()
  
  
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
