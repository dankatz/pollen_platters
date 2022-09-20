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
library(googlesheets4)
library(stringr)
library(tidyr)

setwd("C:/Users/dsk273/Box")
here::i_am("katz_photo.jpg")




### connect to google drive spreadsheet #######################################################################################
deployment_sheet <- read_sheet("https://docs.google.com/spreadsheets/d/1h8XE4uVwhZ4Aez7e9cUUSics7iCskL6rrKBSnsEokdM/edit?usp=sharing") %>% 
  filter(!is.na(scanned_file_name))

platters_to_process_list_rows <- c(51:54)
for(j in 1:length(platters_to_process_list_rows)){

platter_row <- platters_to_process_list_rows[j] #platter_row <- platters_to_process_list_rows[1]
scanned_file_name_focal <- deployment_sheet$scanned_file_name[platter_row]

#sampler number
sampler_n <- deployment_sheet$sampler[platter_row]  #"9"

#time start
time_start <- ymd_hms(deployment_sheet$sampler_start_date_time[platter_row])

#time deployed
time_deploy <- ymd_hms(deployment_sheet$deployment_time[platter_row]) #mdy_hms("6/24/2022 14:59:00")

#time retreived
time_retreived <- ymd_hms(deployment_sheet$retrieval_time[platter_row]) #mdy_hms("6/26/2022 13:57:00")

#angle at retreival 
#angle_retreival <- deployment_sheet$retreival_angle[platter_row] + 90 #the measured angle is from the blue line to the left side of the opening slit. The angles used in this script start at 90 #8 + 360

#platter centroid (manually measured)
platter_centroid_x <- deployment_sheet$image_centroid_x[platter_row] # 4047 #
platter_centroid_y <- deployment_sheet$image_centroid_y[platter_row] # 4044

#programmed step time
step_time_min <- 65.717 
step_time_sec = step_time_min * 60
step_angle <- (12/516) * 365

platter_df <- data.frame(time_period = 1:44, timestep_start = rep(NA, 44), timestep_end = rep(NA,44),
                         angle_start = rep(NA, 44), angle_end = rep(NA, 44))


platter_df <- platter_df %>% 
  mutate(timestep_start = time_deploy + step_time_sec * time_period - step_time_sec,
         timestep_end = step_time_min - (as.numeric(difftime(time_deploy, time_start, units = "mins")) %% step_time_min) + #time left before turn on deploy
                        time_deploy + 
                        step_time_sec * time_period - step_time_min,
         angle_start = -step_angle * time_period - step_angle + (90 + step_angle*2), #the blue line starts on the left hand side of the open slit, so the first slot is just over 90 degrees
         angle_end = -step_angle * time_period + (90 + step_angle *2))
         # seg_obj_n = NA,
         # seg_obj_size = NA)


#scanned file
# file_scanned <- "pp_scan_sampler9_d220624_r.tif"
# file_path <- "C:/Users/dsk273/Box/Cornell/mentoring/student projects/summer 2022/Kent pollen catcher/platter_scans/processed_scan/"
#here("Cornell", "mentoring", "student projects", "summer 2022", "Kent pollen catcher", "platter_scans", "processed_scan", "pp_scan_sampler9_d220624_r.tif")
#rotated_image_r <- raster::stack("C:/Users/dsk273/Desktop/pp_scan_sampler_14_d220510_c.tif")
rotated_image <- terra::rast(here("Cornell", "mentoring", "student projects", "summer 2022", "Kent pollen catcher", "platter_scans", "processed_scan", 
                                  paste0(scanned_file_name_focal, "_r.tif")
                                  #"pp_scan_sampler9_d220624_r.tif"
                                  ))
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

#loop through each individual time chunk and save a file for it
for(i in 1:42){ #max 42 non-overlapping strips
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

  file_save_name <- paste0(gsub(x = scanned_file_name_focal, pattern = ".tif", replacement = ""), "_chunk_", i, ".tif")
  terra::writeRaster(chunk_x, 
                     
                     here("Cornell", "mentoring", "student projects", "summer 2022", "Kent pollen catcher", "platter_scans", "platter_chunks", 
                          file_save_name),
                     overwrite = TRUE)
} #end loop for 
}#end loop for a sample


### run classification macro in FIJI/Labkit ############################################################################################

#run the macro in FIJI:
#C:\Users\danka\Box\Cornell\mentoring\student projects\summer 2022\Kent pollen catcher\Labkit_classifications
#Macro_labkit_plantain2.ijm.ijm

#hacky workaround to save the relevant file names within a .txt file so I can open that within a macro in FIJI
platters_to_process_list_rows
file_name_list_txt <- deployment_sheet$scanned_file_name[min(platters_to_process_list_rows):max(platters_to_process_list_rows)]
write_delim( x = as.data.frame(file_name_list_txt), 
             file = here("Cornell", "mentoring", "student projects", "summer 2022", "Kent pollen catcher", "Labkit_classifications", 
              "hacky_file_list_workaround.txt"),
             col_names = FALSE)



#It's easiest to just open up image j manually and run the macro for all chunks
# system2('C:/Users/danka/Documents/Fiji.app/ImageJ-win64.exe', 
#         'C:/Users/danka/Box/Cornell/mentoring/student projects/summer 2022/Kent pollen catcher/Labkit_classifications/Macro_labkit_plantain4.ijm') #, "pp_scan_samp_22_d220629"

#currently trying it as a batch process; it seems labkit updated so it doesn't play nicely with the macro anymore

test <- raster("C:/Users/dsk273/Box/Cornell/mentoring/student projects/summer 2022/Kent pollen catcher/Labkit_classifications/classified_chunks_ragweed/pp_scan_samp_15_d220916_chunk_10_pm.tif")
plot(test)



### update as of sept 20, 2022, taken from the TX_platters.R script
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

files_to_scan <- dir("C:/Users/dsk273/Box/Cornell/mentoring/student projects/summer 2022/Kent pollen catcher/Labkit_classifications/classified_chunks_ragweed")
setwd("C:/Users/dsk273/Box/Cornell/mentoring/student projects/summer 2022/Kent pollen catcher/Labkit_classifications/classified_chunks_ragweed")

chunk_pm_results <- data.frame(files_to_scan, pol_pix_n = rep(NA, length(files_to_scan)))

#running through each probability map. This could be faster with purrr, but a for loop is still adequate here
for(i in 1:nrow(chunk_pm_results)){
  pol_rast <- raster::raster(chunk_pm_results$files_to_scan[i])
  chunk_pm_results$pol_pix_n[i] <- length(pol_rast[pol_rast[]<0.5])
}

write_csv(chunk_pm_results, "C:/Users/dsk273/Box/Cornell/mentoring/student projects/summer 2022/Kent pollen catcher/Labkit_classifications/classification_chunk_results_ragweed/ragweed_results220920_b.csv")


###NEED TO REWRITE FOR RAGWEED

### read in results from FIJI/Labkit and prepare data ##############################################################################################
#programmed step time
step_time_min <- 65.717 
step_time_sec = step_time_min * 60
step_angle <- (12/516) * 365
n_slots <- round(360/step_angle , 0) #* step_time_min)/60)/24
rotation_time_sec <- (((step_time_min * n_slots)/60)/24) * 24 * 60 * 60 #rotation_time_sec/(60*60)

result_csvs <- dir(here("Cornell", "mentoring", "student projects", "summer 2022", "Kent pollen catcher", "Labkit_classifications", "classification_chunk_results_ragweed"), full.names = TRUE) %>%  
                map_dfr(.x = ., .f = read_csv) 
  
result_csvs <- result_csvs %>% 
              rename_with(make.names) %>% 
              mutate(time_period = gsub(pattern = ".*_", replacement = "", x = Slice),
                     time_period = as.numeric(gsub(pattern = ".tif", replacement = "", x = time_period)),
                     scanned_file_name = substring(Slice, 1, 23),
                     scanned_file_name  = sub("_$", "", scanned_file_name)) %>% #for inconsistent file name length (ie 04 vs 4)
              dplyr::select(scanned_file_name, time_period, Count, Average.Size) %>% 
              arrange(time_period)

#unique(result_csvs$scanned_file_name)
# result_csvs %>% 
#   ggplot(aes(x = time_period, y = Count)) + geom_point() + geom_line() + theme_bw()

pol_dep_raw <- deployment_sheet %>% dplyr::select(scanned_file_name, sampler, sampler_start_date_time, species, deployment_time, retrieval_time, retreival_angle) %>% #time_period, timestep_start, timestep_end) %>% 
            left_join(., result_csvs)


pd <- pol_dep_raw %>% #expand_grid(pol_dep_raw, data.frame(time_chunk = 1:n_slots)) %>% 
  mutate(time_chunk = time_period, 
         chunk_id = paste0(scanned_file_name, "_chunk_", time_period)) %>% 
  #left_join(., chunk_pm_results) %>%  #add in the results from the pollen classification
  mutate(date_deploy_auto = deployment_time,
         date_retreive_auto = retrieval_time,
         chunk_time_start = date_deploy_auto + step_time_sec * time_chunk - step_time_sec,
         chunk_time_end = date_deploy_auto + step_time_sec * time_chunk,
         chunk_hr_med = chunk_time_start + step_time_min/2,
         chunk_hr = hour(chunk_hr_med)) %>% 
  mutate(retreival_angle = as.numeric(retreival_angle),
         retreival_angle_c = case_when(
                                       retreival_angle < 0 ~ 360 + retreival_angle, 
                                       retreival_angle > 360 ~ retreival_angle - 360,
                                       retrieval_time == ymd_hms("2022-07-29 18:00:00") ~ retreival_angle - #the one where Blue didn't take photos...
                                         ((interval(ymd_hms("2022-07-29 18:00:00"),ymd_hms("2022-08-01 10:30:00") )/ seconds(step_time_sec)) * step_angle) + 360*2,
                                       TRUE ~ retreival_angle),
         deploy_duration = date_retreive_auto - date_deploy_auto,
         deploy_duration_num = as.numeric(deploy_duration) * 60 * 60 * 24, #switching to seconds
         time_past_full_rotation = case_when(deploy_duration_num > as.numeric(rotation_time_sec) ~  (deploy_duration_num - as.numeric(rotation_time_sec)),
                                             TRUE ~ 0), #in seconds
         retreival_angle_predicted = case_when(time_past_full_rotation < 1 ~ as.numeric(deploy_duration_num/step_time_sec) * step_angle,
                                               time_past_full_rotation > 1 ~ as.numeric(deploy_duration_num/step_time_sec) * step_angle -360),
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
max_per_platter <- pd %>%  group_by(scanned_file_name) %>% 
  summarize(p_max = max(Count))

pd <- left_join(pd, max_per_platter)

pd <- pd %>% group_by(scanned_file_name) %>% 
  mutate( 
    pol_pix_n = Count,
    pol_pix_ma24 = slider::slide_dbl(pol_pix_n, mean, 
                                     .before = round(12/(step_time_min/60), 0),
                                     .after =  round(12/(step_time_min/60), 0)),
    pol_pix_ma24_rel = pol_pix_n/pol_pix_ma24,
    pol_pix_ma48 = slider::slide_dbl(pol_pix_n, mean, 
                                     .before = round(24/(step_time_min/60), 0),
                                     .after =  round(24/(step_time_min/60), 0)),
    pol_pix_ma48_rel = pol_pix_n/pol_pix_ma48
    #roll_day = slider::slide_index_dbl(pol_pix_n/p_max, chunk_hr_med, mean, .before = days(.5))
  )


pd <- pd %>% 
  filter(species == "plantain" | species == "control") %>% 
  mutate(deploy_time_falsedate = mdy_hm(paste0("6-1-2022 ", hour(deployment_time), ":", minute(deployment_time))),
         time_window_start = deployment_time + lubridate::dminutes(step_time_min) * (time_period - 1),
         time_window_start_hrm = paste(hour(time_window_start), minute(time_window_start), sep = ":"),
         time_window_end = deployment_time + lubridate::dminutes(step_time_min) * (time_period),
         time_window_med = difftime(time_window_end, time_window_start)/2 + time_window_start,
         #time_into_deploy = difftime(time_window_med, time_window_start),
         time_into_deploy = difftime(time_window_med, deployment_time), #+ time_window_med,
         time_into_deploy_falsedate = time_into_deploy + deploy_time_falsedate,
         time_window_hour = hour(time_window_med),
         time_window_min = minute(time_window_med),
         time_window_hm = mdy_hm(paste0("6-1-2022 ", time_window_hour, ":", time_window_min)),
         day_deploy = case_when(time_period < 24 ~ "day 1",
                                time_period >= 24 ~ "day 2"))
#compare observed angle to calculated angle 

pd %>% 
  filter(species == "plantain") %>% 
  ggplot(aes(x = chunk_hr_med , y = Count, group = scanned_file_name, col = chunk_problem)) + geom_line() 




pd %>% 
  filter(deployment_time > mdy_hm("7/20/22 9:00")) %>% 
  filter(species == "plantain") %>% 
  #filter(time_period > 0) %>% 
  ggplot(aes(x = time_window_med, y = Count)) + geom_point() +  theme_bw() + facet_wrap(~scanned_file_name, scales = "free") +
  #scale_color_viridis_d() + 
  geom_line(aes(x = time_window_med, y = zoo::rollmean(Count, 5, align = "center", fill = NA)))




pd %>% 
  filter(deployment_time > mdy_hm("7/17/22 9:00")) %>% 
  filter(deployment_time < mdy_hm("8/17/22 9:00")) %>% 
  filter(species == "plantain") %>% 
  #filter(chunk_problem == "okay") %>% 
  #filter(time_period < 12) %>% 
  ggplot(aes(x = time_into_deploy_falsedate, y = Count)) + geom_point() +  theme_bw() + facet_wrap(~scanned_file_name, scales = "free_y") +
  #scale_color_viridis_c( ) + 
  geom_line(aes(x = time_into_deploy_falsedate, y = zoo::rollmean(Count, 5, align = "center", fill = NA))) +
  #coord_cartesian(ylim= c(0, 1000)) +
  geom_vline(xintercept = mdy_hm("6/1/22 18:00"), color = "blue") +
  geom_vline(xintercept = mdy_hm("6/2/22 6:00"), color = "yellow") +
  geom_vline(xintercept = mdy_hm("6/2/22 18:00"), color = "blue") +
  geom_vline(xintercept = mdy_hm("6/3/22 6:00"), color = "yellow") +
  geom_vline(xintercept = mdy_hm("6/3/22 18:00"), color = "blue") 

