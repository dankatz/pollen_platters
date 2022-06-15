### pollen platter analysis ###

#
library(lubridate)
library(dplyr)

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
step_time_sec = step_time * 60
step_angle <- (12/516) * 365

platter_df <- data.frame(period = 1:48, timestep_start = rep(NA, 48), timestep_end = rep(NA,48),
                         angle_start = rep(NA, 48), angle_end = rep(NA, 48))


platter_df %>% 
  mutate(timestep_start = time_deploy + step_time_sec * period - step_time_sec,
         timestep_end = step_time_min - (as.numeric(difftime(time_deploy, time_start, units = "mins")) %% step_time_min) + #time left before turn on deploy
                        time_deploy + 
                        step_time_sec * period - step_time_min,
         angle_start = step_angle * period - step_angle,
         angle_end = step_angle * period)


#scanned file
file_scanned <- "pp_scan_sampler_14_d220510"

#rotate scanned file 

#divide image up into chunks

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
