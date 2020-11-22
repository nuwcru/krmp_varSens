library(tidyverse)
library(lubridate)

d <- read_csv("data/Clean IVI years 2013-2019.csv") 

d <- d %>% select(-"cam model":-"trigger", -"julian", -"julian2", -"jdate", -"jdate2", -"start_day", -"start_month", -"year", -"end_day", -"end_month", -"sex")
names(d)[7] <- "visit_type"

d <- d %>%
  mutate(date = dmy(date),
         start = as.character(start),
         filled_end = as.character(filled_end),
         date_start = ymd_hms(paste(date, start))) %>%
  mutate(date = as.character(date)) %>%
  filter(visit_type == "provisioning")



# split dataframe into lists based on the yearsite
yearsite_list <- d %>% 
  mutate(ivi = as.double(rep(NA, nrow(d)))) %>%  # empty column to be filled
  mutate(prov_duration = as.double(rep(NA, nrow(d)))) %>%  # empty column to be filled
  group_by(yearsite) %>%
  group_split()

# make sure each yearsite dataframe is descending by date. Very important so we can properly calculate 
# the difference between prov event at time t and prov event at time t-1
for (i in 1:length(yearsite_list)){
  yearsite_list[[i]] <- yearsite_list[[i]] %>% arrange(date_start)
}



# loop through the list of yearsite dataframes, subtract start time from the end of the previous prey delivery    
for (i in 1:length(yearsite_list)){
  for (j in 2:nrow(yearsite_list[[i]])){
    
    # IVI
    current_start <- paste(yearsite_list[[i]][j, "date"], yearsite_list[[i]][j, "start"]) # beginning of current prov event
    prev_end <- paste(yearsite_list[[i]][j-1, "date"],  yearsite_list[[i]][j-1, "start"]) # beginning of previous prov ecent
    yearsite_list[[i]][j, "ivi"] <- abs(as.double(difftime(as.POSIXct(prev_end), 
                                                           as.POSIXct(current_start), 
                                                           units = 'mins')))

    # Provisioning duration following the ivi
    prov_start <- paste(yearsite_list[[i]][j, "date"], yearsite_list[[i]][j, "start"])
    prov_end   <- paste(yearsite_list[[i]][j, "date"],  yearsite_list[[i]][j, "filled_end"])
    yearsite_list[[i]][j, "prov_duration"] <- ifelse(str_detect(prov_end, "NA"), NA,
                                                    as.double(difftime(as.POSIXct(prov_end),
                                                                        as.POSIXct(prov_start),
                                                                        units = 'mins')))
    
    # some errors in the data. End time is before the start time which produces a negative duration. 
    # change these to NA
    yearsite_list[[i]][j, "prov_duration"] <- ifelse(yearsite_list[[i]][j, "prov_duration"] < 0, 
                                                     NA, yearsite_list[[i]][j, "prov_duration"])
  }
}

d <- bind_rows(yearsite_list)



range(d$prov_duration, na.rm = TRUE)
d <- d %>% filter(!prov_duration > 30)
d <- d %>% filter(!ivi > 1440)

hist(log(d$ivi))
