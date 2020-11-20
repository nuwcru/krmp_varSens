library(tidyverse)
library(lubridate)
library(nuwcru)

# come up with a metric for energy delivered
# duration of feeding event? model with data that exists for prey remaining
 # lots of unknown in amount remaining
  # mode

d <- read_csv("data/Clean IVI years 2013-2019.csv") 
glimpse(d)

# convert to dates
d$date <- dmy(d$date)
d$date_start <- ymd_hms(paste(d$date, d$start))
d$date_end <- ymd_hms(paste(d$date, d$end))
d$date_filled_end <- ymd_hms(paste(d$date, d$filled_end))

# make sure dates are descending
d <- d %>% arrange(site, date)

# split dataframe into lists based on the yearsite
yearsite_list <- d %>% mutate(ivi = as.numeric(rep(NA, nrow(d)))) %>% group_split(yearsite)
          
# loop through the list of yearsite dataframes, subtract start time from the end of the previous prey delivery    
for (i in 1:length(yearsite_list)){
  for (j in 2:nrow(yearsite_list[[i]])){
    yearsite_list[[i]][j, "ivi"] <- as.numeric(yearsite_list[[i]][j, "date_start"] - yearsite_list[[i]][j-1, "date_filled_end"])
  }
}

# this is all good, but the amount of time that ellapses between feedings is likely dependant
# on how much was fed. If we a massive brunch, we likely won't have to eat again until supper
# It could appear that some parents aren't able to provision well because they have large IVI's,
# but maybe it's because they're consistently feeding more food each time (larger prey etc.). 

# collapse list
d <- bind_rows(yearsite_list)
d %>% filter(ivi < 0) %>% dplyr::select(site, yearsite, date_start, date_filled_end, ivi)

view <- d %>% filter(site == 1 & month(date) == 8)

View(view)
