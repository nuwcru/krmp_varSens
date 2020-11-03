library(tidyverse)
library(lubridate)

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
x <- d %>% mutate(ivi = as.numeric(rep(NA, nrow(d)))) %>% group_split(yearsite)
          
# loop through the list of yearsite dataframes, subtract start time from the end of the previous prey delivery    
for (j in 1:length(x)){
  for (i in 2:nrow(x[[j]])){
    x[[j]][i, "ivi"] <- as.numeric(x[[j]][i, "date_start"] - x[[j]][i-1, "date_filled_end"])
  }
}

# collapse list
d <- bind_rows(x)
d %>% filter(ivi < 0) %>% dplyr::select(site, yearsite, date_start, date_filled_end, ivi)

view <- d %>% filter(site == 1 & month(date) == 8)

View(view)
