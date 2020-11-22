library(tidyverse)
library(lubridate)

d <- read_csv("Data/erik.csv")
a <- read_csv("Data/chick_weights.csv") 



d <- d %>% 
        mutate(year = year(HatchDate)) %>%
        mutate(id = paste0(CMChickID, "_", year, "_", NestSiteID)) %>%
        mutate(u_measure = paste0(id, "_", ChickMeasurementDate)) %>%
        filter(!duplicated(u_measure)) %>%
        rename("site" = NestSiteID) %>%
  mutate(hatch = yday(HatchDate),
         yday = yday(ChickMeasurementDate)) %>%
  mutate(age = yday - hatch )
d

x <- tibble(id    =d$id, 
           site  =d$site,
           year  =d$year,
           hatch =d$hatch,
           alive =rep(1, nrow(d)),
           yday  =d$yday,
           age   =d$age,
           m_weight = rep(NA, nrow(d))) %>%
  filter(!is.na(age)) %>% filter(year != 2017)
 
x <- x %>% group_by(site, year) %>% tally() %>% 
  mutate(id = paste0(year, "_", site))
  
d <- d %>% filter( year %in% c("2018", "2019")) %>%
  group_by(year, site) %>% tally() %>%
  mutate(id = paste0(year, "_", site))

setdiff(d$id, x$id)  
x %>% filter(id == "2018_30")
