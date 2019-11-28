install.packages("glue")
library(glue)
library(tidyverse)
library(lubridate)
library(stringr)
library(gridExtra)

# Load data ---------------------------------------------------------------
dat <- read_csv("Data/provisioning_KAH.csv")
trt <- read_csv("Data/treatment.csv")

# survival data - only using this to pull info on what sites were supplemented
  trt <- trt[,-1]
  names(trt) <- c("nest", "year", "treatment")
  trt <- filter(trt, year %in% c(2015, 2016, 2017)) # match the years with provisioning data
  trt <- rename(trt, site = nest)

# provisioning data
  dat$Date <- dmy(dat$Date)
  dat$year <- year(dat$Date)
  dat <- rename(dat, site = Site)

dat <- left_join(dat, trt, by = c("site", "year"))

# some NAs in treatment column, because nests failed before treatment could be allocated. ie they default to control
# fill in NAs with 0s
dat$treatment[is.na(dat$treatment)] <- 0
table(is.na(dat$treatment)) # make sure NAs are dealt with

dat$treatment <- as.factor(dat$treatment)

# daily deliveries/chick-------------------------
delDay <- dat %>% group_by(site, Date, Chicks, treatment) %>% tally() %>% arrange(site, Date) %>% mutate(pChick = n/Chicks)

dNested <- delDay %>% mutate(year = year(Date)) %>% group_by(year) %>% nest()
  
  
d_plots <- dNested %>% 
  mutate(plot = map2(data, year, ~ ggplot(data = .x, aes(x = Date, y = pChick, color = treatment)) +
                                         ggtitle(glue("{.y}")) + ylab("deliveries per day per nestling") + 
                                              geom_point()))


# visualize number of deliveries per day per chick from 2015:2017
p1 <- d_plots$plot[[1]] + scale_color_manual(values = c("#999999", "#E69F00")) + geom_smooth(method = "loess") + theme_bw()
p2 <- d_plots$plot[[2]] + scale_color_manual(values = c("#999999", "#E69F00")) + geom_smooth(method = "loess") + theme_bw()
p3 <- d_plots$plot[[3]] + scale_color_manual(values = c("#999999", "#E69F00")) + geom_smooth(method = "loess") + theme_bw()
grid.arrange(p1,p2,p3, nrow = 2)







sum %>% group_by(Site, Date) %>% tally() %>% arrange(Site, Date) %>% filter(Site == "1" & year(Date) == "2015") %>%
  ggplot(aes(x = Date, y = n)) + geom_point() + geom_smooth()

?facet_grid


# visualize guild proportions per day across time -------------------------





