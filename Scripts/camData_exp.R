install.packages("glue")
library(ggthemes)
library(glue)
library(tidyverse)
library(lubridate)
library(stringr)
library(gridExtra)


# plot theme --------------------------------------------------------------

theme_nuwcru <- function(){
  theme_bw() +
    theme(axis.text = element_text(size = 10),
          axis.title = element_text(size = 14),
          axis.text.x = element_text(angle = -45, hjust = -0.05),
          axis.line.x = element_line(color = "black"),
          axis.line.y = element_line(color = "black"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),
          plot.title = element_text(size = 15, vjust = 1, hjust = 0),
          legend.text = element_text(size = 10),
          legend.title = element_blank(),
          legend.position = c(0.9, 0.9),
          legend.key = element_rect(colour = NA, fill = NA),
          legend.background = element_rect(color = "black",
                                           fill = "transparent",
                                           size = 4, linetype = "blank"))
      }

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


# Visualize Prey Deliveries -----------------------------------------------
# **daily deliveries/chick-------------------------
delDay <- dat %>% group_by(site, Date, Chicks, treatment) %>% tally() %>% arrange(site, Date) %>% mutate(pChick = n/Chicks)

dNested <- delDay %>% mutate(year = year(Date)) %>% group_by(year) %>% nest()
  
  
d_plots <- dNested %>% 
  mutate(plot = map2(data, year, ~ ggplot(data = .x, aes(x = Date, y = pChick, color = treatment)) +
                                         ggtitle(glue("{.y}")) + ylab("deliveries per day per nestling") + 
                                              geom_point(alpha = 0.4)))


# visualize number of deliveries per day per chick from 2015:2017
p1 <- d_plots$plot[[1]] + scale_color_manual(values = c("#999999", "#E69F00")) + geom_smooth(method = "loess")  + 
  theme_nuwcru() + theme(axis.title.x = element_blank())
p2 <- d_plots$plot[[2]] + scale_color_manual(values = c("#999999", "#E69F00")) + geom_smooth(method = "loess") + theme_nuwcru() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())
p3 <- d_plots$plot[[3]] + scale_color_manual(values = c("#999999", "#E69F00")) + geom_smooth(method = "loess") + theme_nuwcru() + 
  theme(axis.title.x = element_blank(),axis.title.y = element_blank()) 
grid.arrange(p1,p2,p3, nrow = 1)

# **daily deliveries per nest ---------------------------------------------


delDay <- dat %>% group_by(site, Date,treatment) %>% tally() %>% arrange(site, Date) 

dNested <- delDay %>% mutate(year = year(Date)) %>% group_by(year) %>% nest()


d_plots <- dNested %>% 
  mutate(plot = map2(data, year, ~ ggplot(data = .x, aes(x = Date, y = n, color = treatment)) +
                       ggtitle(glue("{.y}")) + ylab("deliveries per day per nest") + 
                       geom_point()))

# visualize number of deliveries per day per chick from 2015:2017
p1 <- d_plots$plot[[1]] + scale_color_manual(values = c("#999999", "#E69F00")) + geom_smooth(method = "loess")  + 
  theme_nuwcru() + theme(axis.title.x = element_blank())
p2 <- d_plots$plot[[2]] + scale_color_manual(values = c("#999999", "#E69F00")) + geom_smooth(method = "loess") + theme_nuwcru() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())
p3 <- d_plots$plot[[3]] + scale_color_manual(values = c("#999999", "#E69F00")) + geom_smooth(method = "loess") + theme_nuwcru() + 
  theme(axis.title.x = element_blank(),axis.title.y = element_blank()) 
grid.arrange(p1,p2,p3, nrow = 1)


# chick numbers across time -----------------------------------------------
# not the best data set to visualize brood size


#nDay <- dat %>% select(site, Date,treatment, Chicks) %>% group_by(Date, treatment) %>%
#  summarize(mean = mean(Chicks)) %>%arrange(Date) 
#
#dNested <- nDay %>% mutate(year = year(Date)) %>% group_by(year) %>% nest()
#
#d_plots <- dNested %>% 
#  mutate(plot = map2(data, year, ~ ggplot(data = .x, aes(x = Date, y = mean, color = treatment)) +
#                       ggtitle(glue("{.y}")) + ylab("mean brood size") + 
#                       geom_point()))
#
#
## visualize number of deliveries per day per chick from 2015:2017
#p1 <- d_plots$plot[[1]] + scale_color_manual(values = c("#999999", "#E69F00")) + geom_smooth(method = "loess")  + 
#  theme_nuwcru() + theme(axis.title.x = element_blank())
#p2 <- d_plots$plot[[2]] + scale_color_manual(values = c("#999999", "#E69F00")) + geom_smooth(method = "loess") + theme_nuwcru() +
#  theme(axis.title.x = element_blank(), axis.title.y = element_blank())
#p3 <- d_plots$plot[[3]] + scale_color_manual(values = c("#999999", "#E69F00")) + geom_smooth(method = "loess") + theme_nuwcru() + 
#  theme(axis.title.x = element_blank(),axis.title.y = element_blank()) 
#grid.arrange(p1,p2,p3, nrow = 1)
#




# visualize guild proportions per day across time -------------------------



levels(as.factor(dat$BestID))

