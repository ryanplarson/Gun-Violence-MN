#MPD Crime Incidents

library(tidyverse)

crime_2019 <- read_csv("Police_Incidents_2019.csv") %>%
  mutate(date=ymd_hms(reportedDate),
         year=year(date),
         week=isoweek(date)) %>%
  group_by(week) %>%
  tally() %>%
  mutate(year=2019,
         week_id = row_number(),
         begin_date = ISOweek2date(paste(year, paste0("W", sprintf("%02d", week)), 1,sep = "-")))


crime_2020 <- read_csv("Police_Incidents_2020.csv") %>%
  mutate(date=ymd_hms(reportedDate),
         year=year(date),
         week=isoweek(date)) %>%
  group_by(week) %>%
  tally() %>%
  mutate(year=2020) %>%
  filter(week!=53)

crime_2021<- read_csv("Police_Incidents_2021.csv") %>%
  mutate(date=ymd_hms(reportedDate),
         year=year(date),
         week=isoweek(date)) %>%
  group_by(week) %>%
  tally() %>%
  mutate(year=2021) %>%
  filter(week <= 23)

crime_series <- crime_2020 %>%
  rbind(crime_2021) %>%
  arrange(year, week) %>%
  mutate(week_id = row_number(),
         begin_date = ISOweek2date(paste(year, paste0("W", sprintf("%02d", week)), 1,sep = "-")))

ggplot(crime_series)+
  geom_line(aes(x=begin_date, y=n))+ 
  geom_line(data = crime_2019, aes(x = begin_date, y=n), color = "red")+
  scale_x_date(date_labels = "%b-%Y", date_breaks = "6 months")+
  geom_vline(xintercept=crime_series$begin_date[crime_series$year==2020 & 
                                                  crime_series$week==isoweek(date("2020-05-25"))], 
             linetype="dotted", color="red", size=1)+
  geom_label(aes(x=crime_series$begin_date[crime_series$year==2020 & crime_series$week==isoweek(date("2020-05-25"))],
                 y=200), 
             label = "George Floyd", show.legend = FALSE)+
  labs(title = "Weekly Crime Incidents, 2019-2021",
       subtitle = "Source: MPD",
       x = "Week",
       y = "Total Crime Incidents")+
  theme_minimal()+
  theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave('crime_series.png')

#homicide chart
crime_2019_h <- read_csv("Police_Incidents_2019.csv") %>%
  mutate(date=ymd_hms(reportedDate),
         year=year(date),
         week=isoweek(date)) %>%
  filter(offense=="MURDR") %>%
  group_by(week) %>%
  tally() %>%
  mutate(year=2019) %>%
  mutate(week_id = row_number(),
         begin_date = ISOweek2date(paste(year, paste0("W", sprintf("%02d", week)), 1,sep = "-")))


crime_2020_h <- read_csv("Police_Incidents_2020.csv") %>%
  mutate(date=ymd_hms(reportedDate),
         year=year(date),
         week=isoweek(date)) %>%
  filter(offense=="MURDR") %>%
  group_by(week) %>%
  tally() %>%
  mutate(year=2020) %>%
  filter(week!=53)

crime_2021_h <- read_csv("Police_Incidents_2021.csv") %>%
  mutate(date=ymd_hms(reportedDate),
         year=year(date),
         week=isoweek(date)) %>%
  filter(offense=="MURDR") %>%
  group_by(week) %>%
  tally() %>%
  mutate(year=2021) %>%
  filter(week <= 23)

crime_series_h <- crime_2020_h %>%
  rbind(crime_2021_h) %>%
  arrange(year, week) %>%
  mutate(week_id = row_number(),
         begin_date = ISOweek2date(paste(year, paste0("W", sprintf("%02d", week)), 1,sep = "-")))

ggplot(crime_series_h)+
  geom_line(aes(x=begin_date, y=n))+ 
  geom_line(data = crime_2019_h, aes(x = begin_date, y=n), color = "red")+
  geom_smooth(aes(x = begin_date, y=n), se=F)+
  geom_smooth(data = crime_2019_h, aes(x = begin_date, y=n), se=F)+
  scale_x_date(date_labels = "%b-%Y", date_breaks = "6 months")+
  geom_vline(xintercept=crime_series$begin_date[crime_series$year==2020 & 
                                                  crime_series$week==isoweek(date("2020-05-25"))], 
             linetype="dotted", color="red", size=1)+
  geom_label(aes(x=crime_series$begin_date[crime_series$year==2020 & crime_series$week==isoweek(date("2020-05-25"))],
                 y=0.25), 
             label = "George Floyd", show.legend = FALSE)+
  labs(title = "Weekly Minneapolis Homicides, 2019-2021",
       subtitle = "Source: MPD",
       x = "Week",
       y = "Homicides")+
  theme_minimal()+
  theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave('homicide_series.png')