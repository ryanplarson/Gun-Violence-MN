```{r, eval=F, message=F, eval=F, include=F}
#weatherData
#the scraper sometimes doesn't run, so I have saved successful scrapes locally for backup
noaa_token <- "nrEVQsKalyfOauURFnZPpNZQwbAiPaSh"

#list of ncdc stations with coverage over MPLS ZCTAs
stations <- ncdc_stations(extent = as.vector(st_bbox(zcta[zcta$zcta %in% zcta_universe,]))[order(c(2,1,4,3))],
                          token = noaa_token, startdate = "2016-01-01", enddate = "2020-12-31") 

ncdc_datasets(stationid=stations$data$id, token=noaa_token)

datatypes <- c('PRCP','SNWD', 'TMAX') #can't get SNOW to work

begin_date <- as.list(series$begin_date)
end_date <- as.list(series$end_date)

weather <- map2_df(.x = begin_date,
                   .y = end_date,
                   ~ ncdc(datasetid='GHCND', 
                          datatypeid=datatypes,
                          stationid=ncdc_stations(extent = as.vector(st_bbox(zcta[zcta$zcta %in% zcta_universe,]))[order(c(2,1,4,3))],
                                                  token = noaa_token, startdate = .x, enddate = .y)$data$id,  
                          startdate = .x, 
                          enddate = .y, 
                          limit=1000, 
                          token = noaa_token))$data %>%
  group_by(date, datatype) %>%
  summarize(value = mean(value, na.rm = T)) %>%
  pivot_wider(id_cols = date, names_from = datatype, values_from = value))


write_csv(weather, "noaa_scrape.csv")

weather <- read_csv("noaa_scrape.csv")


weather_series <- weather %>% 
  mutate(precip.in = (PRCP/10)*.0393701, #mm ti in
         snow.in = ifelse(is.na(SNWD), 0, SNWD/10)*.0393701, #mm to in
         tmax.f = TMAX/10, 
         week = isoweek(date),
         year = year(date)) %>%
  select(year, week, precip.in, snow.in, tmax.f) 

write_csv(weather_series, "weather_series.csv")

weather_series <- read_csv("weather_series.csv")


#join to series
series <- series %>% left_join(weather_series, by = c("year","weekofyr"="week"))
```