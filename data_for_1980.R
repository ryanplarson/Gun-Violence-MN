#MN tracts
tracts <- get_acs(geography = "tract", 
                  state = "MN",
                  variables = "B01001_001E",
                  output = "wide",
                  survey = "acs5",
                  year = 2020,
                  geometry = T)

#Minneapolis Shapefile
mpls <- st_read("Data/mpls_city-shp/16cdbbfa-ad10-493c-afaf-52b61f2e76e42020329-1-180h9ap.whbo.shp") %>%
  st_transform(st_crs(tracts))

mpls_tract <- tracts %>%
  st_filter(mpls, .predicate = st_intersects) %>%
  mutate(GEOID = as.numeric(GEOID),
         tract_area = as.numeric(st_area(.)),
         tract_area_sqkm = tract_area*.000001, 
         tract_area_sqmi = tract_area_sqkm*.386102,
         intersection_area = as.numeric(st_area(st_intersection(., mpls))),
         perc_intersection = intersection_area/tract_area*100) %>%
  filter(perc_intersection >= 2) %>%
  select(-"B01001_001M")

homicide_tract <- cj_exp_prepost %>%
  group_by(NAME.x) %>%
  summarize(homicide = mean(homicide),
            total_pop = mean(total_pop.x),
            povlevel = mean(povlevel),
            unemp = mean(unemp),
            total_ilf = mean(total_ilf)) %>%
  rename(NAME = NAME.x)


uof_spatial <- read_csv("Data/Police_Use_Of_Force.csv")  %>%
  mutate(date=ymd_hms(ResponseDate),
         year=isoyear(date),
         week=isoweek(date)) %>%
  filter(year >= 2017 & year <= 2022) %>%
  select(OBJECTID, X, Y) %>%
  st_as_sf(coords = c("X", "Y"), crs = "NAD83", remove=F) %>%
  st_join(mpls_tract) %>%
  st_drop_geometry() %>%
  filter(!is.na(NAME)  & NAME %in% mpls_tract$NAME) %>%
  group_by(NAME, .drop=F) %>%
  tally(name = "use_of_force") %>%
  ungroup() %>%
  complete(NAME=mpls_tract$NAME, fill = list(use_of_force = 0))

ois_spatial <- read_csv("Data/Police_Officer_Involved_Shootings.csv") %>%
  mutate(date=ymd_hms(IncidentDate),
         year=isoyear(date),
         week=isoweek(date)) %>%
  filter(year >= 2017 & year <= 2022) %>%
  select(CenterLatitude, CenterLongitude) %>%
  rename(lat = CenterLatitude,
         long = CenterLongitude) %>%
  st_as_sf(coords = c("long", "lat"), crs = "NAD83", remove=F) %>%
  st_join(mpls_tract) %>%
  st_drop_geometry() %>%
  filter(!is.na(NAME)  & NAME %in% mpls_tract$NAME) %>%
  group_by(NAME, .drop=F) %>%
  tally(name = "police_shootings") %>%
  ungroup() %>%
  complete(NAME=mpls_tract$NAME, fill = list(police_shootings = 0))

homicide_tract <- homicide_tract %>%
  left_join(uof_spatial, by = "NAME") %>%
  left_join(ois_spatial, by = "NAME") %>%
  st_drop_geometry() %>%
  mutate(ps = as.factor(ifelse(police_shootings==0, 0, 1)),
         hr = homicide/total_pop*100000,
         povr = povlevel/total_pop,
         uofr = use_of_force/total_pop*100000,
         unempr = unemp/total_ilf*100)

write.csv(homicide_tract, "C:/Users/rlarson21/Documents/Teaching/CJFS-1980-Basic-Quantitative-Criminology-and-Statistics/Spring 2023/Lectures/Week #15/homicide_tract.csv") 
