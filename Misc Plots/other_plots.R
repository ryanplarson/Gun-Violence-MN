```{r firearm combined plot, warning=FALSE}
ggplot(series)+
  geom_line(aes(x=begin_date, y=combined_incid_c))+ 
  scale_x_date(date_labels = "%b-%Y", date_breaks = "6 months")+
  geom_vline(xintercept=series$begin_date[series$year==2020 & series$weekofyr==isoweek(date("2020-05-25"))], 
             linetype="dotted", color="red", size=1)+
  geom_label(aes(x=series$begin_date[series$year==2020 & series$weekofyr==isoweek(date("2020-05-25"))],
                 y=0.033), 
             label = "George Floyd", show.legend = FALSE)+
  labs(title = "Weekly Firearm Cases, 2016-2020",
       subtitle = "Source: Minnesota Hospital Association Discharges",
       x = "Week",
       y = "Firearm Cases per 1,000 Residents")+
  theme_minimal()+
  theme(axis.text.x=element_text(angle=45, hjust=1)) 
```


```{r use of force plot, warning=FALSE}
ggplot(series)+
  geom_line(aes(x=begin_date, y=use_of_force_rate))+ 
  scale_x_date(date_labels = "%b-%Y", date_breaks = "6 months")+
  geom_vline(xintercept=series$begin_date[series$year==2020 & series$weekofyr==isoweek(date("2020-05-25"))], 
             linetype="dotted", color="red", size=1)+
  geom_label(aes(x=series$begin_date[series$year==2020 & series$weekofyr==isoweek(date("2020-05-25"))],
                 y=0.033), 
             label = "George Floyd", show.legend = FALSE)+
  labs(title = "Weekly Uses of Force, 2016-2020",
       subtitle = "Source: MPD",
       x = "Week",
       y = "Uses of Force per 1,000 Residents")+
  theme_minimal()+
  theme(axis.text.x=element_text(angle=45, hjust=1)) 
```

```{r police stops plot, warning=FALSE}
ggplot(series)+
  geom_line(aes(x=begin_date, y=police_stop_rate))+ 
  scale_x_date(date_labels = "%b-%Y", date_breaks = "6 months")+
  geom_vline(xintercept=series$begin_date[series$year==2020 & series$weekofyr==isoweek(date("2020-05-25"))], 
             linetype="dotted", color="red", size=1)+
  geom_label(aes(x=series$begin_date[series$year==2020 & series$weekofyr==isoweek(date("2020-05-25"))],
                 y=0.033), 
             label = "George Floyd", show.legend = FALSE)+
  labs(title = "Weekly Police Stops, 2016-2020",
       subtitle = "Source: MPD",
       x = "Week",
       y = "Stops per 1,000 Residents")+
  theme_minimal()+
  theme(axis.text.x=element_text(angle=45, hjust=1)) 
```


```{r ois plot}
ggplot(series)+
  geom_line(aes(x=begin_date, y=off_inv_shooting_rate))+ 
  scale_x_date(date_labels = "%b-%Y", date_breaks = "6 months")+
  geom_vline(xintercept=series$begin_date[series$year==2020 & series$weekofyr==isoweek(date("2020-05-25"))], 
             linetype="dotted", color="red", size=1)+
  geom_label(aes(x=series$begin_date[series$year==2020 & series$weekofyr==isoweek(date("2020-05-25"))],
                 y=0.033), 
             label = "George Floyd", show.legend = FALSE)+
  labs(title = "Weekly Officer Involved Shootings, 2016-2020",
       subtitle = "Source: MPD",
       x = "Week",
       y = "Stops per 1,000 Residents")+
  theme_minimal()+
  theme(axis.text.x=element_text(angle=45, hjust=1)) 
```
