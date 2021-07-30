#overlayed plot
pre_treat <- series %>%
  filter(year <=2019) %>%
  group_by(weekofyr) %>%
  summarize(y = mean(assault_incid_c, na.rm = T))

post_treat <- series %>%
  filter(year==2020) %>%
  group_by(weekofyr) %>%
  summarize(y = mean(assault_incid_c, na.rm = T))

ggplot()+
  geom_line(data = pre_treat, aes(x=weekofyr, y=y, color="2016-2019 Average"))+ 
  geom_line(data = post_treat, aes(x=weekofyr, y=y, color = "2020"))+
  scale_color_manual(name = "Series", values = c("2016-2019 Average"="darkblue","2020"="red"))+
  geom_vline(xintercept=isoweek("2020-05-25"), 
             linetype="dotted", color="black", size=1)+
  geom_label(aes(x=isoweek("2020-05-25")-7, y=0.050), 
             label = "George Floyd", show.legend = FALSE)+
  labs(title = "Figure X: Weekly Hospital Firearm Assaults, 2016-2020",
       subtitle = "MHA Hospital Data ",
       x = "Week",
       y = "Rate per 1,000 Residents")+
  theme_minimal()