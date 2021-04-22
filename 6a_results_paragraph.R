# Summary values for the Results section
RCS_Data <- read.csv('rcs_results/RCS_table_20210406.csv')
names(RCS_Data)
RCS_Data %>% select(RCS_WS, RCS_buff) %>% 
  summarise(quantile(RCS_WS), quantile(RCS_buff)) 

RCS_Data %>% 
  summarize(median(buffer), range(buffer)[1], range(buffer)[2],
            median(watershed), range(watershed)[1], range(watershed)[2],
            median(AOO_BUF_adj), median(AOO_WS_adj)) #%>%  summarize_all(log10)
#how to describe climate niche breadth?

#taxa that have large residuals between variable CS & calc CS
RCS_Data %>% as_tibble() %>%
  select(species, ends_with('CS')) %>%
  group_by(species,buff_CS,WS_CS) %>%
  summarize(ppt.buf.res=buff_CS-ppt_buf.CS,
            tmax.buf.res=buff_CS-tmax_buf.CS,
            tmin.buf.res=buff_CS-tmin_buf.CS,
            ppt.WS.res=WS_CS-ppt_WS.CS,
            tmax.WS.res=WS_CS-tmax_WS.CS,
            tmin.WS.res=WS_CS-tmin_WS.CS,
            .groups='keep') %>%
  pivot_longer(cols=c(-species, -buff_CS, -WS_CS)) %>%
  group_by(name) %>%
  filter(value==min(value) | value==max(value)) %>%
  arrange(name) %>%
  pivot_wider()
    
RCS_Data %>% as_tibble() %>%
  select(species, ends_with('CS')) %>%
  group_by(species,buff_CS,WS_CS) %>%
  summarize(ppt.buf.res=buff_CS-ppt_buf.CS,
            tmax.buf.res=buff_CS-tmax_buf.CS,
            tmin.buf.res=buff_CS-tmin_buf.CS,
            ppt.WS.res=WS_CS-ppt_WS.CS,
            tmax.WS.res=WS_CS-tmax_WS.CS,
            tmin.WS.res=WS_CS-tmin_WS.CS,
            .groups='keep') %>%  
  pivot_longer(cols=c(-species, -buff_CS, -WS_CS)) %>%
  group_by(name) %>%
  mutate(val_typ=substr(name, 1,4),
         medina=median(value)) %>%
  #summarize(mean(value)) %>%
  ggplot()+geom_density(aes(x=value, fill=name), alpha=0.3) +
  geom_vline(aes(xintercept=medina))+
  facet_wrap(~val_typ)+
  theme(legend.position = 'top')

# mean standard dev of climate variables
climate_ws_raw %>%
  group_by(species, value_type) %>% 
  filter(!(species %in% c('Anaxyrus baxteri', 'Anaxyrus houstonensis'))) %>%
  dplyr::summarize(mean_climate=mean(value), sd_of_env_mean=sd(value)) %>%
  group_by(value_type) %>%
  summarize(range(sd_of_env_mean))

#taxa that had largest and smallest aoo respectively
RCS_Data %>%
  select(species, ends_with('adj')) %>%
  pivot_longer(-species) %>%
  group_by(name) %>%
  filter(value == min(value) | value == max(value)) %>%
  arrange(name, value)

# graph like fig 2 in mims et al 2018
RCS_Data %>% as_tibble() %>%
  select(species, ends_with('adj'), starts_with('RCS')) %>%
  pivot_longer(-species) %>%
  mutate(grain.size=ifelse(grepl('WS', name), 'watershed', 'buffer'),
         new_name=recode(name, AOO_WS_adj='AOO_adj', CS_WS_adj='CS_adj',
                         AOO_BUF_adj='AOO_adj', CS_BUF_adj='CS_adj',
                         RCS_WS='RCS',RCS_buff='RCS')) %>%
  select(-name) %>%
  pivot_wider(names_from=new_name, values_from=value) %>%
  ggplot()+
  geom_point(aes(x=AOO_adj, y=CS_adj))+
  facet_wrap(~grain.size)+
  scale_x_reverse()+
  theme(legend.position='top')

?cor.test
cor.test(RCS_Data$RCS_WS, RCS_Data$RCS_buff)
