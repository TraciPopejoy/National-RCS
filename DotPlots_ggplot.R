# converting dot plot to ggplot

# Read in RCS data
AnuranNatRCS <- read.csv("C:/Users/Owner/Downloads/AOO HUC12 Output_20200720.csv")%>%
  mutate(avgRCS=mean(RCS_WS, RCS_BUF))
head(AnuranNatRCS)

#adding a species name factor that is ordered based on RCS rank
species_ordered_RCS <- AnuranNatRCS[order(AnuranNatRCS$mean_Rank, decreasing = T),]$scientific_name
AnuranNatRCS<-AnuranNatRCS %>% mutate(SpFac=factor(scientific_name,
                                     levels=species_ordered_RCS))
AnuranNatRCS$scientific_name==AnuranNatRCS$SpFac

# RCS Dotplot code
pdf(paste0("RCS_DotPlot_", format(Sys.Date(), '%Y%m%d'),".pdf"), width=10, height=20)

ggplot()+
  geom_point(data=AnuranNatRCS, 
             aes(x=1-StandardizedMean, y=SpFac))+
  scale_x_reverse(name="Area of Occurrence Index")+
  scale_y_discrete(name="")

ggplot()+
  geom_point(data=AnuranNatRCS, 
             aes(x=logHUC12, y=SpFac,color="HUC 12"))+
  geom_point(data=AnuranNatRCS, 
             aes(x=logBuffer1km, y=SpFac, color="1km buffer"))+
  scale_x_continuous(name="log(Area of Occurrence)")+
  scale_y_discrete(name="")+
  scale_color_viridis_d(name="Grain Size",option="plasma")+
  theme_cowplot()+
  theme(panel.grid.major.y = element_line(color="lightgrey"),
        axis.text.y = element_text(size=10),
        axis.title.x = element_text(size=15))
ggsave('C:/Users/Owner/Desktop/Git/National-RCS/rcs_results/area_plot_jpg.jpg',
       width=7, height=15)

# Make sure species are sorted by RCS values small:large (choose only one grain size). 
Fishes <- Fishes[order(Fishes[,22]),] #sorted by RCS_WS
Fishes$Y <- rank(Fishes$RCS_WS)
rarecol <- addalpha("dodgerblue1", 0.4)
rarecol2 <- addalpha("orange", 0.4)
dotchart(Fishes[,22], labels = Fishes[,1], xlim=c(1,-0.03), font=3, cex=0.8, pch=21, bg=rarecol, xlab = "RCS Value")
points(Fishes[1:144,23], Fishes[1:144,24], cex=0.8, pch=21, bg=rarecol2, type = "p")
dev.off()
