# spider_dat<-matrix(c(50,2.19,1,
#                      25,2.19,1,
#                      1,1,10.4),3,3)
# colnames(spider_dat)<-c("HIV Positive","Diabetic","Non-US Born")
# rownames(spider_dat)<-c("TB Progression","Mortality", "LTBI Prevalence")
#
#
#
# spider_dat<-as.data.frame(spider_dat)
# spider_dat <- as.data.frame(t(spider_dat))
#
#
# spider_dat %>% add_rownames( var = "group" ) %>%
#   mutate_each(funs(rescale), -group) -> spider_radar
#
# # ggradar(spider_radar, grid.label.size=0, legend.position = 'none', axis.label.size = 3) + facet_wrap(~group)
#
# HIV<-ggradar(spider_radar[1,], grid.label.size=0, legend.position = 'none', axis.label.size = 3, font.radar = "Helvetica") + facet_wrap(~group)
#
# diabetes<-ggradar(spider_radar[2,], grid.label.size=0, legend.position = 'none', axis.label.size = 3, font.radar = "Helvetica") + facet_wrap(~group)
#
# nusb<-ggradar(spider_radar[3,], grid.label.size=0, legend.position = 'none', axis.label.size = 3, font.radar = "Helvetica") + facet_wrap(~group)
# pdf("~/Documents/Conferences/NTCA 2019/spiders.pdf")
# HIV
# diabetes
# nusb
# dev.off()
#
