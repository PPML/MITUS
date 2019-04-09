library(ggplot2)
library(usmap)


x<-matrix(c(state.abb,"DC",rep("not calibrated",51)),51,2)
for (i in 1:51){
  if (x[i,1] %in% c("MA","CA","TX","NY","NJ","IL", "GA", "VA", "WA", "PA", "FL")){
  x[i,2]<-"calibrated"
}
}
colnames(x)<-c("state","calib")
x<-as.data.frame(x)
plot_usmap(data = x, values = "calib", labels=TRUE) + scale_fill_brewer(name="Calibrated States", type="div", palette = 4)+ theme(legend.position = "right")

# plot_usmap(include=c("MA","CA","TX","NY","NJ","IL", "GA", "VA", "WA", "PA", "FL"), lines="red", )+ labs(title = "MITUS Calibrated States", subtitle = "Top 10 States of TB Incidence + MA")
