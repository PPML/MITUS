# Y <- rnorm(100)
# df <- data.frame(A = rnorm(100), B = runif(100), C = rlnorm(100),
#                  Y = Y)
# colNames <- names(df)[1:3]
# for(i in colNames){
#   plt <- ggplot(df, aes_string(x=i, y = Y)) +
#     geom_point(color="#B20000", size=4, alpha=0.5) +
#     geom_hline(yintercept=0, size=0.06, color="black") +
#     geom_smooth(method=lm, alpha=0.25, color="black", fill="black")
#   print(plt)
#   Sys.sleep(2)
# }
#
#
#
#
# DxUSPlot<-ggplot(outcome_df,aes(x=rownames(outcome_df), y=DxUS))+
#   geom_point(stat="identity", size=4, alpha=0.5)+theme(axis.title.y = element_blank())+coord_flip()
#
#
# #create a new column for popnames
#   # jj<-as.data.frame(outcome_df[,i])
#   rownames(jj)<-rownames(outcome_df)
#  plt<- ggplot(outcome_df, aes(x='pop name',y=outcome_df[,i]))+
#     geom_point(color="#B20000", size=4, alpha=0.5)
#  print(plt)
#
# for (i in 1:10){
# #create a new df
# temp_df<-as.data.frame(outcome_df[,i])
# # colnames(temp_df)<-colnames(outcome_df)[i]
# temp_df$'Population'<-rownames(outcome_df)
# temp_df$type<-ifelse(temp_df$`outcome_df[, i]`< 0, "below", "above")
# temp_df <- temp_df[order(temp_df$type), ]  # sort
# temp_df$'Population'<-'Population'<-factor(temp_df$'Population', levels=temp_df$'Population')
#
# ggplot(temp_df, aes(x='Population', y=outcome_df[,i], label=outcome_df[,i]))+
#   geom_point(stat = 'identity', aes(col=type), size=6)+
#   geom_text(color="white", size=2) +
#   ylim(-limits*1.1,limits*1.1) +
#   coord_flip()
# }
#
#  temp_df<-as.data.frame(outcome_df[,i])
#  # colnames(temp_df)<-colnames(outcome_df)[i]
#  temp_df$'Population'<-rownames(outcome_df)
#  temp_df$type<-ifelse(temp_df$`outcome_df[, i]`< 0, "below", "above")
#  temp_df <- temp_df[order(temp_df$type), ]  # sort
#  temp_df$'Population'<-'Population'<-factor(temp_df$'Population', levels=temp_df$'Population')
#
#  # limits<-max(abs(min(na.omit(temp_df[,1]))), abs(max(na.omit(temp_df[,1]))))
#  ggplot(outcome_df, aes(x=colnames(outcome_df), y="DxUs"))+
#    geom_point(stat = 'identity', size=6)+
#    # geom_text(color="white", size=2) +
#    # ylim(-limits*1.1,limits*1.1) +
#    coord_flip()
#
#
#
#
# outcome_df$Population<-rownames(outcome_df)
# m_outcome_df<-melt(outcome_df,id.vars='Population', variable.name = "outcome")
#
# ggplot(m_outcome_df,aes(x=Population, y=value, label=value))+
#   geom_point(stat="identity", size=6,alpha=1)+theme(axis.title.y = element_blank())+coord_flip()+
#   geom_text(color="white", size=2) +
#   facet_wrap(~outcome,ncol=8, scales = "free_x")
#
#   # ylim(-limits*1.1,limits*1.1)+
#
#
#
#  require(ggplot2)
# require(reshape2)
# df <- data.frame(time = 1:10,
#                  a = cumsum(rnorm(10)),
#                  b = cumsum(rnorm(10)),
#                  c = cumsum(rnorm(10)))
# mdf <- melt(df ,  id.vars = 'time', variable.name = 'series')
#
# # plot on same grid, each series colored differently --
# # good if the series have same scale
# ggplot(df, aes(time,value)) + geom_line(aes(colour = series))
#
# # or plot on different plots
# ggplot(df, aes(time,value)) + geom_line() + facet_grid(series ~ .)
#
#
#
#
#
