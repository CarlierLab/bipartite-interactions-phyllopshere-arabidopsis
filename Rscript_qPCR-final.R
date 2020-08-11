library(readxl)
library(ggplot2)
qPCR_for_R <- read_excel("C:/Users/acarlier/Desktop/R_bipartite/qPCR_for_R4.xlsx")
qPCR_for_R$delta <- qPCR_for_R$Ct_target - qPCR_for_R$Ct_BKL

plot_box <- function(df_in){
df_ref <- subset(df_in,condition=="TSA-bkl1")
df_target <- subset(df_in,condition=="10dTSA-bkl1")
deltas_ref <- aggregate(df_ref$delta,list(df_ref$combination),mean)
names(deltas_ref)[names(deltas_ref)=='Group.1'] <- 'combination'
names(deltas_ref)[names(deltas_ref)=='x'] <- 'delta_ref'
df_results <- merge(df_target,deltas_ref,by='combination')
df_results$deltadelta <- 2^(-(df_results$delta - df_results$delta_ref))
df_results$mean <- ave(df_results$deltadelta, as.factor(df_results$combination), FUN=mean)
df_results$score[df_results$p.value >0.05]  <- "not significant"
df_results$score[df_results$p.value <0.05]  <- "significant" 
return(df_results)
}

MannU <- function(df_in,ntarget){
  df_sub <- subset(df_in,combination==ntarget)
  t <- t.test(delta ~ condition, data=df_sub)
  return(t$p.value)
}

#df_exp1 <- subset(qPCR_for_R,experiment_number==1)
df_exp2 <- subset(qPCR_for_R,experiment_number==2)
df_exp2$p.value = 0

# for (i in unique(df_exp1$combination)){
#   # df_exp1[df_exp1$combination == i,]$p.value <- (MannU(df_exp1,i))
#   df_exp1[df_exp1$combination == i,]$p.value <- 1
# }

# for(i in unique(df_exp1$combination)){
#   for(row in 1:nrow(df_exp1)){
#     if(df_exp1[row,]$combination == i){
#       df_exp1[row,]$p.value <- MannU(df_exp1,i)}}
#   
# }

for(k in unique(df_exp2$combination)){
  for(row in 1:nrow(df_exp2)){
  if(df_exp2[row,]$combination == k){
    df_exp2[row,]$p.value <- MannU(df_exp2,k)}}
    
}



#df_res_exp1 <- plot_box(df_exp1)
df_res_exp2 <- plot_box(df_exp2)

#df_res <- rbind(df_res_exp1,df_res_exp2)
df_res <-df_res_exp2

p <- ggplot(df_exp2,aes(target,delta,color=as.factor(condition))) + theme_bw() + theme(axis.text.x = element_text(angle = -90,vjust=0.5))   + geom_point(size=3, alpha = 0.5) + facet_grid(as.factor(df_exp2$focal)) + scale_shape(solid=FALSE)
plot(p)

