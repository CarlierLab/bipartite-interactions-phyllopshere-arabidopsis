library(ggplot2)
library(knitr)

calc_chi <- function(df_in,caption_in) {
  #calculate PP vs SS
  df_in_PPSS <- subset(df_in,df_in$sources %in% c("P_P","S_S"))
  df_in_PPSS <- droplevels(df_in_PPSS)
  df_in_PPPS <- subset(df_in,df_in$sources %in% c("P_P","P_S"))
  df_in_PPPS <- droplevels(df_in_PPPS)
  df_in_PSSS <- subset(df_in,df_in$sources %in% c("S_S","P_S"))
  df_in_PSSS <- droplevels(df_in_PSSS)
  
  tbl_in_PPSS <- table(df_in_PPSS$dummy,df_in_PPSS$sources)
  print(caption_in)
  print(kable(tbl_in_PPSS, caption = caption_in))
  #print("PP vs SS")
  print(chisq.test(tbl_in_PPSS))
  
  tbl_in_PPPS <- table(df_in_PPPS$dummy,df_in_PPPS$sources)
  #print("PP vs PS")
  print(caption_in)
  print(kable(tbl_in_PPPS,caption=caption_in))
  print(chisq.test(tbl_in_PPPS))
  
  tbl_in_PSSS <- table(df_in_PSSS$dummy,df_in_PSSS$sources)
  #print("PS vs SS")
  print(caption_in)
  print(kable(tbl_in_PSSS,caption = caption_in))
  print(chisq.test(tbl_in_PSSS))
  
}


df <- read.csv("bipartite_data.csv", sep=";")
#removing observations where there is no or ambiguous inhibition score
df <- subset(df, df$inhibition.score != "NaN")

df$dummy[df$inhibition.score >0] <- 1
df$dummy[df$inhibition.score == 0] <- 0

print(kable(table(df$sources,df$test.organism)))
df_bi <- subset(df,monoculture==0)
df_mono <- subset(df,monoculture>0)
df_bi_staph <- subset(df_bi,test.organism == "Staphylococcus")
df_bi_acin <- subset(df_bi,test.organism == "Acinetobacter")
df_mono_staph <- subset(df_mono,test.organism == "Staphylococcus")
df_mono_acin <- subset(df_mono,test.organism == "Acinetobacter")

tbl_bi_staph <- table(df_bi_staph$source,df_bi_staph$dummy)
print(kable(tbl_bi_staph,caption = "tbl_bi_staph"))
mosaicplot(data=df_bi_staph, sources ~inhibition.score, col=c("lightgrey","lightskyblue2","tomato"))
tbl_bi_acin <- table(df_bi_acin$source,df_bi_acin$dummy)
print(kable(tbl_bi_acin,caption = "tbl_bi_acin"))
mosaicplot(data=df_bi_acin, sources ~inhibition.score, col=c("lightgrey","lightskyblue2","tomato"))
tbl_mono_staph <- table(df_mono_staph$source,df_mono_staph$dummy)
print(kable(tbl_mono_staph,caption = "tbl_mono_staph"))
mosaicplot(data=df_mono_staph, sources ~inhibition.score, col=c("lightgrey","lightskyblue2","tomato"))
tbl_mono_acin <- table(df_mono_acin$source,df_mono_acin$dummy)
#plot(mosaicplot(data=df_mono_acin, sources ~inhibition.score, col=c("grey","lightskyblue2","tomato"))
print(kable(tbl_mono_acin,caption = "tbl_mono_acin"))

print("bi staph")
calc_chi(df_bi_staph,"bi_staph")
print("bi acinetobacter")
calc_chi(df_bi_acin,"bi_acin")
print("mono acinetobacter")
calc_chi(df_mono_acin,"mono_acin")
print("mono staph")
calc_chi(df_mono_staph,"mono staph")

# plot inhibition by phylogenetic distance

plot(p <- ggplot(data = df_bi_staph, aes(x= as.factor(inhibition.score), y = as.numeric(as.character(phylogenetic.distance..Tamura.Nei.)))) + geom_boxplot() + geom_jitter(alpha=0.1))
plot(q <- ggplot(data = df_bi_acin, aes(x= as.factor(inhibition.score), y = as.numeric(as.character(phylogenetic.distance..Tamura.Nei.)))) + geom_boxplot() + geom_jitter(alpha=0.1))
plot(s <- ggplot(data = df_mono_staph, aes(x= as.factor(inhibition.score), y = as.numeric(as.character(phylogenetic.distance..Tamura.Nei.)))) + geom_boxplot() + geom_jitter(alpha=0.1))
plot(t <- ggplot(data = df_mono_acin, aes(x= as.factor(inhibition.score), y = as.numeric(as.character(phylogenetic.distance..Tamura.Nei.)))) + geom_boxplot() + geom_jitter(alpha=0.1))

kruskal.test(as.numeric(as.character(phylogenetic.distance..Tamura.Nei.)) ~ inhibition.score, data = df_bi_staph)
kruskal.test(as.numeric(as.character(phylogenetic.distance..Tamura.Nei.)) ~ inhibition.score, data = df_bi_acin)
kruskal.test(as.numeric(as.character(phylogenetic.distance..Tamura.Nei.)) ~ inhibition.score, data = df_mono_staph)
kruskal.test(as.numeric(as.character(phylogenetic.distance..Tamura.Nei.)) ~ inhibition.score, data = df_mono_acin)

# analyzing the taxonomy of interactions

