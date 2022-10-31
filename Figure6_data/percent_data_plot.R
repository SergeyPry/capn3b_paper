library(ggplot2)
library(dplyr)
library(ggsci)

#setwd("c:/Dropbox/capn3b-dmd_paper/Figure6_data")

# Sample	Gene	Expression
data <- read.csv("percent_MC_data.txt", sep = "\t", header = TRUE)
str(data)

data$Genotype <- factor(data$Genotype, levels = c("wildtype", "capn3b rnaless")) 
data$Treatment <-factor(data$Treatment, levels = c("control","0.3uM APM", "0.5uM APM", "0.75uM APM"))


pd = position_dodge(0.8)

# plot trial
ggplot(data, aes(x=Treatment, y= Percent_abnormal, fill = Genotype)) +
  geom_bar(stat = "summary", position="dodge", width = 0.8, alpha = 0.7) +
  geom_point(aes(x= Treatment, y= Percent_abnormal, fill = Genotype), size=2, color = "black", shape = 21, position=position_jitterdodge(0.2)) +
  stat_summary(fun.min=function(x)(mean(x)-sd(x)/sqrt(length(x))),geom="errorbar", 
               fun.max=function(x)(mean(x)+sd(x)/sqrt(length(x))),width = 0.5, position=pd, size = 1, alpha = 0.7) +
  ylab("Abnormality percentage")+
  scale_x_discrete(labels = c("control","0.3µM APM", "0.5µM APM", "0.75µM APM"))+
  theme_bw() + scale_fill_aaas() + 
  theme(
    axis.title.x = element_text(color="black", size= 18, face="bold", margin=margin(t = 20, unit = "pt")),
    axis.title.y = element_text(color="black", size= 18, face="bold", margin=margin(r = 20, unit = "pt")),
    axis.text.x = element_text(angle = 60, vjust=0.5, colour="grey20", size= 18, face="plain"),
    axis.text.y = element_text(colour="grey20",size= 18, face="plain"),
    legend.title = element_text( size = 16, face = "bold"),
    legend.text = element_text( size = 16, face = "plain"),
    legend.key.size = unit(0.5, "cm")
  )

ggsave("MC_percent_results.png", dpi = 600)



# aggregate data for testing

# https://rcompanion.org/handbook/H_06.html

if(!require(psych)){install.packages("psych")}
if(!require(vcd)){install.packages("vcd")}
if(!require(DescTools)){install.packages("DescTools")}
if(!require(rcompanion)){install.packages("rcompanion")}
library(psych)


data2 <- read.csv("table_MC_data.txt", sep = "\t", header = TRUE)
str(data2)

data2$genotype <- factor(data2$genotype, levels = c("wildtype", "capn3b rnaless")) 
data2$treatment <-factor(data2$treatment, levels = c("control","0.3uM APM", "0.5uM APM", "0.75uM APM"))
data2$phenotype <-factor(data2$phenotype, levels = c("normal","abnormal"))

Table = xtabs(count ~ phenotype + genotype + treatment, data=data2)
Table

ftable(Table)
mant_result <- mantelhaen.test(Table)



library(rcompanion)

groupwiseCMH(Table,
             group   = 3,
             fisher  = TRUE,
             gtest   = FALSE,
             chisq   = TRUE,
             method  = "fdr",
             correct = "none",
             digits  = 3)

# Group       Test  p.value    adj.p
# 1    control Chi.square 4.33e-01 4.33e-01
# 2  0.3uM APM Chi.square 4.02e-05 8.04e-05
# 3  0.5uM APM Chi.square 4.23e-09 1.69e-08
# 4 0.75uM APM Chi.square 4.77e-02 6.36e-02

