library(ggplot2)
library(dplyr)
library(ggsci)

#setwd("c:/Dropbox/capn3b-dmd_paper/Figure4_data")

# Sample	Gene	Expression
data <- read.csv("percent_MC_data.txt", sep = "\t", header = TRUE)
str(data)

data$Genotype <- factor(data$Genotype, levels = c("wildtype", "capn3b mut1",	"capn3b rnaless")) 
data$Treatment <-factor(data$Treatment, levels = c("control","0.8% MC"))


pd = position_dodge(0.8)

# plot trial
ggplot(data, aes(x=Genotype, y= Percent_abnormal, fill = Treatment)) +
  geom_bar(stat = "summary", position="dodge", width = 0.8, alpha = 0.7) +
  geom_point(aes(x= Genotype, color = Treatment), size=3, shape = 21, colour = "black", position=position_jitterdodge(0.2)) +
  stat_summary(fun.min=function(x)(mean(x)-sd(x)/sqrt(length(x))),geom="errorbar", 
               fun.max=function(x)(mean(x)+sd(x)/sqrt(length(x))),width = 0.5, position=pd, size = 1.5, alpha = 0.7) +
  ylab("Abnormality percentage")+
  theme_bw() + scale_fill_aaas() + 
  theme(
    axis.title.x = element_text(color="black", size= 18, face="bold", margin=margin(t = 20, unit = "pt")),
    axis.title.y = element_text(color="black", size= 18, face="bold", margin=margin(r = 20, unit = "pt")),
    axis.text.x = element_text(angle = 60, vjust=0.5, colour="grey20", size= 18, face="plain"),
    axis.text.y = element_text(colour="grey20",size= 18, face="plain"),
    legend.title = element_text( size = 16, face = "bold"),
    legend.text = element_text( size = 16, face = "plain"),
    legend.key.size = unit(1, "cm")
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

data2$genotype <- factor(data2$genotype, levels = c("wildtype", "capn3b mut1",	"capn3b rnaless")) 
data2$treatment <-factor(data2$treatment, levels = c("control","MC"))
data2$phenotype <-factor(data2$phenotype, levels = c("normal","abnormal"))

Table = xtabs(count ~ treatment + phenotype + genotype, data=data2)
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
# 1       wildtype Chi.square 1.01e-02 1.01e-02
# 2    capn3b mut1 Chi.square 4.06e-27 1.22e-26
# 3 capn3b rnaless Chi.square 5.33e-26 8.00e-26

