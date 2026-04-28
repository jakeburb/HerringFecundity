rm(list = ls())
library(gulf)
library(FSA)
library(ggplot2)
library(car)
library(dplyr)
library(gridExtra)
library(tidyverse)
library(stringi)
cat("\f")
clg()

current_path = rstudioapi::getActiveDocumentContext()$path
setwd(dirname(current_path ))

fp <- getwd()


##Load Fish details csv
##NOTE ENSURE DATE COLOUMN IS IN NUMBER FORMAT FOR ANALYSIS TO RUN PROPERLY. Will convert date from number format to date format R recognizes later in code
HDet <- read.csv("herring_fish_detail_report_2022.csv")
str(HDet)
##Subset to include only individuals with egg counts
HDet<-subset(HDet, gonad_count != "NA")
str(HDet)

##Compute fecundity from subsample egg count, subsample gonad weight and total gonad weight
Fec<- (HDet$gonad_count*HDet$gonad_weight)/HDet$gonad_sub_sample_weight
Fec
##Add Fecundity to dataframe
HDet<-data.frame(HDet, Fec)
str(HDet)
excel(HDet)

##Wucik plot to spot outliers
plot(HDet$fish_length, HDet$Fec)
plot(log(HDet$fish_length), log(HDet$Fec))
##Identify and remove outlier/measurement errors
HDet[which.min(HDet$Fec),]
HDet<-subset(HDet, gonad_count != 23)
HDet<-subset(HDet, gonad_count != 33)
HDet<-subset(HDet, gonad_photo_id != 13190)
str(HDet)
##Put Date from Number format into date R can reconize
HDet$sample_date<- as.Date(HDet$sample_date, origin = "1899-12-30")

##Filter for Fall and Spring Spwaners
#Using Table 7 from McQuinn to assign spring or fall spwaner based on maturity stage and capture date
HDetS<-HDet %>% filter(sample_date < '2022-07-01') #Spring
str(HDetS)

HDetF<-HDet %>% filter(sample_date > '2022-07-01') #Fall, Can specify any fish captured after July 1st is fall. If captured in October and stage 5 it would be spring, 
str(HDetF)  #Check data to ensure any fish captured in Oct are in fact stage 6, if stage 6 fall, if stage 5, must be assigned as spring.
#All are stage 6 so no further action needed.
##mean, sd, max and min of fecundity
mean(HDetS$Fec)
sd(HDetS$Fec)
max(HDetS$Fec)
min(HDetS$Fec)

mean(HDetF$Fec)
sd(HDetF$Fec)
max(HDetF$Fec)
min(HDetF$Fec)



  
######################################################################
################Global Length-Fecundity###############################
######################################################################
plot(HDet$fish_length, HDet$Fec)
plot(log(HDet$fish_length), log(HDet$Fec))
LFreg<-lm(log(HDet$Fec)~log(HDet$fish_length))
summary(LFreg)
bG<-LFreg$coefficients[2]
aG<-exp(LFreg$coefficients[1])
CIG<-data.frame(confint(LFreg))
bULG<-CIG$X97.5..[2]
aULG<-exp(CIG$X97.5..[1])
bLLG<-CIG$X2.5..[2]
aLLG<-exp(CIG$X2.5..[1])

lengths<-seq(225,330,0.2)
predG<-aG*lengths^bG
LCI<-aLLG*lengths^bULG
UCI<-aULG*lengths^bLLG
PREDIC<-data.frame(lengths,predG,LCI, UCI)
##Messeih (1976) equation and data
#Generate equation 
lengthsM<-seq(255,380,0.2)
MESSEIH<-(exp(-24.945)*1000)*lengthsM^5.127
PREDICM<-data.frame(lengthsM,MESSEIH)
##Load Messeih data
MDet<-read.csv("herringfecundity.csv")
names(MDet) = c("day", "month", "year", "sampleID", "fishnumber", "length", "weight",  "stage", "age",  "gonad weight", "Fec")
str(MDet)
MDet$year<-as.numeric(MDet$year)
MDet<-subset(MDet, year == 1971)
str(MDet)
MDet<-subset(MDet, length != "NA")
str(MDet)

#Spring Messeih data
MDetS<-subset(MDet, month < 7)
str(MDetS)

#Fall Messieh data
MDetF<-subset(MDet, month >= 7)
str(MDetF)



##Nice Plot True Values
fecplotG<-ggplot() + 
  geom_point(data=MDet, aes(x=length, y=Fec), size=2, alpha=0.3, col= "blue")+
  geom_point(data=HDet, aes(x=fish_length, y=Fec), size=2, col= "black")+
  geom_ribbon(data=PREDIC,aes(x=lengths,ymin=LCI,ymax=UCI),alpha=0.25) +
  geom_line(data=PREDIC,aes(x=lengths,y=predG),size=2) +
  geom_line(data=PREDICM,aes(x=lengthsM,y=MESSEIH),size=1.5, linetype = "dashed", col = "blue") +
  ylab("Fecundity") +
  xlab("Fish Length (mm)")+
  scale_x_continuous(limits = c(230,380), breaks = seq(230, 380, by = 30))+
  scale_y_continuous(limits = c(5000,290000), breaks = seq(10000, 290000, by = 40000))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        axis.text=element_text(size=20, colour="black"),
        axis.title=element_text(size=20))
fecplotG
ggsave("Global Length Fecundity Plot 2022.png", dpi=600, width=22, height= 20, units= "cm")

###Nice Plot ln tranformed

LNfecplotG<-ggplot(HDet, aes(x=log(fish_length), y=log(Fec))) + 
  geom_point(size=2, col= "black")+
  geom_smooth(method=lm,se=TRUE,fullrange=TRUE, col="black", size=2)+
  ylab("ln(Fecundity)") +
  xlab("ln(Fish Length (mm))")+
  xlim(5.45,5.8)+
  scale_y_continuous(limits = c(8.75,11.75), breaks = seq(9, 11.5, by = 0.5))+  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        axis.text=element_text(size=20, colour="black"),
        axis.title=element_text(size=20))
LNfecplotG
ggsave("Global Length Fecundity Plot 2022 LN Transformed.png", dpi=600, width=22, height= 20, units= "cm")

######################################################################
################Global Weight-Fecundity###############################
######################################################################
plot(HDet$fish_weight, HDet$Fec)
plot(log(HDet$fish_weight), log(HDet$Fec))
WFreg<-lm(log(HDet$Fec)~log(HDet$fish_weight))
summary(WFreg)
bWG<-WFreg$coefficients[2]
aWG<-exp(WFreg$coefficients[1])
CIWG<-data.frame(confint(WFreg))
bULWG<-CIWG$X97.5..[2]
aULWG<-exp(CIWG$X97.5..[1])
bLLWG<-CIWG$X2.5..[2]
aLLWG<-exp(CIWG$X2.5..[1])

weight<-seq(90,330,0.2)
predWG<-aWG*weight^bWG
LCIW<-aLLWG*weight^bULWG
UCIW<-aULWG*weight^bLLWG
PREDICW<-data.frame(weight,predWG,LCIW, UCIW)

##COmpute estimates for equation in Messeih (1976) and add raw points
#find max and min weight first
min(MDet$weight)
max(MDet$weight)
weightM<-seq(135,495, 0.2)
MESSEIHW<-(exp(-4.181)*1000)*weightM^1.565
PREDICMW<-data.frame(weightM,MESSEIHW)
##Nice Plot True Values
fecplotGW<-ggplot() + 
  geom_point(data=MDet, aes(x=weight, y=Fec), size=2, alpha=0.3, col= "blue")+
  geom_point(data=HDet, aes(x=fish_weight, y=Fec), size=2, col= "black")+
  geom_ribbon(data=PREDICW,aes(x=weight,ymin=LCIW,ymax=UCIW),alpha=0.25) +
  geom_line(data=PREDICW,aes(x=weight,y=predWG),size=2) +
  geom_line(data=PREDICMW,aes(x=weightM,y=MESSEIHW),size=1.5, linetype = "dashed", col = "blue") +
  ylab("Fecundity") +
  xlab("Fish Weight (g)")+
  scale_x_continuous(limits = c(90,495), breaks = seq(100, 490, by = 50))+
  scale_y_continuous(limits = c(5000,260000), breaks = seq(10000, 260000, by = 40000))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        axis.text=element_text(size=20, colour="black"),
        axis.title=element_text(size=20))
fecplotGW
ggsave("Global Weight Fecundity Plot 2022.png", dpi=600, width=22, height= 20, units= "cm")

###Nice Plot ln tranformed

LNfecplotGW<-ggplot(HDet, aes(x=log(fish_weight), y=log(Fec))) + 
  geom_point(size=2, col= "black")+
  geom_smooth(method=lm,se=TRUE,fullrange=TRUE, col="black", size=2)+
  ylab("ln(Fecundity)") +
  xlab("ln(Fish Weight (g))")+
  xlim(4.5,5.8)+
  scale_y_continuous(limits = c(8.75,11.75), breaks = seq(9, 11.5, by = 0.5))+  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        axis.text=element_text(size=20, colour="black"),
        axis.title=element_text(size=20))
LNfecplotGW
ggsave("Global Weight Fecundity Plot 2022 LN Transformed.png", dpi=600, width=22, height= 20, units= "cm")

##Global Length-fecundity and weight-fecundity side by side
ggsave("Global Length and Weight Fecundity Plot 2022.png", arrangeGrob(fecplotG,fecplotGW, nrow=2), dpi=600, width=22, height=40, unit= "cm")

#####################################################
###############Spring Length Fecundity###############
####################################################
plot(HDetS$fish_length, HDetS$Fec)
plot(log(HDetS$fish_length), log(HDetS$Fec))
LFregS<-lm(log(HDetS$Fec)~log(HDetS$fish_length))
summary(LFregS)
bS<-LFregS$coefficients[2]
aS<-exp(LFregS$coefficients[1])
CIS<-data.frame(confint(LFregS))
bULS<-CIS$X97.5..[2]
aULS<-exp(CIS$X97.5..[1])
bLLS<-CIS$X2.5..[2]
aLLS<-exp(CIS$X2.5..[1])

lengths<-seq(225,330,0.2)
predS<-aS*lengths^bS
LCIS<-aLLS*lengths^bULS
UCIS<-aULS*lengths^bLLS
PREDICS<-data.frame(lengths,predS,LCIS, UCIS)

#Messieh (1976) spring equation and predictions
min(MDetS$length)
max(MDetS$length)
max(MDetS$Fec)
lengthsSM<-seq(255,350,0.2)
MESSEIHS<-(exp(-17.30)*1000)*lengthsSM^3.75
PREDICSM<-data.frame(lengthsSM,MESSEIHS)

##Nice Plot True Values
fecplotS<-ggplot() + 
  geom_point(data=MDetS, aes(x=length, y=Fec), size=2, alpha=0.3, col= "blue")+
  geom_point(data=HDetS, aes(x=fish_length, y=Fec), size=2, col= "black")+
  geom_ribbon(data=PREDICS,aes(x=lengths,ymin=LCIS,ymax=UCIS),alpha=0.25) +
  geom_line(data=PREDICS,aes(x=lengths,y=predS),size=2) +
  geom_line(data=PREDICSM,aes(x=lengthsSM,y=MESSEIHS),size=1.5, linetype = "dashed", col = "blue") +
  ylab("Fecundity") +
  xlab("Fish Length (mm)")+
  xlim(230,350)+
  scale_y_continuous(limits = c(5000,110000), breaks = seq(10000, 110000, by = 20000))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        axis.text=element_text(size=20, colour="black"),
        axis.title=element_text(size=20))
fecplotS
ggsave("Spring Length Fecundity Plot 2022.png", dpi=600, width=22, height= 20, units= "cm")

###Nice Plot ln tranformed

LNfecplotS<-ggplot(HDetS, aes(x=log(fish_length), y=log(Fec))) + 
  geom_point(size=2, col= "black")+
  geom_smooth(method=lm,se=TRUE,fullrange=TRUE, col="black", size=2)+
  ylab("ln(Fecundity)") +
  xlab("ln(Fish Length (mm))")+
  xlim(5.45,5.8)+
  scale_y_continuous(limits = c(8.75,11.75), breaks = seq(9, 11.5, by = 0.5))+  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        axis.text=element_text(size=20, colour="black"),
        axis.title=element_text(size=20))
LNfecplotS
ggsave("Spring Length Fecundity Plot 2022 LN Transformed.png", dpi=600, width=22, height= 20, units= "cm")

######################################################################
################Spring Weight-Fecundity###############################
######################################################################
plot(HDetS$fish_weight, HDetS$Fec)
plot(log(HDetS$fish_weight), log(HDetS$Fec))
WFregS<-lm(log(HDetS$Fec)~log(HDetS$fish_weight))
summary(WFregS)
bWS<-WFregS$coefficients[2]
aWS<-exp(WFregS$coefficients[1])
CIWS<-data.frame(confint(WFregS))
bULWS<-CIWS$X97.5..[2]
aULWS<-exp(CIWS$X97.5..[1])
bLLWS<-CIWS$X2.5..[2]
aLLWS<-exp(CIWS$X2.5..[1])

weight<-seq(90,330,0.2)
predWS<-aWS*weight^bWS
LCIWS<-aLLWS*weight^bULWS
UCIWS<-aULWS*weight^bLLWS
PREDICWS<-data.frame(weight,predWS,LCIWS, UCIWS)

#Messieh (1976) spring equation and predictions
min(MDetS$weight)
max(MDetS$weight)
max(MDetS$Fec)
weightSM<-seq(135, 370,0.2)
MESSEIHWSW<-(exp(-2.47)*1000)*weightSM^1.22
PREDICSMW<-data.frame(weightSM,MESSEIHWSW)

##Nice Plot True Values
  fecplotWS<-ggplot() + 
  geom_point(data=MDetS, aes(x=weight, y=Fec), size=2, alpha=0.3, col= "blue")+
  geom_point(data=HDetS, aes(x=fish_weight, y=Fec), size=2, col= "black")+
  geom_ribbon(data=PREDICWS,aes(x=weight,ymin=LCIWS,ymax=UCIWS),alpha=0.25) +
  geom_line(data=PREDICWS,aes(x=weight,y=predWS),size=2) +
  geom_line(data=PREDICSMW,aes(x=weightSM,y=MESSEIHWSW),size=1.5, linetype = "dashed", col = "blue") +
  ylab("Fecundity") +
  xlab("Fish Weight (g)")+
  scale_x_continuous(limits = c(90,370), breaks = seq(100, 370, by = 40))+
  scale_y_continuous(limits = c(5000,110000), breaks = seq(10000, 110000, by = 20000))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        axis.text=element_text(size=20, colour="black"),
        axis.title=element_text(size=20))
fecplotWS
ggsave("Spring Weight Fecundity Plot 2022.png", dpi=600, width=22, height= 20, units= "cm")

###Nice Plot ln tranformed

LNfecplotWS<-ggplot(HDetS, aes(x=log(fish_weight), y=log(Fec))) + 
  geom_point(size=2, col= "black")+
  geom_smooth(method=lm,se=TRUE,fullrange=TRUE, col="black", size=2)+
  ylab("ln(Fecundity)") +
  xlab("ln(Fish Weight (g))")+
  xlim(4.5,5.8)+
  scale_y_continuous(limits = c(8.75,11.75), breaks = seq(9, 11.5, by = 0.5))+  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        axis.text=element_text(size=20, colour="black"),
        axis.title=element_text(size=20))
LNfecplotWS
ggsave("Spring Weight Fecundity Plot 2022 LN Transformed.png", dpi=600, width=22, height= 20, units= "cm")

##Spring Length-fecundity and weight-fecundity side by side
ggsave("Spring Length and Weight Fecundity Plot 2022.png", arrangeGrob(fecplotS,fecplotWS, nrow=2), dpi=600, width=22, height=40, unit= "cm")

#########################################################
####################Fall Length Fecundity################
#########################################################
plot(HDetF$fish_length, HDetF$Fec)
plot(log(HDetF$fish_length), log(HDetF$Fec))
LFregF<-lm(log(HDetF$Fec)~log(HDetF$fish_length))
summary(LFregF)
bF<-LFregF$coefficients[2]
aF<-exp(LFregF$coefficients[1])
CIF<-data.frame(confint(LFregF))
bULF<-CIF$X97.5..[2]
aULF<-exp(CIF$X97.5..[1])
bLLF<-CIF$X2.5..[2]
aLLF<-exp(CIF$X2.5..[1])

lengths<-seq(225,330,0.2)
predF<-aF*lengths^bF
LCIF<-aLLF*lengths^bULF
UCIF<-aULF*lengths^bLLF
PREDICF<-data.frame(lengths,predF,LCIF, UCIF)

###Create predictions for Messeih (1976) equation across the length range observed for fall herring 
max(MDetF$length)
min(MDetF$length)
max(MDetF$Fec)
lengthsFM<-seq(280,380,0.2)
MESSEIHF<-(exp(-17.87)*1000)*lengthsFM^3.93
PREDICFM<-data.frame(lengthsFM,MESSEIHF)
##Nice Plot True Values
fecplotF<-ggplot() + 
  geom_point(data=MDetF, aes(x=length, y=Fec), size=2, alpha=0.3, col= "blue")+
  geom_point(data=HDetF, aes(x=fish_length, y=Fec), size=2, col= "black")+
  geom_ribbon(data=PREDICF,aes(x=lengths,ymin=LCIF,ymax=UCIF),alpha=0.25) +
  geom_line(data=PREDICF,aes(x=lengths,y=predF),size=2) +
  geom_line(data=PREDICFM,aes(x=lengthsFM,y=MESSEIHF),size=1.5, linetype = "dashed", col = "blue") +
  ylab("Fecundity") +
  xlab("Fish Length (mm)")+
  scale_x_continuous(limits = c(230,380), breaks = seq(230, 380, by = 30))+
  scale_y_continuous(limits = c(5000,260000), breaks = seq(10000, 260000, by = 40000))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        axis.text=element_text(size=20, colour="black"),
        axis.title=element_text(size=20))
fecplotF
ggsave("Fall Length Fecundity Plot 2022.png", dpi=600, width=22, height= 20, units= "cm")

###Nice Plot ln tranformed

LNfecplotF<-ggplot(HDetF, aes(x=log(fish_length), y=log(Fec))) + 
  geom_point(size=2, col= "black")+
  geom_smooth(method=lm,se=TRUE,fullrange=TRUE, col="black", size=2)+
  ylab("ln(Fecundity)") +
  xlab("ln(Fish Length (mm))")+
  xlim(5.45,5.8)+
  scale_y_continuous(limits = c(8.75,11.75), breaks = seq(9, 11.5, by = 0.5))+  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        axis.text=element_text(size=20, colour="black"),
        axis.title=element_text(size=20))
LNfecplotF
ggsave("Fall Length Fecundity Plot 2022 LN Transformed.png", dpi=600, width=22, height= 20, units= "cm")

######################################################################
################Fall Weight-Fecundity###############################
######################################################################
plot(HDetF$fish_weight, HDetF$Fec)
plot(log(HDetF$fish_weight), log(HDetF$Fec))
WFregF<-lm(log(HDetF$Fec)~log(HDetF$fish_weight))
summary(WFregF)
bWF<-WFregF$coefficients[2]
aWF<-exp(WFregF$coefficients[1])
CIWF<-data.frame(confint(WFregF))
bULWF<-CIWF$X97.5..[2]
aULWF<-exp(CIWF$X97.5..[1])
bLLWF<-CIWF$X2.5..[2]
aLLWF<-exp(CIWF$X2.5..[1])

weight<-seq(90,330,0.2)
predWF<-aWF*weight^bWF
LCIWF<-aLLWF*weight^bULWF
UCIWF<-aULWF*weight^bLLWF
PREDICWF<-data.frame(weight,predWF,LCIWF, UCIWF)
###Create predictions for Messeih (1976) equation across the length range observed for fall herring 
min(MDet$weight)
max(MDet$weight)
weightFWM<-seq(200,495,0.2)
MESSEIHWF<-(exp(-2.28)*1000)*weightFWM^1.26
PREDICWFM<-data.frame(weightFWM, MESSEIHWF)
##Nice Plot True Values
fecplotWF<-ggplot() + 
  geom_point(data=MDetF, aes(x=weight, y=Fec), size=2, alpha=0.3, col= "blue")+
  geom_point(data=HDetF, aes(x=fish_weight, y=Fec), size=2, col= "black")+
  geom_ribbon(data=PREDICWF,aes(x=weight,ymin=LCIWF,ymax=UCIWF),alpha=0.25) +
  geom_line(data=PREDICWF,aes(x=weight,y=predWF),size=2) +
  geom_line(data=PREDICWFM,aes(x=weightFWM,y=MESSEIHWF),size=1.5, linetype = "dashed", col = "blue") +
  ylab("Fecundity") +
  xlab("Fish Weight (g)")+
  xlim(90,495)+
  scale_y_continuous(limits = c(5000,260000), breaks = seq(10000, 260000, by = 40000))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        axis.text=element_text(size=20, colour="black"),
        axis.title=element_text(size=20))
fecplotWF
ggsave("Fall Weight Fecundity Plot 2022.png", dpi=600, width=22, height= 20, units= "cm")

###Nice Plot ln tranformed

LNfecplotWF<-ggplot(HDetF, aes(x=log(fish_weight), y=log(Fec))) + 
  geom_point(size=2, col= "black")+
  geom_smooth(method=lm,se=TRUE,fullrange=TRUE, col="black", size=2)+
  ylab("ln(Fecundity)") +
  xlab("ln(Fish Weight (g))")+
  xlim(4.5,5.8)+
  scale_y_continuous(limits = c(8.75,11.75), breaks = seq(9, 11.5, by = 0.5))+  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        axis.text=element_text(size=20, colour="black"),
        axis.title=element_text(size=20))
LNfecplotWF
ggsave("Fall Weight Fecundity Plot 2022 LN Transformed.png", dpi=600, width=22, height= 20, units= "cm")

##Fall Length-fecundity and weight-fecundity side by side
ggsave("Fall Length and Weight Fecundity Plot 2022.png", arrangeGrob(fecplotF,fecplotWF, nrow=2), dpi=600, width=22, height=40, unit= "cm")


########################Fall and Spring on one Plot####################
#########################Length############################
##Nice Plot True Values
fecplotSF<-ggplot() + 
  geom_point(data=HDetS, aes(x=fish_length, y=Fec), size=2, col= "green4")+
  geom_ribbon(data=PREDICS,aes(x=lengths,ymin=LCIS,ymax=UCIS), fill= "green4",alpha=0.25) +
  geom_line(data=PREDICS,aes(x=lengths,y=predS),size=1.5, col= "green4") +
  geom_point(data=HDetF, aes(x=fish_length, y=Fec), size=2, col= "orange")+
  geom_ribbon(data=PREDICF,aes(x=lengths,ymin=LCIF,ymax=UCIF),alpha=0.25, fill = "orange") +
  geom_line(data=PREDICF,aes(x=lengths,y=predF),size=1.5, col = "orange") +
  ylab("Fecundity") +
  xlab("Fish Length (mm)")+
  xlim(230,330)+
  scale_y_continuous(limits = c(5000,130000), breaks = seq(10000, 130000, by = 20000))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        axis.text=element_text(size=20, colour="black"),
        axis.title=element_text(size=20))
fecplotSF
ggsave("Spring and Fall Length Fecundity Plot 2022.png", dpi=600, width=22, height= 20, units= "cm")


###Spring and Fall average fecundity boxplots and statistical test
str(HDetS)
SC<- replicate(86, "Spring")
HDetS<-cbind(HDetS,SC)
str(HDetS)

str(HDetF)
SC<- replicate(84, "Fall")
HDetF<-cbind(HDetF,SC)
str(HDetF)

HDetN<- rbind(HDetS, HDetF)
str(HDetN)

ggplot(data=HDetN, aes(x=SC, y=Fec, colour = SC))+
  geom_boxplot(alpha=0, size=1.25)+
  geom_jitter(alpha=0.25)+
  ylab("Fecundity")+
  xlab("Population")+
  scale_y_continuous(limits = c(5000,120000), breaks = seq(10000, 120000, by = 20000))+
  scale_x_discrete(limits = c("Spring", "Fall"))+
  scale_colour_manual(breaks = c("Spring", "Fall"), values = c("green4", "orange"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        legend.position = "none",
        axis.text=element_text(size=14, colour="black"),
        axis.title=element_text(size=20))
ggsave("BoxPlot Spring and Fall Length Fecundity 2022.png", dpi=600, width=12, height= 20, units= "cm")

#Check if fecundity of fall and spring spawners is significantly different in 2022
SCtest<-t.test(HDetN$Fec~HDetN$SC)
SCtest
shapiro.test(HDetN$Fec)##Check model assumptions


##Spring and Fall Herring Average Fecundity 70s, 80s, 2022 by season and year
YDet<-read.csv("FecundityByYear.csv")
str(YDet)
YDet$SC<-factor(YDet$SC, levels = c("Spring", "Fall"))

ggplot(data=YDet, aes(x=Year, y=Fec, colour= SC))+
  geom_boxplot(alpha=0, size=1.25)+
  ylab("Fecundity")+
  xlab("Year")+
  scale_y_continuous(limits = c(5000,270000), breaks = seq(10000, 270000, by = 20000))+
  scale_colour_manual(breaks = c("Spring", "Fall"), values = c("green4", "orange"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        legend.title = element_blank(),
        legend.text = element_text(size=14, colour="black"),
        axis.text=element_text(size=14, colour="black"),
        axis.title=element_text(size=20))
ggsave("BoxPlot Spring and Fall Length Fecundity By Year.png", dpi=600, width=16, height= 22, units= "cm")

anovatest<-aov(YDet$Fec~ YDet$Year + YDet$SC + YDet$Year*YDet$SC)
summary(anovatest)
TukeyHSD(anovatest)

FecS<-subset(YDet, SC == "Spring")
FecF<-subset(YDet, SC =="Fall")

Stest<-aov(FecS$Fec~FecS$Year)
summary(Stest)
TukeyHSD(Stest)

Ftest<-aov(FecF$Fec~FecF$Year)
summary(Ftest)
TukeyHSD(Ftest)

##Compute mean, SD, max and min of fecundity across 1970s, 80s and 2022
Fec1970<-subset(YDet, Year == "1970s" )
Fec1970S<-subset(Fec1970, SC == "Spring")
Fec1970F<- subset(Fec1970, SC == "Fall")

Fec1980<-subset(YDet, Year == "1980s" )
Fec1980S<-subset(Fec1980, SC == "Spring")
Fec1980F<- subset(Fec1980, SC == "Fall")

#Spring 1970s
mean(Fec1970S$Fec)
sd(Fec1970S$Fec)
max(Fec1970S$Fec)
min(Fec1970S$Fec)

#Fall 1970s
mean(Fec1970F$Fec)
sd(Fec1970F$Fec)
max(Fec1970F$Fec)
min(Fec1970F$Fec)


#Spring 1980s
mean(Fec1980S$Fec)
sd(Fec1980S$Fec)
max(Fec1980S$Fec)
min(Fec1980S$Fec)

#Fall 1980s
mean(Fec1980F$Fec)
sd(Fec1980F$Fec)
max(Fec1980F$Fec)
min(Fec1980F$Fec)

#Spring 2022 
mean(HDetS$Fec)
sd(HDetS$Fec)
max(HDetS$Fec)
min(HDetS$Fec)

#Fall 2022
mean(HDetF$Fec)
sd(HDetF$Fec)
max(HDetF$Fec)
min(HDetF$Fec)


#####Estimate potential reproductive output using hypothetical situation##
###First just try 1million individuals with average fecundity in 1970s, 1980s and 2022
#####Spring#########
n<-1000000
###1970s
MH1970S<-n*(mean(Fec1970S$Fec))
MH1970S
MH1970ULS<-MH1970S + (1.96*sd(Fec1970S$Fec))
MH1970ULS
MH1970LLS<-MH1970S - (1.96*sd(Fec1970S$Fec))
MH1970LLS
###1980s
MH1980S<-n*(mean(Fec1980S$Fec))
MH1980S
MH1980ULS<-MH1980S + (1.96*sd(Fec1980S$Fec))
MH1980ULS
MH1980LLS<-MH1980S - (1.96*sd(Fec1980S$Fec))
MH1980LLS
##2022
MH2022S<-n*(mean(HDetS$Fec))
MH2022S
MH2022ULS<-MH2022S + (1.96*sd(HDetS$Fec))
MH2022ULS
MH2022LLS<-MH2022S - (1.96*sd(HDetS$Fec))
MH2022LLS


####Fall######
###1970s
MH1970F<-n*(mean(Fec1970F$Fec))
MH1970F
MH1970ULF<-MH1970F + (1.96*sd(Fec1970F$Fec))
MH1970ULF
MH1970LLF<-MH1970F - (1.96*sd(Fec1970F$Fec))
MH1970LLF
###1980s
MH1980F<-n*(mean(Fec1980F$Fec))
MH1980F
MH1980ULF<-MH1980F + (1.96*sd(Fec1980F$Fec))
MH1980ULF
MH1980LLF<-MH1980F - (1.96*sd(Fec1980F$Fec))
MH1980LLF
##2022
MH2022F<-n*(mean(HDetF$Fec))
MH2022F
MH2022ULF<-MH2022F + (1.96*sd(HDetF$Fec))
MH2022ULF
MH2022LLF<-MH2022F - (1.96*sd(HDetF$Fec))
MH2022LLF

M<-(c(MH1970S,MH1980S,MH2022S,MH1970F,MH1980F,MH2022F))/100000000
LL<-c(MH1970LLS,MH1980LLS,MH2022LLS,MH1970LLF,MH1980LLF,MH2022LLF)/100000000
UL<-c(MH1970ULS,MH1980ULS,MH2022ULS,MH1970ULF,MH1980ULF,MH2022ULF)/100000000
years<-c("1970s", "1980s", "2022", "1970s", "1980s", "2022")
spawningC<-c("Spring", "Spring","Spring", "Fall", "Fall", "Fall")

HypRO<-data.frame(years,spawningC,M,LL,UL)
HypRO

##Simulation of potential reproductive output from 1 million individuals for fall and spring 2022 and fall and spring 1970s based on Messaih 1976
###1 million individuls ranging from 230 - 325mm for 2022 spring and fall
###1 million individuals ranging from 260-350 for spring and 260-380 for fall for 1970s
n2022O<-sample(300:325, 200000, replace =TRUE)
n2022y<-sample(230:299, 800000, replace =TRUE)
n2022<-c(n2022y, n2022O)

n1970O<-sample(320:350, 200000, replace =TRUE)
n1970y<-sample(260:319, 800000, replace =TRUE)
n1970S<-c(n1970y,n1970O)
n1970FO<-sample(350:380, 200000, replace =TRUE)
n1970Fy<-sample(260:349, 800000, replace =TRUE)
n1970F<-c(n1970Fy,n1970FO)

##Spring 2022 simulation
fec2022S<-aS*n2022^bS
Mfec2022S<-sum(fec2022S)
Mfec2022LLS<-sum(fec2022S) - 1.96*(sd(fec2022S))
Mfec2022ULS<-sum(fec2022S) + 1.96*(sd(fec2022S))

#Spring 1976 simulation
fecMesS<-(exp(-17.30)*1000)*n1970S^3.75
MfecMesS<-sum(fecMesS)
MfecMesLLS<-sum(fecMesS) - 1.96*(sd(fecMesS))
MfecMesULS<-sum(fecMesS) + 1.96*(sd(fecMesS))

##Fall 2022 simulation
fec2022F<-aF*n2022^bF
Mfec2022F<-sum(fec2022F)
Mfec2022LLF<-sum(fec2022F) - 1.96*(sd(fec2022F))
Mfec2022ULF<-sum(fec2022F) + 1.96*(sd(fec2022F))

##Fall 1976 simulation
fecMesF<-(exp(-17.87)*1000)*n1970F^3.93
MfecMesF<-sum(fecMesF)
MfecMesLLF<-sum(fecMesF) - 1.96*(sd(fecMesF))
MfecMesULF<-sum(fecMesF) + 1.96*(sd(fecMesF))


ML<-(c(MfecMesS,Mfec2022S,MfecMesF,Mfec2022F))/1000000000
LLl<-c(MfecMesLLS,Mfec2022LLS,MfecMesLLF, Mfec2022LLF)/1000000000
ULl<-c(MfecMesULS,Mfec2022ULS,MfecMesULF,Mfec2022ULF)/1000000000
yearsl<-c("1976", "2022", "1976", "2022")
spawningCl<-c("Spring", "Spring", "Fall", "Fall")

SimF<-data.frame(yearsl,spawningCl,ML,LLl,ULl)
SimF


ggplot(data=SimF, aes(x=yearsl, y=ML, colour= spawningCl))+
  geom_point(size=5)+
  ylab("Potential Reproductive Output (Billions)")+
  xlab("Year")+
  scale_colour_manual(breaks = c("Spring", "Fall"), values = c("green4", "orange"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        legend.title = element_blank(),
        legend.text = element_text(size=14, colour="black"),
        axis.text=element_text(size=14, colour="black"),
        axis.title=element_text(size=20))
ggsave("Simulated Reproductive Output.png", dpi=600, width=16, height= 22, units= "cm")


############Do this 1000 times using replicate#######
nu<-1000
#####Spring 2022##########
S2022<-replicate(nu, {n2022O<-sample(300:325, 200000, replace =TRUE) 
                      n2022y<-sample(230:299, 800000, replace =TRUE)
                      n2022<-c(n2022y, n2022O)
                      fec2022S<-aS*n2022^bS
                      Mfec2022S<-sum(fec2022S)
                      Mfec2022S})
#Take mean of Spring 2022 and compute confidence limits
hist(S2022)
MeanS2022<-mean(S2022)/1000000000
stddevS2022<-sd(S2022)
LLS2022<-(mean(S2022)- 1.96*(sd(S2022)))/1000000000
ULS2022<-(mean(S2022) + 1.96*(sd(S2022)))/1000000000
medS2022<-median(S2022)/1000000000

####Spring 1970###### 
S1970<-replicate(nu, {n1970O<-sample(320:350, 200000, replace =TRUE)
                      n1970y<-sample(260:319, 800000, replace =TRUE)
                      n1970S<-c(n1970y,n1970O)
                      fecMesS<-(exp(-17.30)*1000)*n1970S^3.75
                      MfecMesS<-sum(fecMesS)
                      MfecMesS})
                      
#Take mean of Spring 1970 and compute confidence limits
hist(S1970)
MeanS1970<-mean(S1970)/1000000000
stddevS1970<-sd(S1970)
LLS1970<-(mean(S1970) - 1.96*(sd(S1970)))/1000000000
ULS1970<-(mean(S1970) + 1.96*(sd(S1970)))/1000000000
medS1970<-median(S1970)/1000000000


#######Fall 2022##########
F2022<-replicate(nu, {n2022O<-sample(300:325, 200000, replace =TRUE) 
                      n2022y<-sample(230:299, 800000, replace =TRUE)
                      n2022<-c(n2022y, n2022O)
                      fec2022F<-aF*n2022^bF
                      Mfec2022F<-sum(fec2022F)
                      Mfec2022F})
#Take mean, median and sd of Fall 2022 and compute confidence limits
hist(F2022)
MeanF2022<-mean(F2022)/1000000000
stddevF2022<-sd(S2022)
LLF2022<-(mean(F2022) - 1.96*(sd(S2022)))/1000000000
ULF2022<-(mean(F2022) + 1.96*(sd(S2022)))/1000000000
medF2022<-median(F2022)/1000000000


######Fall 1970###################
F1970<-replicate(nu, {n1970FO<-sample(350:380, 200000, replace =TRUE)
                      n1970Fy<-sample(260:349, 800000, replace =TRUE)
                      n1970F<-c(n1970Fy,n1970FO)
                      fecMesF<-(exp(-17.87)*1000)*n1970F^3.93
                      MfecMesF<-sum(fecMesF)
                      MfecMesF})
hist(F1970)
MeanF1970<-mean(F1970)/1000000000
stddevF1970<-sd(F1970)
LLF1970<-(mean(F1970) - 1.96*(sd(F1970)))/1000000000
ULF1970<-(mean(F1970) + 1.96*(sd(F1970)))/1000000000
medf1970<-median(F1970)/1000000000


##Create dataframe
ML<-c(medS1970,medS2022,medf1970, medF2022)
LLl<-c(LLS1970,LLS2022,LLF1970, LLF2022)
ULl<-c(ULS1970,ULS2022,ULF1970,Mfec2022ULF)
yearsl<-c("1976", "2022", "1976", "2022")
spawningCl<-c("Spring", "Spring", "Fall", "Fall")

SimFT<-data.frame(yearsl,spawningCl,ML,LLl,ULl)
SimFT


###Examining Relative Fecundity (standerdized by length) as per Tim Wards SUggestions###
fecplotF<-ggplot() + 
  geom_point(data=MDetF, aes(x=length, y=(Fec/weight), size=2, alpha=0.3, col= "blue"))+
  geom_point(data=HDetF, aes(x=fish_length, y=(Fec/fish_weight), size=2, col= "black"))+
  ylab("Fecundity") +
  xlab("Fish Length (mm)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        axis.text=element_text(size=20, colour="black"),
        axis.title=element_text(size=20))
fecplotF

ggsave("Fall Length RELATIVE Fecundity Plot 2022_TEST.png", dpi=600, width=22, height= 20, units= "cm")