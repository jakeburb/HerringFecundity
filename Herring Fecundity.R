rm(list = ls())
library(gulf)
library(FSA)
library(ggplot2)
library(car)
library(dplyr)
library(gridExtra)
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
##Wucik plot to spot outliers
plot(HDet$fish_length, HDet$Fec)
plot(log(HDet$fish_length), log(HDet$Fec))
##Identify and remove outlier/measurement errors
HDet[which.min(HDet$Fec),]
HDet<-subset(HDet, gonad_count != 23)
HDet<-subset(HDet, gonad_photo_id != 13190)
str(HDet)
##Put Date from Number format into date R can reconize
HDet$sample_date<- as.Date(HDet$sample_date, origin = "1899-12-30")

##Filter for Fall and Spring Spwaners
#Using Table 7 from McQuinn to assign spring or fall spwaner based on maturity stage and capture date
HDetS<-HDet %>% filter(sample_date < '2022-07-01') #Spring
str(HDetS)

HDetF<-HDet %>% filter(sample_date > '2022-07-01') #Fall, Can specify any fish captured after July 1st is fall. If captured in October and stage 5 it would be spring, 
str(HDetF)  #Check data to ensure any fish captured in Oct are in fact stage 6, if stage 6 fall, if stage 5, must be assigned as spring. #All are stage 6 so no further action needed.
  
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
MESSEIH<-(exp(-24.945)*1000)*lengths^5.127
PREDIC<-data.frame(lengths,predG,LCI, UCI, MESSEIH)
##Nice Plot True Values
fecplotG<-ggplot() + 
  geom_point(data=HDet, aes(x=fish_length, y=Fec), size=2, col= "black")+
  geom_ribbon(data=PREDIC,aes(x=lengths,ymin=LCI,ymax=UCI),alpha=0.25) +
  geom_line(data=PREDIC,aes(x=lengths,y=predG),size=2) +
  geom_line(data=PREDIC,aes(x=lengths,y=MESSEIH),size=1.5, linetype = "dashed", col = "blue") +
  ylab("Fecundity") +
  xlab("Fish Length (mm)")+
  xlim(230,330)+
  scale_y_continuous(limits = c(5000,120000), breaks = seq(10000, 120000, by = 20000))+
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
bG<-WFreg$coefficients[2]
aG<-exp(WFreg$coefficients[1])
CIG<-data.frame(confint(WFreg))
bULG<-CIG$X97.5..[2]
aULG<-exp(CIG$X97.5..[1])
bLLG<-CIG$X2.5..[2]
aLLG<-exp(CIG$X2.5..[1])

weight<-seq(90,330,0.2)
predG<-aG*weight^bG
LCI<-aLLG*weight^bULG
UCI<-aULG*weight^bLLG
MESSEIH<-(exp(-4.181)*1000)*weight^1.565
PREDIC<-data.frame(weight,predG,LCI, UCI, MESSEIH)
##Nice Plot True Values
fecplotGW<-ggplot() + 
  geom_point(data=HDet, aes(x=fish_weight, y=Fec), size=2, col= "black")+
  geom_ribbon(data=PREDIC,aes(x=weight,ymin=LCI,ymax=UCI),alpha=0.25) +
  geom_line(data=PREDIC,aes(x=weight,y=predG),size=2) +
  geom_line(data=PREDIC,aes(x=weight,y=MESSEIH),size=1.5, linetype = "dashed", col = "blue") +
  ylab("Fecundity") +
  xlab("Fish Weight (g)")+
  xlim(90,330)+
  scale_y_continuous(limits = c(5000,120000), breaks = seq(10000, 120000, by = 20000))+
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
MESSEIHS<-(exp(-17.30)*1000)*lengths^3.75
PREDICS<-data.frame(lengths,predS,LCIS, UCIS, MESSEIHS)
##Nice Plot True Values
fecplotS<-ggplot() + 
  geom_point(data=HDetS, aes(x=fish_length, y=Fec), size=2, col= "black")+
  geom_ribbon(data=PREDICS,aes(x=lengths,ymin=LCIS,ymax=UCIS),alpha=0.25) +
  geom_line(data=PREDICS,aes(x=lengths,y=predS),size=2) +
  geom_line(data=PREDICS,aes(x=lengths,y=MESSEIHS),size=1.5, linetype = "dashed", col = "blue") +
  ylab("Fecundity") +
  xlab("Fish Length (mm)")+
  xlim(230,330)+
  scale_y_continuous(limits = c(5000,100000), breaks = seq(10000, 100000, by = 20000))+
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
MESSEIHWS<-(exp(-2.47)*1000)*weight^1.22
PREDICWS<-data.frame(weight,predWS,LCIWS, UCIWS, MESSEIHWS)
##Nice Plot True Values
fecplotWS<-ggplot() + 
  geom_point(data=HDetS, aes(x=fish_weight, y=Fec), size=2, col= "black")+
  geom_ribbon(data=PREDICWS,aes(x=weight,ymin=LCIWS,ymax=UCIWS),alpha=0.25) +
  geom_line(data=PREDICWS,aes(x=weight,y=predWS),size=2) +
  geom_line(data=PREDICWS,aes(x=weight,y=MESSEIHWS),size=1.5, linetype = "dashed", col = "blue") +
  ylab("Fecundity") +
  xlab("Fish Weight (g)")+
  xlim(90,300)+
  scale_y_continuous(limits = c(5000,100000), breaks = seq(10000, 100000, by = 20000))+
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

##Global Length-fecundity and weight-fecundity side by side
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
MESSEIHF<-(exp(-17.87)*1000)*lengths^3.93
PREDICF<-data.frame(lengths,predF,LCIF, UCIF, MESSEIHF)
##Nice Plot True Values
fecplotF<-ggplot() + 
  geom_point(data=HDetF, aes(x=fish_length, y=Fec), size=2, col= "black")+
  geom_ribbon(data=PREDICF,aes(x=lengths,ymin=LCIF,ymax=UCIF),alpha=0.25) +
  geom_line(data=PREDICF,aes(x=lengths,y=predF),size=2) +
  geom_line(data=PREDICF,aes(x=lengths,y=MESSEIHF),size=1.5, linetype = "dashed", col = "blue") +
  ylab("Fecundity") +
  xlab("Fish Length (mm)")+
  xlim(230,330)+
  scale_y_continuous(limits = c(5000,120000), breaks = seq(10000, 120000, by = 20000))+
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
MESSEIHWF<-(exp(-2.28)*1000)*weight^1.26
PREDICWF<-data.frame(weight,predWF,LCIWF, UCIWF, MESSEIHWF)
##Nice Plot True Values
fecplotWF<-ggplot() + 
  geom_point(data=HDetF, aes(x=fish_weight, y=Fec), size=2, col= "black")+
  geom_ribbon(data=PREDICWF,aes(x=weight,ymin=LCIWF,ymax=UCIWF),alpha=0.25) +
  geom_line(data=PREDICWF,aes(x=weight,y=predWF),size=2) +
  geom_line(data=PREDICWF,aes(x=weight,y=MESSEIHWF),size=1.5, linetype = "dashed", col = "blue") +
  ylab("Fecundity") +
  xlab("Fish Weight (g)")+
  xlim(90,330)+
  scale_y_continuous(limits = c(5000,120000), breaks = seq(10000, 120000, by = 20000))+
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

##Global Length-fecundity and weight-fecundity side by side
ggsave("Fall Length and Weight Fecundity Plot 2022.png", arrangeGrob(fecplotF,fecplotWF, nrow=2), dpi=600, width=22, height=40, unit= "cm")