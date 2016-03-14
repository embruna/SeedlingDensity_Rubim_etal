
#R CODE FOR IMPORTING, MANIPULATING, AND ANALYZING THE DATASETS USED IN: Rubim et al. Seedling Density MS

#These are the last data entry and manipulation" things that need to be dealth with prior to analyses

#Email of August 28

#1 Me explica de novo o que querem dizer os favores nas colunas morre1, morre2, etc.  Acredito que era que se uma planta competidora corria era trocada por outra do viveiro, eh isso? E se era uma planta focal?  
#RESPOSTA: Foi as vezes que houve mortalidade no balde. quando era focal outra foi sorteada  e substitituida.

#In other words, 
#1.a.4 5.a.9 means the first plant was alive peiods intervals 1-4, then replaced with another from 6-9
#DECSION: key is focal plant - so should we do analyes only with 1-4? My gut is YES because it meant that one died.  Or can do for both, and average the RGR?
#NOTE: must make sure #days corresponds to dates of intervals, which durrently don't.
#DECSION 2: #Change ht and la from dead to zero?

#Momentos  Experimento_1	dias

#1	3/8/08	0
#2	4/11/08	34
#3	5/30/08	84
#4	8/31/08	186
#5	12/5/08	282
#6	3/29/09	396
#7	7/16/09	495
#8	1/21/10	684
#9	4/23/10	776


#Experimento_2  dias

#4/15/09	0
#7/16/09	91
#1/21/10	280
#4/23/10	372





##################
###Set WD and load packages you need.
##################

setwd("/Users/emiliobruna/Dropbox/SHARED FOLDERS/Rubim Manuscripts/Seedling density[3]/EB Revision/Data-Rubim-CSV/READY FOR ANALYSIS")
library(reshape2)
library(dplyr)
library(tidyr)
library(ggplot2)


#################
###STEP 1: LOAD THE RAW DATA
##################
#CLear out everything from the environment 
rm(list=ls())
# load the  CSV files and save them as dataframes

canopy<-read.csv("canopy_cover.csv", dec=".", header = TRUE, sep = ",", check.names=FALSE )
bmass<-read.csv("final_sdlg_biomass.csv", dec=".", header = TRUE, sep = ",", check.names=FALSE )
# Make treatment an ordered factor (i.e., 1<2<4)
bmass$trt <- ordered(bmass$trt, levels = c("one", "two", "four"))

# COHORT 1
Exp_Data_C1<-read.csv("EXP_DATA_COHORT1_26nov2014.csv", dec=".", header = TRUE, sep = ",", check.names=FALSE )
# Make treatment an ordered factor (i.e., 1<2<4)
Exp_Data_C1$trt <- ordered(Exp_Data_C1$trt, levels = c("one", "two", "four"))
#Rubim included three more blocks after the experiment started to increase the sample siz - Blocks 17, 18, 19 - but they weren't in the ground
#nearly as long as the others. I would suggest excluding them from the analyses. The following line does that.
Exp_Data_C1<-Exp_Data_C1[Exp_Data_C1$block<17,]
summary(Exp_Data_C1)

# COHORT 2
Exp_Data_C2<-read.csv("EXP_DATA_COHORT2_20nov2014.csv", dec=".", header = TRUE, sep = ",", check.names=FALSE )
# Make treatment an ordered factor (i.e., 1<2<4)
Exp_Data_C2$trt <- ordered(Exp_Data_C2$trt, levels = c("one", "two", "four"))
summary(Exp_Data_C2)



##################
### 
### CALCULATE THE TOTAL LEAF AREA OF EACH SEEDLING
### NOTE THAT STEP 2 and STEP 3 are essentially the same thing, should be made into a function!!!!
###
##################

##################
### COHORT 1
### step 2: reshape data - wide to long - and lineup adjacent columns with 
### leaf lengths and percent of each leaf missing
##################

#select the columns for leaf length 
mdata.la.c1<-Exp_Data_C1[,c(1:82,159)]

#Now melt them with reshape2 package to create a dataframe with 
#leaf-lengths in long form 
mdata.la.c1 <- melt(mdata.la.c1, id.var=c("cohort", "seedling.id.no","block","trt", "sdlg.no", "sdlg.type", "days"))
#rename the column "value" as leaf length
names(mdata.la.c1)[names(mdata.la.c1)=="value"] <- "leaf.length"
#now split the column with leaf number and time interval into 
#two columns to make it easier to sum leaf areas for individual plants 
names(mdata.la.c1)[names(mdata.la.c1)=="variable"] <- "leaf"
split.c1 <- mdata.la.c1$leaf
split.c1<-colsplit(split.c1, ".t", c("leaf", "interval"))
mdata.la.c1<-cbind(mdata.la.c1,split.c1)
#and to keep thinks clean delete the column that you just split in two
mdata.la.c1$leaf <- NULL
#now add a column with % of leaf area missing to that dataframe
#to do so create a dataframe with leaf-areas missing in long form
mdata.miss.c1<-Exp_Data_C1[,c(2,83:158)]  
mdata.miss.c1 <- melt(mdata.miss.c1, id.var=c("seedling.id.no"))
#then rename the column "value" for as "leaf.percentage.missing"
names(mdata.miss.c1)[names(mdata.miss.c1)=="value"] <- "leaf.percentage.missing"
#add that column to the dataframe of leaf lengths and rename that column
cohort1.la<-cbind(mdata.la.c1,mdata.miss.c1$leaf.percentage.missing)
names(cohort1.la)[names(cohort1.la)=="mdata.miss.c1$leaf.percentage.missing"] <- "leaf.percentage.missing"
#the final steps are cleanup: remove NA, which leaves for each plant with records for only the leaves it actually had
#and to add a column identifying which cohort this is (which you will need when you rbind cohort 1 and 2 for analyses)

##################
###step 3: calclulate the area of each leaf, subtract the area missing from each leaf
###and sum to find the total leaf area of each plant in each time interval
##################
#first need to select only the rows that don't have NA in leaf area
cohort1.la<-filter(cohort1.la, leaf.length >=0)
#Adds a column with the leaf area (uncorrected). LA is calculated using the formula in Bruna 2002 Oecologia 
cohort1.la[, "uncorrected.leaf.area"] <- 0.53+(0.831*cohort1.la$leaf.length)
#first removes NA from the column of percet of each leaf missing (in some cases NA because there was no measurment was taken
#then adds a column with the "corrected" leaf area (i.e., corrected leaf area - % missing) 
cohort1.la$leaf.percentage.missing[is.na(cohort1.la$leaf.percentage.missing)] <- 0
cohort1.la[, "corrected.leaf.area"] <-cohort1.la$uncorrected.leaf.area - (cohort1.la$uncorrected.leaf.area*(cohort1.la$leaf.percentage.missing/100))



##################
#### COHORT 2
### step 3: reshape data - wide to long - and lineup adjacent columns with 
### leaf lengths and percent of each leaf missing
##################

#select the columns for leaf length and number of days plant was alive (need to calclulate rgr)
mdata.la<-Exp_Data_C2[,c(1:34,63)]

#Now melt them with reshape2 package to create a dataframe with 
#leaf-lengths in long form 
mdata.la <- melt(mdata.la, id.var=c("cohort","seedling.id.no","block","trt", "sdlg.no", "sdlg.type", "days"))
#rename the column "value" as leaf length
names(mdata.la)[names(mdata.la)=="value"] <- "leaf.length"
#now split the column with leaf number and time interval into 
#two columns to make it easier to sum leaf areas for individual plants 
names(mdata.la)[names(mdata.la)=="variable"] <- "leaf"
split <- mdata.la$leaf
split<-colsplit(split, ".t", c("leaf", "interval"))
mdata.la<-cbind(mdata.la,split)
#and to keep thinks clean delete the column that you just split in two
mdata.la$leaf <- NULL
#now add a column with % of leaf area missing to that dataframe
#to do so create a dataframe with leaf-areas missing in long form
mdata.miss<-Exp_Data_C2[,c(2,35:62)]
mdata.miss <- melt(mdata.miss, id.var=c("seedling.id.no"))
#then rename the column "value" for as "leaf.percentage.missing"
names(mdata.miss)[names(mdata.miss)=="value"] <- "leaf.percentage.missing"
#add that column to the dataframe of leaf lengths and rename that column
cohort2.la<-cbind(mdata.la,mdata.miss$leaf.percentage.missing)
names(cohort2.la)[names(cohort2.la)=="mdata.miss$leaf.percentage.missing"] <- "leaf.percentage.missing"
#the final steps are cleanup: remove NA, which leaves for each plant with records for only the leaves it actually had
#cohort2.la<-na.omit(cohort2.la)

##################
###step 3: calclulate the area of each leaf, subtract the area missing from each leaf
###and sum to find the total leaf area of each plant in each time interval
##################
#Adds a column with the leaf area (uncorrected). LA is calculated using the formula in Bruna 2002 Oecologia 
cohort2.la<-filter(cohort2.la, leaf.length >=0)
cohort2.la$leaf.length<-as.numeric(as.character(cohort2.la$leaf.length))
cohort2.la<-na.omit(cohort2.la)
summary(cohort2.la)

cohort2.la[, "uncorrected.leaf.area"] <- 0.53+(0.831*cohort2.la$leaf.length)




#first removes NA from the column of percet of each leaf missing (in some cases NA because there was no measurment was taken
#then adds a column with the "corrected" leaf area (i.e., corrected leaf area - % missing) 
cohort2.la$leaf.percentage.missing[is.na(cohort2.la$leaf.percentage.missing)] <- 0
cohort2.la[, "corrected.leaf.area"] <-cohort2.la$uncorrected.leaf.area - (cohort2.la$uncorrected.leaf.area*(cohort2.la$leaf.percentage.missing/100))

str(cohort2.la)



##################
###step 4: AGGREGATE DATA AS NEEDED FOR ANALYSIS  
##################

#Aggregate data: summed the leaf-areas for each plant in each sampling interval
COHORT1<-aggregate(cohort1.la$corrected.leaf.area, by=list(cohort1.la$cohort, cohort1.la$block, 
                    cohort1.la$trt, cohort1.la$sdlg.type,cohort1.la$seedling.id.no, cohort1.la$interval,cohort1.la$days), 
                    FUN=sum, na.rm=TRUE)
#rename the columns
names(COHORT1)[1] <- "cohort"
names(COHORT1)[2] <- "block"
names(COHORT1)[3] <- "trt"
names(COHORT1)[4] <- "sdlg.type"
names(COHORT1)[5] <- "sdlg.id.no"
names(COHORT1)[6] <- "interval"
names(COHORT1)[7] <- "days"
names(COHORT1)[8] <- "leaf.area"

COHORT1


COHORT2<-aggregate(cohort2.la$corrected.leaf.area, by=list(cohort2.la$cohort, cohort2.la$block, 
                    cohort2.la$trt, cohort2.la$sdlg.type,cohort2.la$seedling.id.no, cohort2.la$interval, cohort2.la$days), 
                    FUN=sum, na.rm=TRUE)
names(COHORT2)[1] <- "cohort"
names(COHORT2)[2] <- "block"
names(COHORT2)[3] <- "trt"
names(COHORT2)[4] <- "sdlg.type"
names(COHORT2)[5] <- "sdlg.id.no"
names(COHORT2)[6] <- "interval"
names(COHORT2)[7] <- "days"
names(COHORT2)[8] <- "leaf.area"

COHORT2

#using reshape2 to cast the data back into wide form to calclulate RGR (its a data frame, hence dcast)
COHORT1.long<- dcast(COHORT1, cohort + block + trt + sdlg.type + sdlg.id.no + days ~ interval)
COHORT2.long<- dcast(COHORT2, cohort + block + trt + sdlg.type + sdlg.id.no + days ~ interval)


#Sigh. Have to rename the columes of the intervals, it's just easier than remembering to use "1" in all formulas
#renaming columsn cohort 1
names(COHORT1.long)[7] <- "t1"
names(COHORT1.long)[8] <- "t2"
names(COHORT1.long)[9] <- "t3"
names(COHORT1.long)[10] <- "t4"
names(COHORT1.long)[11] <- "t5"
names(COHORT1.long)[12] <- "t6"
names(COHORT1.long)[13] <- "t7"
names(COHORT1.long)[14] <- "t8"
names(COHORT1.long)[15] <- "t9"
#renaming columns cohort 2
names(COHORT2.long)[7] <- "t1"
names(COHORT2.long)[8] <- "t2"
names(COHORT2.long)[9] <- "t3"
names(COHORT2.long)[10] <- "t4"


#A little sorting, just to make it easier to visualize
COHORT1.long<-COHORT1.long[with(COHORT1.long, order(sdlg.id.no,sdlg.type, block, trt)), ]
COHORT2.long<-COHORT2.long[with(COHORT2.long, order(sdlg.id.no,sdlg.type, block, trt)), ]

#################
###CALCLULATIONS OF RGR BAsed on Leaf Area
##################
#Formula for RGR is rgr=[ln(W2)-ln(W1)]/t2-t1, where W is weegith at time t 

# NEED TO FIGURE OUT WHIHC OF THESE IS TRUE - ARE PLANTS WITH ZERO IN LAST COLUMN ALIVE OR DEAD???
# Some of the plants had a leaf area of "0" in the last measuremnt
# If this means the plant is dead, then use this to convert those values to NA so that they don't mess up the calclulation of RGR
COHORT1.long$t9[COHORT1.long$t9== 0] <- NA
COHORT2.long$t4[COHORT2.long$t4== 0] <- NA

# BUT IF THEY ARE STIL ALIVE AND HAVE NO LEAVES, then first calclulate RGR below, then convert the Inf answer to NA using the do.call lines below 

COHORT1.long[, "rgr1.9"] <-(log(COHORT1.long$t9)-log(COHORT1.long$t1))/COHORT1.long$days
#COHORT1.long[, "rgr5.9"] <-(log(COHORT1.long$t9)-log(COHORT1.long$t5))/(COHORT1.long$days-365) #Cohort 1 year 2-3
COHORT1.long[, "rgr1.4"] <-(log(COHORT1.long$t4)-log(COHORT1.long$t1))/COHORT1.long$days
COHORT2.long[, "rgr1.4"] <-(log(COHORT2.long$t4)-log(COHORT2.long$t1))/COHORT2.long$days


COHORT1.long[, "perc.chng.1.9"] <-((COHORT1.long$t9-COHORT1.long$t1)/COHORT1.long$t1)*100
#COHORT1.long[, "perc.chng.5.9"] <-((COHORT1.long$t9-COHORT1.long$t5)/COHORT1.long$t5)*100 #Cohort 1 year 2-3
COHORT1.long[, "perc.chng.1.4"] <-((COHORT1.long$t4-COHORT1.long$t1)/COHORT1.long$t1)*100
COHORT2.long[, "perc.chng.1.4"] <-((COHORT2.long$t4-COHORT2.long$t1)/COHORT2.long$t1)*100




# YOU ONLY NEED THESE IF THE ZEROS FOR PLANT LEAF AREA IN THE LAST TIME INTERVAL MEANS PLANTS SURVIVED BUT HAD NO LEAVES. 
# COHORT1.long<-do.call(data.frame,lapply(COHORT1.long, function(x) replace(x, is.infinite(x),NA)))
# COHORT2.long<-do.call(data.frame,lapply(COHORT1.long, function(x) replace(x, is.infinite(x),NA)))

str(COHORT1.long)
summary(COHORT1.long)

# Box plot of median leaf area of plants in each block at T1 
boxplot(t1~block,data=COHORT1.long) #Cohort 1
boxplot(t1~block,data=COHORT2.long) #Cohort2


hist(COHORT1.long$rgr1.9)

hist(COHORT1.long$rgr1.4)

hist(COHORT2.long$rgr1.4)

#BOX PLOT INCLUDING ALL PLANTS 

rgrALL1.4 <- ggplot(COHORT1.long, aes(x=trt, y=rgr1.4)) + 
  geom_boxplot()
rgrALL1.4


rgrALL1.9 <- ggplot(COHORT1.long, aes(x=trt, y=rgr1.9)) + 
  geom_boxplot()
rgrALL1.9


rgrALL2 <- ggplot(COHORT2.long, aes(x=trt, y=rgr1.4)) + 
  geom_boxplot()
rgrALL2




#######PUTTING BOTH COHORTS TOGETHER

LFAREA1<-COHORT1.long


# DELETE UNECESSARY COLUMNS
LFAREA1.1 <- LFAREA1[ -c(7:15, 18:19)]
LFAREA1.2 <- LFAREA1[ -c(7:17)]


str(LFAREA1.1)

#COnvert each to long 
foo<-gather(LFAREA1.1, "interval", "rgr.la", 7:8)  
foo2<-gather(LFAREA1.2, "interval", "perc.chng.la", 7:8)
foo3<-cbind(foo,foo2)
foo3<- foo3[ -c(9:15)]
#chnage the values of the cells
foo3$interval <-as.character(foo3$interval) #need these because they are factors, covert to chars
foo3$interval <- replace(foo3$interval, foo3$interval=="rgr1.9", "24 mos")
foo3$interval <- replace(foo3$interval, foo3$interval=="rgr1.4", "12 mos")
#delete interval 2
str(foo3)

# 
# Cohort 2
LFAREA2<-COHORT2.long

# DELETE UNECESSARY COLUMNS
LFAREA2.1 <- LFAREA2[ -c(7:10,12)]
LFAREA2.2 <- LFAREA2[ -c(7:11)]
str(LFAREA2.1)
str(LFAREA2.2)

#COnvert each to long 
foo4<-gather(LFAREA2.1, "interval", "rgr.la", 7)  
foo5<-gather(LFAREA2.2, "interval", "perc.chng.la", 7)
foo6<-cbind(foo4,foo5)
foo6<- foo6[ -c(9:15)]
#chnage the values of the cells
foo6$interval <-as.character(foo6$interval) #need these because they are factors, covert to chars
foo6$interval <- replace(foo6$interval, foo6$interval=="rgr1.4", "12 mos")
#delete interval 2
str(foo3)
str(foo6)
foo7<-rbind(foo3,foo6)
#add canopy
foo7<-left_join(foo7,canopy[1:2], by = "block")
foo7$cohort<-as.factor(foo7$cohort)
#Reduce dataset: only include "focal" seedlings.


Focal<-filter(foo7, sdlg.type == "focal")
Focal<-droplevels(Focal)
Focal<-filter(Focal, interval == "12 mos")
Focal<-droplevels(Focal)
aov1<-aov(rgr.la ~ trt+cohort+block, data = Focal)
hist(foo7$rgr.la)
# hist(foo7$perc.chng.la)
summary(aov1)
bmass1Fig <- ggplot(Focal, aes(x=cohort, y=rgr.la)) + 
  geom_boxplot()
bmass1Fig

bmass1Fig <- ggplot(Focal, aes(x=cohort, y=perc.chng.la)) + 
  geom_boxplot()
bmass1Fig



#################
###CALCLULATIONS OF RGR BAsed on seedling height
##################
#rgr.ht=(log(Exp_Data$ht.final)-log(Exp_Data$ht.initial))/Exp_Data$days

# COHORT 1
str(Exp_Data_C1)
ht.cohort1<-Exp_Data_C1

# DELETE UNECESSARY COLUMNS
ht.cohort1 <- ht.cohort1[ -c(7:158,169:181)]
str(ht.cohort1)

# CALC of RGR 
ht.cohort1[, "rgr.ht.1.4"] <-(log(ht.cohort1$ht.4)-log(ht.cohort1$ht.1))/ht.cohort1$days
ht.cohort1[, "rgr.ht.1.9"] <-(log(ht.cohort1$ht.9)-log(ht.cohort1$ht.1))/ht.cohort1$days

ht.cohort1<-ht.cohort1[with(ht.cohort1, order(seedling.id.no, block, trt, sdlg.type)), ]
str(ht.cohort1)
summary(ht.cohort1)

# COHORT 2
str(Exp_Data_C2)
ht.cohort2<-Exp_Data_C2

# DELETE UNECESSARY COLUMNS
ht.cohort2 <- ht.cohort2[ -c(7:62,68:73)]

# CALC of RGR 
ht.cohort2[, "rgr.ht.1.4"] <-(log(ht.cohort2$ht.4)-log(ht.cohort2$ht.1))/ht.cohort2$days
str(ht.cohort2)
ht.cohort2<-ht.cohort2[with(ht.cohort2, order(seedling.id.no, block, trt, sdlg.type)), ]
str(ht.cohort2)

## LONG FORM THEN JOIN THEM
# DELETE UNECESSARY COLUMNS
ht.cohort1.1 <- ht.cohort1[ -c(8:16, 18)]
ht.cohort1.2 <- ht.cohort1[ -c(8:17)]
ht.cohort2 <- ht.cohort2[ -c(8:11)]
#COnvert each to long 

ht.cohort1.1<-gather(ht.cohort1.1, "interval", "rgr.ht", 8)  
ht.cohort1.2<-gather(ht.cohort1.2, "interval", "rgr.ht", 8)  
ht.cohort2<-gather(ht.cohort2, "interval", "rgr.ht", 8)  

#chnage the values of the cells
ht.cohort1.1$interval <-as.character(ht.cohort1.1$interval) #need these because they are factors, covert to chars
ht.cohort1.1$interval <- replace(ht.cohort1.1$interval, ht.cohort1.1$interval=="rgr.ht.1.4", "12 mos")

#chnage the values of the cells
ht.cohort1.2$interval <-as.character(ht.cohort1.2$interval) #need these because they are factors, covert to chars
ht.cohort1.2$interval <- replace(ht.cohort1.2$interval, ht.cohort1.2$interval=="rgr.ht.1.9", "24 mos")


#chnage the values of the cells
ht.cohort2$interval <-as.character(ht.cohort2$interval) #need these because they are factors, covert to chars
ht.cohort2$interval <- replace(ht.cohort2$interval, ht.cohort2$interval=="rgr.ht.1.4", "12 mos")
str(ht.cohort1.1)
str(ht.cohort1.2)
str(ht.cohort2)
all.height.rgr<-rbind(ht.cohort1.2,ht.cohort1.1,ht.cohort2)


dim(all.height.rgr)
dim(foo7)
str(foo7)

all.rgr<-cbind(foo7,all.height.rgr)
all.rgr <- all.rgr[ -c(11:18)]

plot(all.rgr$rgr.ht,all.rgr$rgr.la)






###################################################
######      ANALYSES - RGR  ###
######################################################
#Reduce dataset: only include "focal" seedlings.
rgr.test<-filter(all.rgr, sdlg.type == "focal")
rgr.test<-droplevels(rgr.test)


#SIMPLER AS ANOVA
aov.la<-aov(rgr.la ~ trt*cohort+block, data = rgr.test)
summary(aov.la)

aov.ht<-aov(rgr.ht ~ trt*cohort+block, data = rgr.test)
summary(aov.ht)

plot(rgr.test$canopy.openess.percent,rgr.test$rgr.ht)

# GLM LA

glm.1<-glm(rgr.la ~ block, family = gaussian, data = rgr.test)
summary(glm.1)
glm.2<-glm(rgr.la ~ trt+block, family = gaussian, data = rgr.test)
summary(glm.2)
glm.3<-glm(rgr.la ~ cohort+block, family = gaussian, data = rgr.test)
summary(glm.3)
glm.4<-glm(rgr.la ~ trt+cohort+block, family = gaussian, data = rgr.test)
summary(glm.4)

anova(glm.1,glm.2, test = "Chisq")  #no imp fit by adding trt over just block
anova(glm.1,glm.3, test = "Chisq")  #sig improve fit when adding cohort
anova(glm.3,glm.4, test = "Chisq")  #adding trt means nothing
anova(glm.3,glm.2, test = "Chisq")  ##this just emphasizes the importance of cohort, not trt

AIC(glm.1, glm.2, glm.3,glm.4)

# GLM - Height

glm.1<-glm(rgr.ht ~ block, family = gaussian, data = rgr.test)
summary(glm.1)
glm.2<-glm(rgr.ht ~ trt+block, family = gaussian, data = rgr.test)
summary(glm.2)
glm.3<-glm(rgr.ht ~ cohort+block, family = gaussian, data = rgr.test)
summary(glm.3)
glm.4<-glm(rgr.ht ~ trt+cohort+block, family = gaussian, data = rgr.test)
summary(glm.4)

anova(glm.1,glm.2, test = "Chisq")  #no imp fit by adding trt over just block
anova(glm.1,glm.3, test = "Chisq")  #sig improve fit when adding cohort
anova(glm.3,glm.4, test = "Chisq")  #adding trt means nothing
anova(glm.3,glm.2, test = "Chisq")  ##this just emphasizes the importance of cohort, not trt

AIC(glm.1, glm.2, glm.3,glm.4)


###################################################
######      ANALYSES - BMASS  ###
######################################################

bmass[, "RSratio"] <-(bmass$root.biomass/(bmass$lf.biomass+bmass$stem.biomass))
bmass[, "Total.bmass"] <-bmass$root.biomass+bmass$lf.biomass+bmass$stem.biomass
bmass<-full_join(bmass, canopy, by = "block")
bmass$location<-NULL
bmass$site.openess.percent<-NULL
bmass$mask.openess.percent<-NULL
bmass$sky.area.percent<-NULL
bmass$LAI.4.ring<-NULL 
bmass$LAI.5.ring<-NULL
str(bmass)

glm.1<-glm(Total.bmass ~ block, family = identity, data = bmass)
summary(glm.1)
glm.2<-glm(Total.bmass ~ trt+block, family = log, data = bmass)
summary(glm.2)
glm.3<-glm(Total.bmass ~ cohort+block, family = log, data = bmass)
summary(glm.3)
glm.4<-glm(Total.bmass ~ trt+cohort+block, family = log, data = bmass)
summary(glm.4)

anova(glm.1,glm.2, test = "Chisq")  #no imp fit by adding trt over just block
anova(glm.1,glm.3, test = "Chisq")  #sig improve fit when adding cohort
anova(glm.3,glm.4, test = "Chisq")  #adding trt means nothing
# anova(glm.3,glm.2, test = "Chisq")  ##this just emphasizes the importance of cohort, not trt

AIC(glm.1, glm.2, glm.3,glm.4)

#SIMPLER AS ANOVA
aov.bm<-aov(RSratio ~ trt*cohort+block, data = bmass)
summary(aov.bm)














# logit trasnform, but actually don't do this...
# bmass<-bmass[, "logitRS"] <-log10(bmass$RSratio+/(1-bmass$RSratio))



# SPlit the cohorts 
bmass.1<-filter(bmass, cohort == "1")
bmass.1<-droplevels(bmass.1)
str(bmass.1)
hist(bmass.1$RSratio)
aov.bmass.1<-aov(RSratio ~ trt+block, data = bmass.1)
summary(aov.bmass.1)


bmass.2<-filter(bmass, cohort == "2")
bmass.2<-droplevels(bmass.2)
hist(bmass.2$RSratio)
aov.bmass.2<-aov(RSratio ~ trt+block, data = bmass.2)
summary(aov.bmass.2)


# Graph 
bmass1Fig <- ggplot(bmass.1, aes(x=trt, y=RSratio)) + 
  geom_boxplot()
bmass1Fig

bmass2Fig <- ggplot(bmass.2, aes(x=trt, y=RSratio)) + 
  geom_boxplot()
bmass2Fig

# BOX PLOTS
boxplot(rgr.ht.1.4~trt,data=ht.cohort2) #Cohort 2 1.4
boxplot(ht.4~trt,data=ht.cohort2) #Cohort 1 ht final after 4























###################################################
###
###      ANALYSES - RGR Based on LEAF AREA     ###
###
###################################################

# COHORT ONE
#Reduce dataset: only include "focal" seedlings.
Focal.One<-filter(COHORT1.long, sdlg.type == "focal")
Focal.One<-droplevels(Focal.One)
Focal.One<-full_join(Focal.One, canopy, by = "block")
summary(Focal.One)

glm.1<-glm(rgr1.9 ~ block, family = gaussian, data = Focal.One)
summary(glm.1)
glm.2<-glm(rgr1.9 ~ trt+block, family = gaussian, data = Focal.One)
summary(glm.2)
anova(glm.1,glm.2, test = "Chisq")
AIC(glm.1, glm.2)

#SIMPLER AS ANOVA
aov1<-aov(rgr1.9 ~ trt+block, data = Focal.One)
summary(aov1)
plot(Focal.One$canopy.openess.percent, Focal.One$rgr1.9)

# repeated measures ANOVA with growth after 1 year, then after 2 years.
#
#
#
# NEED TO CODE THE RM ANOVA
#
#
#




# COHORT TWO
#Reduce dataset: only include "focal" seedlings.
Focal.Two<-filter(COHORT2.long, sdlg.type == "focal")
Focal.Two<-droplevels(Focal.Two)
Focal.Two<-full_join(Focal.Two, canopy, by = "block")
summary(Focal.Two)

glm.3<-glm(rgr1.4 ~ block, family = gaussian, data = Focal.Two)
summary(glm.3)
glm.4<-glm(rgr1.4 ~ trt+block, family = gaussian, data = Focal.Two)
summary(glm.4)
anova(glm.3,glm.4, test = "Chisq")
AIC(glm.3, glm.4)

#SIMPLER AS ANOVA
aov2<-aov(rgr1.4 ~ trt+block, data = Focal.Two)
summary(aov2)

#####################
# BIOMASS

bmass[, "RSratio"] <-(bmass$root.biomass/(bmass$lf.biomass+bmass$stem.biomass))
bmass<-full_join(bmass, canopy, by = "block")
bmass<-select(bmass, block:canopy.openess.percent)

# logit trasnform, but actually don't do this...
# bmass<-bmass[, "logitRS"] <-log10(bmass$RSratio+/(1-bmass$RSratio))

# SPlit the cohorts 
bmass.1<-filter(bmass, cohort == "1")
bmass.1<-droplevels(bmass.1)
str(bmass.1)
hist(bmass.1$RSratio)
aov.bmass.1<-aov(RSratio ~ trt+block, data = bmass.1)
summary(aov.bmass.1)


bmass.2<-filter(bmass, cohort == "2")
bmass.2<-droplevels(bmass.2)
hist(bmass.2$RSratio)
aov.bmass.2<-aov(RSratio ~ trt+block, data = bmass.2)
summary(aov.bmass.2)


# Graph 
bmass1Fig <- ggplot(bmass.1, aes(x=trt, y=RSratio)) + 
  geom_boxplot()
bmass1Fig

bmass2Fig <- ggplot(bmass.2, aes(x=trt, y=RSratio)) + 
  geom_boxplot()
bmass2Fig


# BOX PLOTS
boxplot(rgr.ht.1.4~trt,data=ht.cohort2) #Cohort 2 1.4
boxplot(ht.4~trt,data=ht.cohort2) #Cohort 1 ht final after 4


# ANOVA
ht.cohort2.focal<-filter(ht.cohort2, sdlg.type == "focal")
ht.cohort2.focal<-droplevels(ht.cohort2.focal)

boxplot(ht.4~trt,data=ht.cohort2.focal) #Cohort 1 ht focal final after 4

aov.C2.14<-aov(rgr.ht.1.4 ~ trt+block, data = ht.cohort2.focal)
summary(aov.C2.14)

# BOX PLOTS
boxplot(rgr.ht.1.4~block,data=ht.cohort1) #Cohort 1 1.4
boxplot(rgr.ht.1.9~block,data=ht.cohort1) #Cohort 1 1.9
boxplot(ht.4~trt,data=ht.cohort1) #Cohort 1 ht final after 4
boxplot(ht.9~trt,data=ht.cohort1) #Cohort 1 ht final after 9

# ANOVA
# First reduce to only "focal' seedlings
ht.cohort1.focal<-filter(ht.cohort1, sdlg.type == "focal")
ht.cohort1.focal<-droplevels(ht.cohort1.focal)


boxplot(ht.4~trt,data=ht.cohort1.focal) #Cohort 1 focal ht final after 4
boxplot(ht.9~trt,data=ht.cohort1.focal) #Cohort 1 ht focal final after 9


aov.C1.14<-aov(rgr.ht.1.4 ~ trt+block, data = ht.cohort1.focal)
summary(aov.C1.14)

aov.C1.19<-aov(rgr.ht.1.9 ~ trt+block, data = ht.cohort1.focal)
summary(aov.C1.19)




############


# SEEDLINGS FROM THE DEMOG PLOTS
# #NEED TO IDENTIFY WHICH IS THE ONE THAT IS THE FOCAL ONE IN THE PAIR!!!


plot.sdlgs<-read.csv("demog_plot_sdlgs.csv", dec=".", header = TRUE, sep = ",", check.names=FALSE )
plot.sdlgs[, "rgr.ht.0809"] <-(log(plot.sdlgs$ht.2009)-log(plot.sdlgs$ht.2008))/365
plot.sdlgs[, "rgr.ht.0810"] <-(log(plot.sdlgs$ht.2010)-log(plot.sdlgs$ht.2008))/365

boxplot(rgr.ht.0810~trt,data=plot.sdlgs) #08 sdeedlings rgr 0810
boxplot(ht.2010~trt,data=plot.sdlgs) #08 sflgs final ht 2010



aov.plot.0810<-aov(rgr.ht.0810 ~ trt+light, data = plot.sdlgs)
summary(aov.plot.0810)

plot(plot.sdlgs$light,plot.sdlgs$rgr.ht.0810)




