
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
library("lme4")
library("nlme")

#################
###STEP 1: LOAD THE RAW DATA
##################
# Clear the environment 
rm(list=ls())

# load the data on canopy cover and biomass
canopy<-read.csv("canopy_cover.csv", dec=".", header = TRUE, sep = ",", check.names=FALSE )
bmass<-read.csv("final_sdlg_biomass.csv", dec=".", header = TRUE, sep = ",", check.names=FALSE )
# Make treatment an ordered factor (i.e., 1<2<4)
# Exclude blocks 17-19 from Cohort 1. These were added late in the experiment and are not a valid comparison
bmass$trt <- ordered(bmass$trt, levels = c("one", "two", "four"))
bmass<-bmass[!(bmass$cohort == "1" & bmass$block>16), ]

# load data on COHORT 1
Exp_Data_C1<-read.csv("EXP_DATA_COHORT1_26nov2014.csv", dec=".", header = TRUE, sep = ",", check.names=FALSE )
# Make treatment an ordered factor (i.e., 1<2<4)
Exp_Data_C1$trt <- ordered(Exp_Data_C1$trt, levels = c("one", "two", "four"))
# Exclude blocks 17-19 from Cohort 1. These were added late in the experiment and are not a valid comparison
Exp_Data_C1<-Exp_Data_C1[Exp_Data_C1$block<17,]
summary(Exp_Data_C1)

# load data on COHORT 2
Exp_Data_C2<-read.csv("EXP_DATA_COHORT2_20nov2014.csv", dec=".", header = TRUE, sep = ",", check.names=FALSE )
# Make treatment an ordered factor (i.e., 1<2<4)
Exp_Data_C2$trt <- ordered(Exp_Data_C2$trt, levels = c("one", "two", "four"))
# summary(Exp_Data_C2)


##################
### STEP 2: CALCULATE THE TOTAL LEAF AREA OF EACH SEEDLING
### NOTE THAT STEP 2 and STEP 3 are essentially the same thing, should be made into a function!!!!
##################

##################
### COHORT 1: reshape data (wide to long) & lineup adjacent columns with leaf lengths & % of each leaf missing
##################

# select the columns for leaf length 
mdata.la.c1<-Exp_Data_C1[,c(1:82,159)]
# melt them with reshape2 package to create a dataframe with leaf-lengths in long form 
mdata.la.c1 <- melt(mdata.la.c1, id.var=c("cohort", "seedling.id.no","block","trt", "sdlg.no", "sdlg.type", "days"))
# rename the column "value" as leaf length
names(mdata.la.c1)[names(mdata.la.c1)=="value"] <- "leaf.length"
# split column with leaf number & time interval into 2 columns to make it easier to sum leaf areas of individual plants 
names(mdata.la.c1)[names(mdata.la.c1)=="variable"] <- "leaf"
split.c1 <- mdata.la.c1$leaf
split.c1<-colsplit(split.c1, ".t", c("leaf", "interval"))
mdata.la.c1<-cbind(mdata.la.c1,split.c1)
# delete the column that you just split in two
mdata.la.c1$leaf <- NULL
# add a column with % of leaf area missing to that dataframe
# to do so create a dataframe with leaf-areas missing in long form
mdata.miss.c1<-Exp_Data_C1[,c(2,83:158)]  
mdata.miss.c1 <- melt(mdata.miss.c1, id.var=c("seedling.id.no"))
# rename the column "value" for as "leaf.percentage.missing"
names(mdata.miss.c1)[names(mdata.miss.c1)=="value"] <- "leaf.percentage.missing"
# add that column to the dataframe of leaf lengths and rename that column
cohort1.la<-cbind(mdata.la.c1,mdata.miss.c1$leaf.percentage.missing)
names(cohort1.la)[names(cohort1.la)=="mdata.miss.c1$leaf.percentage.missing"] <- "leaf.percentage.missing"
# the final steps are cleanup: remove NA, which leaves for each plant with records for only the leaves it actually had
# and to add a column identifying which cohort this is (which you will need when you rbind cohort 1 and 2 for analyses)

##################
### STEP 3 (cohort 1): calclulate the area of each leaf, subtract the area missing from each leaf
### and sum to find the total leaf area of each plant in each time interval
##################
# first need to select only the rows that don't have NA in leaf area
cohort1.la<-filter(cohort1.la, leaf.length >=0)
# Adds a column with the leaf area (uncorrected). LA is calculated using the formula in Bruna 2002 Oecologia 
cohort1.la[, "uncorrected.leaf.area"] <- 0.53+(0.831*cohort1.la$leaf.length)
# first removes NA from the column of percet of each leaf missing (in some cases NA because there was no measurment was taken
# then adds a column with the "corrected" leaf area (i.e., corrected leaf area - % missing) 
cohort1.la$leaf.percentage.missing[is.na(cohort1.la$leaf.percentage.missing)] <- 0
cohort1.la[, "corrected.leaf.area"] <-cohort1.la$uncorrected.leaf.area - (cohort1.la$uncorrected.leaf.area*(cohort1.la$leaf.percentage.missing/100))




##################
### COHORT 2: reshape data (wide to long) & lineup adjacent columns with leaf lengths & % of each leaf missing
##################

# select the columns for leaf length and number of days plant was alive (need to calclulate rgr)
mdata.la<-Exp_Data_C2[,c(1:34,63)]

# melt them with reshape2 package to create a dataframe with leaf-lengths in long form 
mdata.la <- melt(mdata.la, id.var=c("cohort","seedling.id.no","block","trt", "sdlg.no", "sdlg.type", "days"))
# rename the column "value" as leaf length
names(mdata.la)[names(mdata.la)=="value"] <- "leaf.length"
# split column with leaf number & time interval into 2 columns to make it easier to sum leaf areas for individual plants 
names(mdata.la)[names(mdata.la)=="variable"] <- "leaf"
split <- mdata.la$leaf
split<-colsplit(split, ".t", c("leaf", "interval"))
mdata.la<-cbind(mdata.la,split)
# to keep thinks clean delete the column that you just split in two
mdata.la$leaf <- NULL
# add a column with % of leaf area missing to that dataframe by creating dataframe with leaf-areas missing in long form
mdata.miss<-Exp_Data_C2[,c(2,35:62)]
mdata.miss <- melt(mdata.miss, id.var=c("seedling.id.no"))
# rename the column "value" for as "leaf.percentage.missing"
names(mdata.miss)[names(mdata.miss)=="value"] <- "leaf.percentage.missing"
# add that column to the dataframe of leaf lengths and rename that column
cohort2.la<-cbind(mdata.la,mdata.miss$leaf.percentage.missing)
names(cohort2.la)[names(cohort2.la)=="mdata.miss$leaf.percentage.missing"] <- "leaf.percentage.missing"
# the final steps are cleanup: remove NA, which leaves for each plant with records for only the leaves it actually had
# cohort2.la<-na.omit(cohort2.la)

##################
###STEP 3 (Cohort 2): calclulate the area of each leaf, subtract the area missing from each leaf and sum to find the total leaf area of each plant in each time interval
##################

# Add a column with the leaf area (uncorrected). LA is calculated using the formula in Bruna 2002 Oecologia 
cohort2.la<-filter(cohort2.la, leaf.length >=0)
cohort2.la$leaf.length<-as.numeric(as.character(cohort2.la$leaf.length))
cohort2.la<-na.omit(cohort2.la)
# summary(cohort2.la)
cohort2.la[, "uncorrected.leaf.area"] <- 0.53+(0.831*cohort2.la$leaf.length)

# first removes NA from the column of percet of each leaf missing (in some cases NA because there was no measurment was taken
# then adds a column with the "corrected" leaf area (i.e., corrected leaf area - % missing) 
cohort2.la$leaf.percentage.missing[is.na(cohort2.la$leaf.percentage.missing)] <- 0
cohort2.la[, "corrected.leaf.area"] <-cohort2.la$uncorrected.leaf.area - (cohort2.la$uncorrected.leaf.area*(cohort2.la$leaf.percentage.missing/100))
# str(cohort2.la)


######################################################
### STEP 4: AGGREGATE DATA FOR ANALYSIS  
######################################################

# Aggregate data: summed the leaf-areas for each plant in each sampling interval
COHORT1<-aggregate(cohort1.la$corrected.leaf.area, by=list(cohort1.la$cohort, cohort1.la$block, 
                    cohort1.la$trt, cohort1.la$sdlg.type,cohort1.la$seedling.id.no,  cohort1.la$interval,cohort1.la$days), 
                    FUN=sum, na.rm=TRUE)
# rename the columns
names(COHORT1)[1] <- "cohort"
names(COHORT1)[2] <- "block"
names(COHORT1)[3] <- "trt"
names(COHORT1)[4] <- "sdlg.type"
names(COHORT1)[5] <- "sdlg.id.no"
names(COHORT1)[6] <- "interval"
names(COHORT1)[7] <- "days"
names(COHORT1)[8] <- "leaf.area"

# COHORT1


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

# COHORT2

# using reshape2 cast the data back into wide form to calclulate RGR (its a data frame, hence dcast)
COHORT1.wide<- dcast(COHORT1, cohort + block + trt + sdlg.type + sdlg.id.no + days ~ interval)
COHORT2.wide<- dcast(COHORT2, cohort + block + trt + sdlg.type + sdlg.id.no + days ~ interval)


# Sigh. Have to rename the columes of the intervals, it's just easier than remembering to use "1" in all formulas
# renaming columns for Cohort 1
names(COHORT1.wide)[7] <- "LAt1"
names(COHORT1.wide)[8] <- "LAt2"
names(COHORT1.wide)[9] <- "LAt3"
names(COHORT1.wide)[10] <- "LAt4"
names(COHORT1.wide)[11] <- "LAt5"
names(COHORT1.wide)[12] <- "LAt6"
names(COHORT1.wide)[13] <- "LAt7"
names(COHORT1.wide)[14] <- "LAt8"
names(COHORT1.wide)[15] <- "LAt9"

# renaming columns for Cohort 2
names(COHORT2.wide)[7] <- "LAt1"
names(COHORT2.wide)[8] <- "LAt2"
names(COHORT2.wide)[9] <- "LAt3"
names(COHORT2.wide)[10] <- "LAt4"

# A little sorting to make it easier to visualize
COHORT1.wide<-COHORT1.wide[with(COHORT1.wide, order(sdlg.id.no,sdlg.type, block, trt)), ]
COHORT2.wide<-COHORT2.wide[with(COHORT2.wide, order(sdlg.id.no,sdlg.type, block, trt)), ]

###################################################
### CALCLULATIONS OF RGR BAsed on Leaf Area
####################################################

# Formula for RGR is rgr=[ln(W2)-ln(W1)]/t2-t1, where W is weegith at time t 

# NEED TO FIGURE OUT WHIHC OF THESE IS TRUE - ARE PLANTS WITH ZERO IN LAST COLUMN ALIVE OR DEAD???
# Some of the plants had a leaf area of "0" in the last measuremnt
# If this means the plant is dead, then use this to convert those values to NA so that they don't mess up the calclulation of RGR
COHORT1.wide$LAt9[COHORT1.wide$LAt9== 0] <- NA
COHORT2.wide$LAt4[COHORT2.wide$LAt4== 0] <- NA

# BUT IF THEY ARE STIL ALIVE AND HAVE NO LEAVES, then first calclulate RGR below, then convert the Inf answer to NA using the do.call lines below 

COHORT1.wide[, "rgrLA_wet1"] <-(log(COHORT1.wide$LAt3)-log(COHORT1.wide$LAt1))/84 #t1-t3
COHORT1.wide[, "rgrLA_dry1"] <-(log(COHORT1.wide$LAt5)-log(COHORT1.wide$LAt3))/(282-84) #t3-t5
COHORT1.wide[, "rgrLA_wet2"] <-(log(COHORT1.wide$LAt7)-log(COHORT1.wide$LAt5))/(495-282) #t5-t7
COHORT1.wide[, "rgrLA_dry2"] <-(log(COHORT1.wide$LAt8)-log(COHORT1.wide$LAt7))/(684-495) #t7-t8
COHORT1.wide[, "rgrLA_wet3"] <-(log(COHORT1.wide$LAt9)-log(COHORT1.wide$LAt8))/(COHORT1.wide$days-684) #t7-t8
COHORT1.wide[, "rgrLA_yr1"] <-(log(COHORT1.wide$LAt6)-log(COHORT1.wide$LAt1))/369    #5 08 march 08 - 29  march 09 
COHORT1.wide[, "rgrLA_yr2"] <-(log(COHORT1.wide$LAt9)-log(COHORT1.wide$LAt6))/(COHORT1.wide$days-396) #March 8 08 
COHORT1.wide[, "rgrLA_yrs1-2"] <-(log(COHORT1.wide$LAt9)-log(COHORT1.wide$LAt1))/COHORT1.wide$days 

COHORT2.wide[, "rgrLA_wet2"] <-(log(COHORT2.wide$LAt2)-log(COHORT2.wide$LAt1))/91
COHORT2.wide[, "rgrLA_dry2"] <-(log(COHORT2.wide$LAt3)-log(COHORT2.wide$LAt1))/(280-91)
COHORT2.wide[, "rgrLA_wet3"] <-(log(COHORT2.wide$LAt4)-log(COHORT2.wide$LAt3))/(COHORT2.wide$days-280)
COHORT2.wide[, "rgrLA_yr1"] <-(log(COHORT2.wide$LAt4)-log(COHORT2.wide$LAt1))/COHORT2.wide$days

# YOU ONLY NEED THESE IF THE ZEROS FOR PLANT LEAF AREA IN THE LAST TIME INTERVAL MEANS PLANTS SURVIVED BUT HAD NO LEAVES. 
# COHORT1.wide<-do.call(data.frame,lapply(COHORT1.wide, function(x) replace(x, is.infinite(x),NA)))
# COHORT2.wide<-do.call(data.frame,lapply(COHORT1.wide, function(x) replace(x, is.infinite(x),NA)))

# str(COHORT1.wide)
# summary(COHORT1.wide)

# #################
# ###CALCLULATIONS OF RGR Based on seedling height
# ##################
# #rgr.ht=(log(Exp_Data$ht.final)-log(Exp_Data$ht.initial))/Exp_Data$days
# 
# # ########
# # COHORT 1
# # ########
#
ht.cohort1<-Exp_Data_C1

# DELETE UNECESSARY COLUMNS
ht.cohort1 <- ht.cohort1[ -c(7:158,169:181)]

ht.cohort2<-Exp_Data_C2
# DELETE UNECESSARY COLUMNS
ht.cohort2 <- ht.cohort2[ -c(7:62,68:73)]

# CALC of RGR
ht.cohort1[, "rgrHT_wet1"] <-(log(ht.cohort1$ht.3)-log(ht.cohort1$ht.1))/84 #t1-t3
ht.cohort1[, "rgrHT_dry1"] <-(log(ht.cohort1$ht.5)-log(ht.cohort1$ht.3))/(282-84) #t3-t5
ht.cohort1[, "rgrHT_wet2"] <-(log(ht.cohort1$ht.7)-log(ht.cohort1$ht.5))/(495-282) #t5-t7
ht.cohort1[, "rgrHT_dry2"] <-(log(ht.cohort1$ht.8)-log(ht.cohort1$ht.7))/(684-495) #t7-t8
ht.cohort1[, "rgrHT_wet3"] <-(log(ht.cohort1$ht.9)-log(ht.cohort1$ht.8))/(ht.cohort1$days-684) #t7-t8
ht.cohort1[, "rgrHT_yr1"] <-(log(ht.cohort1$ht.6)-log(ht.cohort1$ht.1))/369    #5 08 march 08 - 29  march 09
ht.cohort1[, "rgrHT_yr2"] <-(log(ht.cohort1$ht.9)-log(ht.cohort1$ht.6))/(ht.cohort1$days-396) #March 8 08
ht.cohort1[, "rgrHT_yrs1-2"] <-(log(ht.cohort1$ht.9)-log(ht.cohort1$ht.1))/ht.cohort1$days
ht.cohort2[, "rgrHT_wet2"] <-(log(ht.cohort2$ht.2)-log(ht.cohort2$ht.1))/91
ht.cohort2[, "rgrHT_dry2"] <-(log(ht.cohort2$ht.3)-log(ht.cohort2$ht.1))/(280-91)
ht.cohort2[, "rgrHT_wet3"] <-(log(ht.cohort2$ht.4)-log(ht.cohort2$ht.3))/(ht.cohort2$days-280)
ht.cohort2[, "rgrHT_yr1"] <-(log(ht.cohort2$ht.4)-log(ht.cohort2$ht.1))/ht.cohort2$days

ht.cohort1.1 <- ht.cohort1[ -c(8:16)]
ht.cohort1.1<-gather(ht.cohort1.1, "interval", "rgr.ht", 8:15)

ht.cohort2.1 <- ht.cohort2[ -c(8:11)]
ht.cohort2.1<-gather(ht.cohort2.1, "interval", "rgr.ht", 8:11)

str(ht.cohort1.1)
str(ht.cohort2.1)


###################################################
######      PREP DATA FOR ANALYSES  --- LEAF AREA
######################################################

# Cohort 1

COHORT1.wide<-COHORT1.wide[with(COHORT1.wide, order(sdlg.id.no, sdlg.type, block, trt)), ]
ht.cohort1<-ht.cohort1[with(ht.cohort1, order(seedling.id.no, sdlg.type, block, trt, days)), ]
sdlg_size_start<-as.data.frame(COHORT1.wide$LAt1)
sdlg_size_middle<-as.data.frame(COHORT1.wide$LAt6)
sdlg_size_end<-as.data.frame(COHORT1.wide$LAt9)
sdlg_size_end<-cbind(sdlg_size_start,sdlg_size_middle, sdlg_size_end, as.data.frame(ht.cohort1$ht.9),as.data.frame(ht.cohort1$trt),as.data.frame(ht.cohort1$block),as.data.frame(ht.cohort1$cohort),as.data.frame(ht.cohort1$days))
names(sdlg_size_end)[1] <- "LA_initial"
names(sdlg_size_end)[2] <- "LA_middle"
names(sdlg_size_end)[3] <- "LA_final"
names(sdlg_size_end)[4] <- "HT_final"
names(sdlg_size_end)[5] <- "trt"
names(sdlg_size_end)[6] <- "block"
names(sdlg_size_end)[7] <- "cohort"
names(sdlg_size_end)[8] <- "days"
na.omit(sdlg_size_end)

# Cohort 2
COHORT2.wide<-COHORT2.wide[with(COHORT2.wide, order(sdlg.id.no, sdlg.type, block, trt)), ]
ht.cohort2<-ht.cohort2[with(ht.cohort2, order(seedling.id.no, sdlg.type, block, trt, days)), ]

# Add the initial seedling size
as.factor(ht.cohort2$sdlg.id.no)
full_join(ht.cohort2,COHORT2.wide,by="sdlg.id.no")
sdlg_size_C2start<-as.data.frame(COHORT2.wide$LAt1)
sdlg_size_C2end<-as.data.frame(COHORT2.wide$LAt4)
cbind(sdlg_size_C2start,sdlg_size_C2end)

# Bind together and rename columns
sdlg_size_C2end<-cbind(sdlg_size_C2end,as.data.frame(ht.cohort2$ht.4), as.data.frame(ht.cohort2$trt), as.data.frame(ht.cohort2$block),as.data.frame(ht.cohort2$cohort),as.data.frame(ht.cohort2$days))
names(sdlg_size_C2end)[1] <- "LA_final"
names(sdlg_size_C2end)[2] <- "HT_final"
names(sdlg_size_C2end)[3] <- "trt"
names(sdlg_size_C2end)[4] <- "block"
names(sdlg_size_C2end)[5] <- "cohort"
names(sdlg_size_C2end)[6] <- "days"
na.omit(sdlg_size_C2end)
sdlg_size_end$block<-as.factor(sdlg_size_end$block) #set factor as a block
str(sdlg_size_C2end)
summary(sdlg_size_C2end)


# ######################################################
# Analysis: Is initial Leaf Area same accross treatments? Cohort 1
# ######################################################
# Visualizations of the data and ANOVA
boxplot(LAt1~trt,data=COHORT1.wide) #Cohort 1 LA initial by trt
boxplot(LAt1~block,data=COHORT1.wide) #Cohort 1 LA initial by trt
aov.LAt1C1<-aov(LAt1 ~ trt+block/trt, data = COHORT1.wide)
summary(aov.LAt1C1) # INITIAL LEAF AREA - NO SINIFICANT DIFFERENCE AMONG TREATMENTS, ONLY AMONG BLOCKS.1)

# Calclulate the mean & SD of initial LA by treatment 
na.omit(COHORT1.wide)%>% group_by(trt) %>%summarise(avg=mean(LAt1))
na.omit(COHORT1.wide)%>% group_by(trt) %>%summarise(sd=sd(LAt1))

# ######################################################
# Analysis: Is initial Leaf Area same accross treatments? Cohort 2
# ######################################################
# Visualizations of the data and ANOVA
boxplot(LAt1~trt,data=COHORT2.wide) #Cohort 1 LA initial by trt
boxplot(LAt1~block,data=COHORT2.wide) #Cohort 1 LA initial by trt
aov.LAt1C2<-aov(LAt1 ~ trt+block/trt, data = COHORT2.wide)
summary(aov.LAt1C2) # INITIAL LEAF AREA - NO SINIFICANT DIFFERENCE AMONG TREATMENTS, ONLY AMONG BLOCKS.

# Calclulate the mean & SD of initial LA by treatment 
na.omit(COHORT2.wide)%>% group_by(trt) %>%summarise(avg=mean(LAt1))
na.omit(COHORT2.wide)%>% group_by(trt) %>%summarise(sd=sd(LAt1))

# ######################################################
# Analysis: total leaf area after 12 months (Cohorts 1 and 2)
# ######################################################

ONE_YR_LA_C1<-dplyr::select(COHORT1.wide,cohort, block, trt, sdlg.type, LAt1,LAt6, days)
names(ONE_YR_LA_C1)[5] <- "LA_initial"
names(ONE_YR_LA_C1)[6] <- "LA_final"

ONE_YR_LA_C2<-dplyr::select(COHORT2.wide,cohort, block, trt, sdlg.type, LAt1,LAt4, days)
names(ONE_YR_LA_C2)[5] <- "LA_initial"
names(ONE_YR_LA_C2)[6] <- "LA_final"
ONE_YEAR_LA<-rbind(ONE_YR_LA_C1,ONE_YR_LA_C2)
ONE_YEAR_LA$cohort<-as.factor(ONE_YEAR_LA$cohort)
ONE_YEAR_LA$block<-as.factor(ONE_YEAR_LA$block)
summary(ONE_YEAR_LA)
str(ONE_YEAR_LA)

# Visualization
la1yr_plot<-ggplot(ONE_YEAR_LA, aes(x=LA_initial, y=LA_final, color=trt))+
  geom_point(shape=16, position=position_jitter(width=1,height=.5), size=3)+ #shape 1 = hollow circles, jitter plot to seperate overlapping points
  geom_smooth(method=lm, se=FALSE)+    #Add linear regression lines and Don't add shaded confidence region
  ylab("Final Leaf Area (cm2)") +
  xlab("Final leaf Area (cm2)")+
  ggtitle("A) 12 months")+
  scale_colour_manual(values=c("gray54", "orangered2","darkblue"))

la1yr_plot<-la1yr_plot+ theme_classic()+theme(legend.direction = 'vertical', 
                                              legend.position = c(0.1,0.85),
                                              plot.title = element_text(color="black", size=18, hjust=0.05, vjust=-.2),
                                              legend.background = element_rect(size=0.5, linetype="solid", color="black"),
                                              axis.title.x=element_text(colour="black", size = 18, vjust=-0.5),            #sets x axis title size, style, distance from axis #add , face = "bold" if you want bold
                                              axis.title.y=element_text(colour="black", size = 18, vjust=2),            #sets y axis title size, style, distance from axis #add , face = "bold" if you want bold
                                              axis.text=element_text(colour="black", size = 14),                             #sets size and style of labels on axes
                                              legend.title = element_blank(), #remove title of legend
                                              legend.text = element_text(color="black", size=14))   #plot margin - top, right, bottom, left
    

print(la1yr_plot)




boxplot(LA_final~cohort,data=ONE_YEAR_LA) #Cohort 1 LA initial by trt

# Calclulate means and SD of initial and final by trt
na.omit(ONE_YEAR_LA)%>% summarise(avg=mean(LA_final)) #leaf area 
na.omit(ONE_YEAR_LA)%>% summarise(avg=mean(LA_initial))
na.omit(ONE_YEAR_LA)%>% summarise(sd=sd(LA_final)) #leaf area 
na.omit(ONE_YEAR_LA)%>% summarise(sd=sd(LA_initial))

na.omit(ONE_YEAR_LA)%>% group_by(cohort) %>%summarise(avg=mean(LA_final))
na.omit(ONE_YEAR_LA)%>% group_by(cohort) %>%summarise(sd=sd(LA_final))

# GLMs
GLM_LA1.1<-glm(LA_final ~ 1, data = ONE_YEAR_LA, family = gaussian) #INTERCEPT
summary(GLM_LA1.1)
GLM_LA1.2<-glm(LA_final ~ LA_initial, data = ONE_YEAR_LA, family = gaussian) #INITIAL LEAF AREA 
summary(GLM_LA1.2)
GLM_LA1.3<-glm(LA_final ~ block, data = ONE_YEAR_LA, family = gaussian) #BLOCK
summary(GLM_LA1.3)
GLM_LA1.4<-glm(LA_final ~ trt, data = ONE_YEAR_LA, family = gaussian) #TREATMENT
summary(GLM_LA1.4)
GLM_LA1.5<-glm(LA_final ~ days, data = ONE_YEAR_LA, family = gaussian) #days
summary(GLM_LA1.5)
GLM_LA1.6<-glm(LA_final ~ cohort, data = ONE_YEAR_LA, family = gaussian) #cohort
summary(GLM_LA1.6)
GLM_LA1.7<-glm(LA_final ~ block+LA_initial, data = (ONE_YEAR_LA), family = gaussian) #INITIAL LA + BLOCK
summary(GLM_LA1.7)
GLM_LA1.8<-glm(LA_final ~ trt+LA_initial, data = ONE_YEAR_LA, family = gaussian) # INITIAL LA + TRT
summary(GLM_LA1.8)
GLM_LA1.9<-glm(LA_final ~ trt+cohort, data = ONE_YEAR_LA, family = gaussian) # INITIAL LA + TRT
summary(GLM_LA1.9)
GLM_LA1.10<-glm(LA_final ~ trt*cohort, data = ONE_YEAR_LA, family = gaussian) #TRT*cohort
summary(GLM_LA1.10)
GLM_LA1.11<-glm(LA_final ~ trt+LA_initial+block, data = ONE_YEAR_LA, family = gaussian) #INITIAL LA+TRT+BLOCK
summary(GLM_LA1.11)
GLM_LA1.12<-glm(LA_final ~ trt+LA_initial+cohort, data = ONE_YEAR_LA, family = gaussian) #INITIAL LA+TRT+cohort
summary(GLM_LA1.12)
GLM_LA1.13<-glm(LA_final ~ trt*cohort+LA_initial, data = ONE_YEAR_LA, family = gaussian) #INITIAL LA+TRT*cohort
summary(GLM_LA1.13)
GLM_LA1.14<-glm(LA_final ~ trt+LA_initial+cohort+block, data = ONE_YEAR_LA, family = gaussian) #INITIAL LA+TRT+BLOCK+cohort
summary(GLM_LA1.14)
GLM_LA1.15<-glm(LA_final ~ trt*cohort+block, data = ONE_YEAR_LA, family = gaussian) #INITIAL LA+TRT*cohort
summary(GLM_LA1.15)
GLM_LA1.16<-glm(LA_final ~ trt*cohort+LA_initial+block, data = ONE_YEAR_LA, family = gaussian) #INITIAL LA+TRT*cohort
summary(GLM_LA1.16)

model.namesLA1 <- c("1 Intercept", "2 Initial Leaf Area", "3 Block", "4 Density",
                 "5 Days",
                 "6 Cohort",
                 "7 Block+Initial Leaf Area",
                 "8 Density+Initial Leaf Area",
                 "9 Density+Cohort",
                 "10 Density*Cohort",
                 "11 Density+Initial Leaf Area+Block",
                 "12 Density+Initial Leaf Area+Cohort",
                 "13 Density*Cohort+Initial Leaf Area",
                 "14 Density+Cohort+Initial Leaf Area+Block",
                 "15 Density*Cohort+Block",
                 "16 Density*Cohort+Initial Leaf Area+Block")
                 
#OLD GLMS# 
# GLM_LA1.1<-glm(LA_final ~ 1, data = ONE_YEAR_LA, family = gaussian)
# summary(GLM_LA1.1)
# GLM_LA1.2<-glm(LA_final ~ LA_initial, data = ONE_YEAR_LA, family = gaussian)
# summary(GLM_LA1.2)
# GLM_LA1.3<-glm(LA_final ~ trt, data = ONE_YEAR_LA, family = gaussian)
# summary(GLM_LA1.3)
# GLM_LA1.4<-glm(LA_final ~ cohort, data = ONE_YEAR_LA, family = gaussian)
# summary(GLM_LA1.4)
# GLM_LA1.5<-glm(LA_final ~ block, data = ONE_YEAR_LA, family = gaussian)
# summary(GLM_LA1.5)
# GLM_LA1.6<-glm(LA_final ~ LA_initial+trt, data = ONE_YEAR_LA, family = gaussian)
# summary(GLM_LA1.6)
# GLM_LA1.7<-glm(LA_final ~ LA_initial+block, data = ONE_YEAR_LA, family = gaussian)
# summary(GLM_LA1.7)
# GLM_LA1.8<-glm(LA_final ~ LA_initial+cohort, data = ONE_YEAR_LA, family = gaussian)
# summary(GLM_LA1.8)
# GLM_LA1.9<-glm(LA_final ~ LA_initial*cohort*trt+block, data = ONE_YEAR_LA, family = gaussian)
# summary(GLM_LA1.9)

anova(GLM_LA1.1,GLM_LA1.2, test = "Chisq")  #imp fit by adding LA t1  over just intercept
anova(GLM_LA1.2,GLM_LA1.3, test = "Chisq")  # imp fit by adding trt over just intercept 
anova(GLM_LA1.3,GLM_LA1.4, test = "Chisq")  # imp fit by adding cohort instead of intercept
anova(GLM_LA1.4,GLM_LA1.5, test = "Chisq")  #imp fit by adding block instead of intercept
anova(GLM_LA1.4,GLM_LA1.6, test = "Chisq")  
anova(GLM_LA1.7,GLM_LA1.2, test = "Chisq")  
anova(GLM_LA1.2,GLM_LA1.7, test = "Chisq")  
anova(GLM_LA1.2,GLM_LA1.9, test = "Chisq")  

AIC(GLM_LA1.1,GLM_LA1.2,GLM_LA1.3,GLM_LA1.4,GLM_LA1.5,GLM_LA1.6,GLM_LA1.7,GLM_LA1.8,GLM_LA1.9,GLM_LA1.10,
    GLM_LA1.11,GLM_LA1.12,GLM_LA1.13,GLM_LA1.14,GLM_LA1.15,GLM_LA1.16) #for cohort 2 best model fit is just ht

# TO REPORT: http://www.ashander.info/posts/2015/10/model-selection-glms-aic-what-to-report/
summ.tableLA12 <- do.call(rbind, lapply(list(GLM_LA1.1,GLM_LA1.2,GLM_LA1.3,GLM_LA1.4,GLM_LA1.5,GLM_LA1.6,
                                             GLM_LA1.7,GLM_LA1.8,GLM_LA1.9,GLM_LA1.10,GLM_LA1.11,GLM_LA1.12,
                                             GLM_LA1.13,GLM_LA1.14,GLM_LA1.15,GLM_LA1.16), broom::glance))
summ.tableLA12
table.cols12 <- c("df.residual", "deviance", "AIC")
reported.tableLA12 <- summ.tableLA12[table.cols]
names(reported.tableLA12) <- c("Resid. Df", "Resid. Dev", "AIC")

reported.tableLA12[['dAIC']] <-  with(reported.tableLA12, AIC - min(AIC))
reported.tableLA12[['weight']] <- with(reported.tableLA12, exp(- 0.5 * dAIC) / sum(exp(- 0.5 * dAIC)))
reported.tableLA12$AIC <- NULL
reported.tableLA12$weight <- round(reported.tableLA12$weight, 2)
reported.tableLA12$dAIC <- round(reported.tableLA12$dAIC, 1)
row.names(reported.tableLA12) <- model.namesLA1
reported.tableLA12

write.table(reported.tableLA12, file = "LA_1yr_GLM.csv", , sep = ",", na = "NA", row.names = TRUE)


# ######################################################
# total leaf area after 24 months (Cohort 1)
# ######################################################
# 

# Are final height and length correlated?

plot(sdlg_size_end)
hist(sdlg_size_end$LA_final)
hist(sdlg_size_end$HT_final)
str(sdlg_size_end)
summary(sdlg_size_end)
cor(sdlg_size_end[2:3], method="spearman", use="complete.obs")  #Storng correlation between final leaf area and final height

# Calclulate LA and SD by trt
na.omit(sdlg_size_end)%>% group_by(trt) %>%summarise(avg=mean(LA_final)) #final seedling la by trt
na.omit(sdlg_size_end)%>% group_by(trt) %>%summarise(sd=sd(LA_final))

# do same but only seedlings aive the who experiment to see if the smaller size is due to density or transplanting in smaller ones when some died.
all.exp.sdlgs<-filter(sdlg_size_end, days >=773)
na.omit(all.exp.sdlgs)%>% group_by(trt) %>%summarise(avg=mean(LA_final)) #final seedling la by trt
na.omit(all.exp.sdlgs)%>% group_by(trt) %>%summarise(sd=sd(LA_final))
boxplot(LA_final~trt,data=all.exp.sdlgs) 

# Visualizations

la2yr_plot<-ggplot(sdlg_size_end, aes(x=LA_initial, y=LA_final, color=trt))+
  geom_point(shape=16, position=position_jitter(width=1,height=.5), size=3)+ #shape 1 = hollow circles, jitter plot to seperate overlapping points
  geom_smooth(method=lm, se=FALSE)+    #Add linear regression lines and Don't add shaded confidence region
  ylab("Final Leaf Area (cm2)") +
  xlab("Final leaf Area (cm2)")+
  ggtitle("B) 24 months")+
  scale_colour_manual(values=c("gray54", "orangered2","darkblue"))

la2yr_plot<-la2yr_plot+ theme_classic()+theme(legend.direction = 'vertical', 
                                              legend.position = c(0.1,0.85),
                                              plot.title = element_text(color="black", size=18, hjust=0.05, vjust=-.2),
                                              legend.background = element_rect(size=0.5, linetype="solid", color="black"),
                                              axis.title.x=element_text(colour="black", size = 18, vjust=-0.5),            #sets x axis title size, style, distance from axis #add , face = "bold" if you want bold
                                              axis.title.y=element_text(colour="black", size = 18, vjust=2),            #sets y axis title size, style, distance from axis #add , face = "bold" if you want bold
                                              axis.text=element_text(colour="black", size = 14),                             #sets size and style of labels on axes
                                              legend.title = element_blank(), #remove title of legend
                                              legend.text = element_text(color="black", size=14))   #plot margin - top, right, bottom, left


print(la2yr_plot)



ggplot(sdlg_size_end, aes(x=days, y=LA_final, color=trt)) + geom_point(shape=1)+
  geom_point(shape=1, position=position_jitter(width=1,height=.5))+ #shape 1 = hollow circles, jitter plot to seperate overlapping points
  geom_smooth(method=lm, se=FALSE)    #Add linear regression lines and Don't add shaded confidence region

boxplot(LA_final~trt,data=sdlg_size_end) #Cohort 1 LA initial by trt
boxplot(LA_final~block,data=sdlg_size_end) #Cohort 1 LA initial by block
boxplot(LA_final~days,data=sdlg_size_end) #Cohort 1 LA initial by block

hist(sdlg_size_end$LA_final)
hist(sdlg_size_end$LA_final)

GLM_LA2.1<-glm(LA_final ~ 1, data = sdlg_size_end, family = gaussian) #INTERCEPT
  summary(GLM_LA2.1)
GLM_LA2.2<-glm(LA_final ~ LA_initial, data = sdlg_size_end, family = gaussian) #INITIAL LEAF AREA 
  summary(GLM_LA2.2)
GLM_LA2.3<-glm(LA_final ~ block, data = sdlg_size_end, family = gaussian) #BLOCK
  summary(GLM_LA2.3)
GLM_LA2.4<-glm(LA_final ~ trt, data = sdlg_size_end, family = gaussian) #TREATMENT
  summary(GLM_LA2.4)
GLM_LA2.5<-glm(LA_final ~ days, data = sdlg_size_end, family = gaussian) #days
  summary(GLM_LA2.5)
GLM_LA2.6<-glm(LA_final ~ block+LA_initial, data = (sdlg_size_end), family = gaussian) #INITIAL LA + BLOCK
  summary(GLM_LA2.6)
GLM_LA2.7<-glm(LA_final ~ trt+LA_initial, data = sdlg_size_end, family = gaussian) # INITIAL LA + TRT
  summary(GLM_LA2.7)
GLM_LA2.8<-glm(LA_final ~ trt*days, data = sdlg_size_end, family = gaussian) #TRT*DAYS
  summary(GLM_LA2.8)
GLM_LA2.9<-glm(LA_final ~ trt+LA_initial+block, data = sdlg_size_end, family = gaussian) #INITIAL LA+TRT+BLOCK
  summary(GLM_LA2.9)
GLM_LA2.10<-glm(LA_final ~ trt+LA_initial+days, data = sdlg_size_end, family = gaussian) #INITIAL LA+TRT+DAYS
  summary(GLM_LA2.10)
GLM_LA2.11<-glm(LA_final ~ trt*days+LA_initial, data = sdlg_size_end, family = gaussian) #INITIAL LA+TRT*DAYS
  summary(GLM_LA2.11)
GLM_LA2.12<-glm(LA_final ~ trt+LA_initial+days+block, data = sdlg_size_end, family = gaussian) #INITIAL LA+TRT+BLOCK+DAYS
  summary(GLM_LA2.12)
GLM_LA2.13<-glm(LA_final ~ trt*days+block, data = sdlg_size_end, family = gaussian) #INITIAL LA+TRT*DAYS
  summary(GLM_LA2.13)
GLM_LA2.14<-glm(LA_final ~ trt*days+LA_initial+block, data = sdlg_size_end, family = gaussian) #INITIAL LA+TRT*DAYS
  summary(GLM_LA2.14)
GLM_LA2.15<-glm(LA_final ~ trt+days, data = sdlg_size_end, family = gaussian) #TRT*DAYS
  summary(GLM_LA2.15)


model.names <- c("1 Intercept", "2 Initial Leaf Area", "3 Block", "4 Density","5 Days Since Transplant",
                 "6 Initial Leaf Area+Block",
                 "7 Density+Initial Leaf Area",
                 "8 Density*Days Since Transplant",
                 "9 Density+Initial Leaf Area+Block",
                 "10 Density+Initial Leaf Area+Days Since Transplant",
                 "11 Density*Days Since Transplant+Initial Leaf Area",
                 "12 Density+Days Since Transplant+Initial Leaf Area+Block",
                 "13 Density*Days Since Transplant+Block",
                 "14 Density*Days Since Transplant+Initial Leaf Area+Block",
                 "15 Density+Days Since Transplanting")


AIC(GLM_LA2.1,GLM_LA2.2,GLM_LA2.3,GLM_LA2.4,GLM_LA2.5,GLM_LA2.6,GLM_LA2.7,GLM_LA2.8, GLM_LA2.9, 
    GLM_LA2.10, GLM_LA2.11, GLM_LA2.12, GLM_LA2.13, GLM_LA2.14, GLM_LA2.15) #adding block doesn't result in lower AIC


# TO REPORT: http://www.ashander.info/posts/2015/10/model-selection-glms-aic-what-to-report/
summ.tableLA24 <- do.call(rbind, lapply(list(GLM_LA2.1,GLM_LA2.2,GLM_LA2.3,GLM_LA2.4,GLM_LA2.5,GLM_LA2.6,GLM_LA2.7,
                                             GLM_LA2.8, GLM_LA2.9, GLM_LA2.10, GLM_LA2.11, GLM_LA2.12, GLM_LA2.13, GLM_LA2.14, GLM_LA2.15), broom::glance))
summ.tableLA24
table.colsLA24 <- c("df.residual", "deviance", "AIC")
reported.tableLA24 <- summ.tableLA24[table.colsLA24]
names(reported.tableLA24) <- c("Resid. Df", "Resid. Dev", "AIC")

reported.tableLA24[['dAIC']] <-  with(reported.tableLA24, AIC - min(AIC))
reported.tableLA24[['weight']] <- with(reported.tableLA24, exp(- 0.5 * dAIC) / sum(exp(- 0.5 * dAIC)))
reported.tableLA24$AIC <- NULL
reported.tableLA24$weight <- round(reported.tableLA24$weight, 2)
reported.tableLA24$dAIC <- round(reported.tableLA24$dAIC, 1)
row.names(reported.tableLA24) <- model.names
reported.tableLA24

write.table(reported.tableLA24, file = "LA_2yr_GLM.csv", , sep = ",", na = "NA", row.names = TRUE)

anova(GLM_LA2.1,GLM_LA2.2, test = "Chisq")  #do not imp fit by adding INITIAL LA over just intercept
anova(GLM_LA2.1,GLM_LA2.3, test = "Chisq")  #imp fit by adding BLOCK over just intercept 
anova(GLM_LA2.1,GLM_LA2.4, test = "Chisq")  #imp fit by adding TRT P=0.056
anova(GLM_LA2.1,GLM_LA2.5, test = "Chisq")  #adding DAYS SINCE TPLANT imporives fit
anova(GLM_LA2.3,GLM_LA2.6, test = "Chisq")  #adding trt and trt*ht interaction best
anova(GLM_LA2.10,GLM_LA2.12, test = "Chisq")  #adding trt and trt*ht interaction and days
# ######################################################
# ######################################################
# BIOMASS ANALYSES - SET UP
# ######################################################
# ####################################################### ######################################################

bmass[, "RSratio"] <-(bmass$root.biomass/(bmass$lf.biomass+bmass$stem.biomass))
bmass[, "Total.bmass"] <-bmass$root.biomass+bmass$lf.biomass+bmass$stem.biomass
bmass<-full_join(bmass, canopy, by = "block")
bmass$location<-NULL
bmass$site.openess.percent<-NULL
bmass$mask.openess.percent<-NULL
bmass$sky.area.percent<-NULL
bmass$LAI.4.ring<-NULL 
bmass$LAI.5.ring<-NULL
bmass$cohort<-as.factor(bmass$cohort)
bmass$block<-as.integer(bmass$block)

bmass<-arrange(bmass,cohort,block,trt)

LA_for_Bmass<-dplyr::filter(ONE_YEAR_LA,sdlg.type=="focal")
LA_for_Bmass<-droplevels(LA_for_Bmass)
LA_for_Bmass<-arrange(LA_for_Bmass,cohort,block,trt)

str(LA_for_Bmass)
str(bmass)

bmass<-cbind(bmass,LA_for_Bmass$LA_initial,LA_for_Bmass$LA_final)

names(bmass)[11] <- "LA_initial"
names(bmass)[12] <- "LA_final"

bmassC1<-filter(bmass,  cohort == "1") #cohort 1 after 2 years
bmassC1$cohort<-as.factor(bmassC1$cohort)

bmassC2<-filter(bmass,  cohort == "2") #cohort 2 after 1 year
bmassC2$cohort<-as.factor(bmassC2$cohort)

# Corr of leaf area and biomass
corr1<-bmass$Total.bmass
corr2<-bmass$LA_final
CORR<-cbind(corr1,corr2)
cor(CORR, method="spearman", use="complete.obs")  #Storng correlation between final leaf area and final height

# Visualizations

ggplot(bmass, aes(x=LA_initial, y=Total.bmass, color=trt)) + geom_point(shape=1)+
  geom_point(shape=1, position=position_jitter(width=1,height=.5))+ #shape 1 = hollow circles, jitter plot to seperate overlapping points
  geom_smooth(method=lm, se=FALSE)    #Add linear regression lines and Don't add shaded confidence region


# ######################################################
# total biomass after 12 months (Cohort 2)
# ######################################################
# 
# Histogram
hist(bmassC2$Total.bmass)
# BOX PLOTS
boxplot(Total.bmass~trt,data=bmassC2) #Cohort bmass final

# GLMs
GLM_BM_1<-glm(Total.bmass ~ 1, data = bmassC2, family = gaussian)
summary(GLM_BM_1)
GLM_BM_2<-glm(Total.bmass ~ LA_initial, data = bmassC2, family = gaussian)
summary(GLM_BM_2)
GLM_BM_3<-glm(Total.bmass ~ trt, data = bmassC2, family = gaussian)
summary(GLM_BM_3)
GLM_BM_4<-glm(Total.bmass ~ block, data = bmassC2, family = gaussian)
summary(GLM_BM_4)
GLM_BM_5<-glm(Total.bmass ~ trt+block, data = bmassC2, family = gaussian)
summary(GLM_BM_5)
GLM_BM_6<-glm(Total.bmass ~ trt+LA_initial+block, data = bmassC2, family = gaussian)
summary(GLM_BM_6)

anova(GLM_BM_1,GLM_BM_1, test = "Chisq")  #imp fit by adding HT over just intercept
anova(GLM_BM_1,GLM_BM_1, test = "Chisq")  #no imp fit by adding trt over just intercept (but close)
anova(GLM_BM_1,GLM_BM_1, test = "Chisq")  #no imp fit by adding trt instead of ht
anova(GLM_BM_1,GLM_BM_1, test = "Chisq")  #adding trt and ht better than just ht 
anova(GLM_BM_1,GLM_BM_1, test = "Chisq")  #adding trt and trt*ht interaction best
anova(GLM_BM_1,GLM_BM_1, test = "Chisq")  #adding trt and trt*ht interaction and days

AIC(GLM_BM_1,GLM_BM_2,GLM_BM_3,GLM_BM_4,GLM_BM_5,GLM_BM_6) #adding block doesn't result in lower AIC

# TO REPORT: http://www.ashander.info/posts/2015/10/model-selection-glms-aic-what-to-report/



bm1yr_plot<-ggplot(bmassC2, aes(x=LA_final, y=Total.bmass, color=trt))+
  geom_point(shape=16, position=position_jitter(width=1,height=.5), size=3)+ #shape 1 = hollow circles, jitter plot to seperate overlapping points
  geom_smooth(method=lm, se=FALSE)+    #Add linear regression lines and Don't add shaded confidence region
  ylab("Total biomass (g)") +
  xlab("Final leaf Area (cm2)")+
  ggtitle("A) 12 months")+
  scale_colour_manual(values=c("gray54", "orangered2","darkblue"))

bm1yr_plot<-bm1yr_plot+ theme_classic()+theme(legend.direction = 'vertical', 
                                              legend.position = c(0.1,0.85),
                                              plot.title = element_text(color="black", size=18, hjust=0.05, vjust=-.2),
                                              legend.background = element_rect(size=0.5, linetype="solid", color="black"),
                                              axis.title.x=element_text(colour="black", size = 18, vjust=-0.5),            #sets x axis title size, style, distance from axis #add , face = "bold" if you want bold
                                              axis.title.y=element_text(colour="black", size = 18, vjust=2),            #sets y axis title size, style, distance from axis #add , face = "bold" if you want bold
                                              axis.text=element_text(colour="black", size = 14),                             #sets size and style of labels on axes
                                              legend.title = element_blank(), #remove title of legend
                                              legend.text = element_text(color="black", size=14))   #plot margin - top, right, bottom, left


print(bm1yr_plot)



# ######################################################
# total biomass after 24 months (Cohort 1)
# ######################################################
# 



bm2yr_plot<-ggplot(bmassC1, aes(x=LA_final, y=Total.bmass, color=trt))+
  geom_point(shape=16, position=position_jitter(width=1,height=.5), size=3)+ #shape 1 = hollow circles, jitter plot to seperate overlapping points
  geom_smooth(method=lm, se=FALSE)+    #Add linear regression lines and Don't add shaded confidence region
  ylab("Total biomass (g)") +
  xlab("Final leaf Area (cm2)")+
  ggtitle("B) 24 months")+
  scale_colour_manual(values=c("gray54", "orangered2","darkblue"))

bm2yr_plot<-bm2yr_plot+ theme_classic()+theme(legend.direction = 'vertical', 
                                              legend.position = c(0.1,0.85),
                                              plot.title = element_text(color="black", size=18, hjust=0.05, vjust=-.2),
                                              legend.background = element_rect(size=0.5, linetype="solid", color="black"),
                                              axis.title.x=element_text(colour="black", size = 18, vjust=-0.5),            #sets x axis title size, style, distance from axis #add , face = "bold" if you want bold
                                              axis.title.y=element_text(colour="black", size = 18, vjust=2),            #sets y axis title size, style, distance from axis #add , face = "bold" if you want bold
                                              axis.text=element_text(colour="black", size = 14),                             #sets size and style of labels on axes
                                              legend.title = element_blank(), #remove title of legend
                                              legend.text = element_text(color="black", size=14))   #plot margin - top, right, bottom, left


print(bm2yr_plot)





# ######################################################
# RS after 12 months (Cohort 2)
# ######################################################
# 

rs1yr_plot<-ggplot(bmassC2, aes(x=LA_final, y=RSratio, color=trt))+
  geom_point(shape=16, position=position_jitter(width=1,height=.5), size=3)+ #shape 1 = hollow circles, jitter plot to seperate overlapping points
  geom_smooth(method=lm, se=FALSE)+    #Add linear regression lines and Don't add shaded confidence region
  ylab("root:shoot ratio") +
  xlab("Final leaf Area (cm2)")+
  ggtitle("C) 12 months")+
  scale_colour_manual(values=c("gray54", "orangered2","darkblue"))

rs1yr_plot<-rs1yr_plot+ theme_classic()+theme(legend.direction = 'vertical', 
                                              legend.position = c(0.1,0.85),
                                              plot.title = element_text(color="black", size=18, hjust=0.05, vjust=-.2),
                                              legend.background = element_rect(size=0.5, linetype="solid", color="black"),
                                              axis.title.x=element_text(colour="black", size = 18, vjust=-0.5),            #sets x axis title size, style, distance from axis #add , face = "bold" if you want bold
                                              axis.title.y=element_text(colour="black", size = 18, vjust=2),            #sets y axis title size, style, distance from axis #add , face = "bold" if you want bold
                                              axis.text=element_text(colour="black", size = 14),                             #sets size and style of labels on axes
                                              legend.title = element_blank(), #remove title of legend
                                              legend.text = element_text(color="black", size=14))   #plot margin - top, right, bottom, left


print(rs1yr_plot)











# ######################################################
# RS after 24 months (Cohort 1)
# ######################################################
# 



rs2yr_plot<-ggplot(bmassC1, aes(x=LA_final, y=RSratio, color=trt))+
  geom_point(shape=16, position=position_jitter(width=1,height=.5), size=3)+ #shape 1 = hollow circles, jitter plot to seperate overlapping points
  geom_smooth(method=lm, se=FALSE)+    #Add linear regression lines and Don't add shaded confidence region
  ylab("root:shoot ratio") +
  xlab("Final leaf Area (cm2)")+
  ggtitle("D) 24 months")+
  scale_colour_manual(values=c("gray54", "orangered2","darkblue"))

rs2yr_plot<-rs2yr_plot+ theme_classic()+theme(legend.direction = 'vertical', 
                                              legend.position = c(0.1,0.85),
                                              plot.title = element_text(color="black", size=18, hjust=0.05, vjust=-.2),
                                              legend.background = element_rect(size=0.5, linetype="solid", color="black"),
                                              axis.title.x=element_text(colour="black", size = 18, vjust=-0.5),            #sets x axis title size, style, distance from axis #add , face = "bold" if you want bold
                                              axis.title.y=element_text(colour="black", size = 18, vjust=2),            #sets y axis title size, style, distance from axis #add , face = "bold" if you want bold
                                              axis.text=element_text(colour="black", size = 14),                             #sets size and style of labels on axes
                                              legend.title = element_blank(), #remove title of legend
                                              legend.text = element_text(color="black", size=14))   #plot margin - top, right, bottom, left


print(rs2yr_plot)























############


# SEEDLINGS FROM THE DEMOG PLOTS
# #NEED TO IDENTIFY WHICH IS THE ONE THAT IS THE FOCAL ONE IN THE PAIR!!!


plot.sdlgs<-read.csv("demog_plot_sdlgs.csv", dec=".", header = TRUE, sep = ",", check.names=FALSE )
plot.sdlgs[, "rgr.ht.0809"] <-(log(plot.sdlgs$ht.2009)-log(plot.sdlgs$ht.2008))/365
plot.sdlgs[, "rgr.ht.0810"] <-(log(plot.sdlgs$ht.2010)-log(plot.sdlgs$ht.2008))/365

boxplot(rgr.ht.0810~trt,data=plot.sdlgs) #08 sdeedlings rgr 0810
boxplot(ht.2008~trt,data=plot.sdlgs) #08 sflgs ht 08
boxplot(ht.2009~trt,data=plot.sdlgs) #08 sflgs ht 09
boxplot(ht.2010~trt,data=plot.sdlgs) #08 sflgs final ht 2010



aov.plot.0810<-aov(rgr.ht.0810 ~ trt+light+ht.2008, data = plot.sdlgs)
summary(aov.plot.0810)



plot(plot.sdlgs$ht.2008,plot.sdlgs$rgr.ht.0810)
plot(plot.sdlgs$ht.2008,plot.sdlgs$ht.2010)


ggplot(plot.sdlgs, aes(x=ht.2008, y=ht.2010, color=trt)) + geom_point(shape=1)+
  geom_point(shape=1, position=position_jitter(width=1,height=.5))+ #shape 1 = hollow circles, jitter plot to seperate overlapping points
  geom_smooth(method=lm, se=FALSE)    #Add linear regression lines and Don't add shaded confidence region



