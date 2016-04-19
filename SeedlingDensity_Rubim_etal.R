
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
                    cohort1.la$trt, cohort1.la$sdlg.type,cohort1.la$seedling.id.no,  cohort1.la$interval,cohort1.la$days), 
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
COHORT1.wide<- dcast(COHORT1, cohort + block + trt + sdlg.type + sdlg.id.no + days ~ interval)
COHORT2.wide<- dcast(COHORT2, cohort + block + trt + sdlg.type + sdlg.id.no + days ~ interval)


#Sigh. Have to rename the columes of the intervals, it's just easier than remembering to use "1" in all formulas
#renaming columsn cohort 1
names(COHORT1.wide)[7] <- "LAt1"
names(COHORT1.wide)[8] <- "LAt2"
names(COHORT1.wide)[9] <- "LAt3"
names(COHORT1.wide)[10] <- "LAt4"
names(COHORT1.wide)[11] <- "LAt5"
names(COHORT1.wide)[12] <- "LAt6"
names(COHORT1.wide)[13] <- "LAt7"
names(COHORT1.wide)[14] <- "LAt8"
names(COHORT1.wide)[15] <- "LAt9"
#renaming columns cohort 2
names(COHORT2.wide)[7] <- "LAt1"
names(COHORT2.wide)[8] <- "LAt2"
names(COHORT2.wide)[9] <- "LAt3"
names(COHORT2.wide)[10] <- "LAt4"


#A little sorting, just to make it easier to visualize
COHORT1.wide<-COHORT1.wide[with(COHORT1.wide, order(sdlg.id.no,sdlg.type, block, trt)), ]
COHORT2.wide<-COHORT2.wide[with(COHORT2.wide, order(sdlg.id.no,sdlg.type, block, trt)), ]

#################
###CALCLULATIONS OF RGR BAsed on Leaf Area
##################
#Formula for RGR is rgr=[ln(W2)-ln(W1)]/t2-t1, where W is weegith at time t 

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

# Box plot of median leaf area of plants in each block at T1 
boxplot(LAt1~block,data=COHORT1.wide) #Cohort 1
boxplot(LAt1~block,data=COHORT2.wide) #Cohort2

# hist(COHORT1.wide$rgrLA_2yrs)
# hist(COHORT1.wide$rgrLA_yr1)
# hist(COHORT2.wide$rgrLA_yr2)

#BOX PLOT INCLUDING ALL PLANTS 

# rgrALL1.6 <- ggplot(COHORT1.wide, aes(x=trt, y=rgrLA_yr1)) + 
#   geom_boxplot()
# rgrALL1.6
# 
# rgrALL1.9 <- ggplot(COHORT1.wide, aes(x=trt, y=rgrLA_2yrs)) + 
#   geom_boxplot()
# rgrALL1.9
# 
# rgrALL2 <- ggplot(COHORT2.wide, aes(x=trt, y=rgrLA_yr2)) + 
#   geom_boxplot()
# rgrALL2


####### PUTTING BOTH COHORTS TOGETHER - RGR #######
# Cohort 1
LFAREA_1_RGR<-COHORT1.wide
LFAREA_1_RGR <- LFAREA_1_RGR[ -c(7:15)]
LFAREA_1_RGR<-gather(LFAREA_1_RGR, "interval", "rgr.la", 7:14)  

# Cohort 2
LFAREA_2_RGR<-COHORT2.wide
LFAREA_2_RGR <- LFAREA_2_RGR[ -c(7:10)]
LFAREA_2_RGR<-gather(LFAREA_2_RGR, "interval", "rgr.la", 7:10)  
# str(LFAREA_1_RGR)
# str(LFAREA_2_RGR)

#Bind the 2
LARGR_BOTH<-rbind(LFAREA_1_RGR,LFAREA_2_RGR)
LARGR_BOTH[, "season"] <-0
LARGR_BOTH[, "calendar_year"] <-0
LARGR_BOTH[, "cohort_year"] <-0
LARGR_BOTH[, "duration"] <-0

LARGR_BOTH$interval <- as.character(LARGR_BOTH$interval)
LARGR_BOTH <- within(LARGR_BOTH, season[interval == 'rgrLA_wet1' ] <- 'rainy')
LARGR_BOTH <- within(LARGR_BOTH, season[interval == 'rgrLA_wet2' ] <- 'rainy')
LARGR_BOTH <- within(LARGR_BOTH, season[interval == 'rgrLA_wet3' ] <- 'rainy')
LARGR_BOTH <- within(LARGR_BOTH, season[interval == 'rgrLA_dry1' ] <- 'dry')
LARGR_BOTH <- within(LARGR_BOTH, season[interval == 'rgrLA_dry2' ] <- 'dry')
LARGR_BOTH <- within(LARGR_BOTH, season[interval == 'rgrLA_yr1' ] <- '1R1D')
LARGR_BOTH <- within(LARGR_BOTH, season[interval == 'rgrLA_yr2' ] <- '2R1D')
LARGR_BOTH <- within(LARGR_BOTH, season[interval == 'rgrLA_yrs1-2' ] <- '3R2D')
LARGR_BOTH <- within(LARGR_BOTH, season[interval == 'rgrLA_1yr' ] <- '2R1D')

LARGR_BOTH <- within(LARGR_BOTH, calendar_year[interval == 'rgrLA_wet1' ] <- '2008')
LARGR_BOTH <- within(LARGR_BOTH, calendar_year[interval == 'rgrLA_dry1' ] <- '2008')
LARGR_BOTH <- within(LARGR_BOTH, calendar_year[interval == 'rgrLA_wet2' ] <- '2009')
LARGR_BOTH <- within(LARGR_BOTH, calendar_year[interval == 'rgrLA_dry2' ] <- '2009')
LARGR_BOTH <- within(LARGR_BOTH, calendar_year[interval == 'rgrLA_wet3' ] <- '2010')
LARGR_BOTH <- within(LARGR_BOTH, calendar_year[interval == 'rgrLA_yr1'  & cohort == 1] <- '2008')
LARGR_BOTH <- within(LARGR_BOTH, calendar_year[interval == 'rgrLA_yr2'  & cohort == 1] <- '2009-2010')
LARGR_BOTH <- within(LARGR_BOTH, calendar_year[interval == 'rgrLA_yrs1-2'] <- '2008-2010')
LARGR_BOTH <- within(LARGR_BOTH, calendar_year[interval == 'rgrLA_yr1'  & cohort == 2] <- '2009-2010')

LARGR_BOTH <- within(LARGR_BOTH, duration[interval == 'rgrLA_wet1'] <- "seasonal")
LARGR_BOTH <- within(LARGR_BOTH, duration[interval == 'rgrLA_wet1'] <- "seasonal")
LARGR_BOTH <- within(LARGR_BOTH, duration[interval == 'rgrLA_wet2'] <- "seasonal")
LARGR_BOTH <- within(LARGR_BOTH, duration[interval == 'rgrLA_wet3'] <- "seasonal")
LARGR_BOTH <- within(LARGR_BOTH, duration[interval == 'rgrLA_dry1'] <- "seasonal")
LARGR_BOTH <- within(LARGR_BOTH, duration[interval == 'rgrLA_dry2'] <- "seasonal")
LARGR_BOTH <- within(LARGR_BOTH, duration[interval == 'rgrLA_yr1' & cohort == 1  ] <- 'annual')
LARGR_BOTH <- within(LARGR_BOTH, duration[interval == 'rgrLA_yr2' & cohort == 1  ] <- 'annual')
LARGR_BOTH <- within(LARGR_BOTH, duration[interval == 'rgrLA_yr1' & cohort == 2  ] <- 'annual')
LARGR_BOTH <- within(LARGR_BOTH, duration[interval == 'rgrLA_yrs1-2' & cohort == 1] <- 'multiyear')

LARGR_BOTH <- within(LARGR_BOTH, cohort_year[interval == 'rgrLA_wet1' & cohort == 1  ] <- 'cohort_yr_1')
LARGR_BOTH <- within(LARGR_BOTH, cohort_year[interval == 'rgrLA_dry1' & cohort == 1  ] <- 'cohort_yr_1')
LARGR_BOTH <- within(LARGR_BOTH, cohort_year[interval == 'rgrLA_wet2' & cohort == 1  ] <- 'cohort_yr_2')
LARGR_BOTH <- within(LARGR_BOTH, cohort_year[interval == 'rgrLA_dry2' & cohort == 1  ] <- 'cohort_yr_2')
LARGR_BOTH <- within(LARGR_BOTH, cohort_year[interval == 'rgrLA_wet3' & cohort == 1  ] <- 'cohort_yr_2')
LARGR_BOTH <- within(LARGR_BOTH, cohort_year[interval == 'rgrLA_wet2' & cohort == 2  ] <- 'cohort_yr_1')
LARGR_BOTH <- within(LARGR_BOTH, cohort_year[interval == 'rgrLA_dry2' & cohort == 2  ] <- 'cohort_yr_1')
LARGR_BOTH <- within(LARGR_BOTH, cohort_year[interval == 'rgrLA_wet3' & cohort == 2  ] <- 'cohort_yr_1')
LARGR_BOTH <- within(LARGR_BOTH, cohort_year[interval == 'rgrLA_yr1' & cohort == 1  ] <- 'cohort_yr_1')
LARGR_BOTH <- within(LARGR_BOTH, cohort_year[interval == 'rgrLA_yr2' & cohort == 1] <- 'cohort_yr_2')
LARGR_BOTH <- within(LARGR_BOTH, cohort_year[interval == 'rgrLA_yr1' & cohort == 2] <- 'cohort_yr_1')
LARGR_BOTH <- within(LARGR_BOTH, cohort_year[interval == 'rgrLA_yrs1-2' & cohort == 1] <- 'cohort_yr_1+2')

# COnvert back to factors
LARGR_BOTH$interval <- as.factor(LARGR_BOTH$interval)
LARGR_BOTH$season <- as.factor(LARGR_BOTH$season)
LARGR_BOTH$calendar_year <- as.factor(LARGR_BOTH$calendar_year)
LARGR_BOTH$cohort_year <- as.factor(LARGR_BOTH$cohort_year)
LARGR_BOTH$duration <- as.factor(LARGR_BOTH$duration)
summary(LARGR_BOTH)
str(LARGR_BOTH)


#################
###CALCLULATIONS OF RGR BAsed on seedling height
##################
#rgr.ht=(log(Exp_Data$ht.final)-log(Exp_Data$ht.initial))/Exp_Data$days

# COHORT 1

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

HTRGR_BOTH<-rbind(ht.cohort1.1,ht.cohort2.1)
HTRGR_BOTH[, "season"] <-0
HTRGR_BOTH[, "calendar_year"] <-0
HTRGR_BOTH[, "cohort_year"] <-0
HTRGR_BOTH[, "duration"] <-0

HTRGR_BOTH$interval <- as.character(HTRGR_BOTH$interval)

HTRGR_BOTH <- within(HTRGR_BOTH, season[interval == 'rgrHT_wet1' ] <- 'rainy')
HTRGR_BOTH <- within(HTRGR_BOTH, season[interval == 'rgrHT_wet2' ] <- 'rainy')
HTRGR_BOTH <- within(HTRGR_BOTH, season[interval == 'rgrHT_wet3' ] <- 'rainy')
HTRGR_BOTH <- within(HTRGR_BOTH, season[interval == 'rgrHT_dry1' ] <- 'dry')
HTRGR_BOTH <- within(HTRGR_BOTH, season[interval == 'rgrHT_dry2' ] <- 'dry')
HTRGR_BOTH <- within(HTRGR_BOTH, season[interval == 'rgrHT_yr1' ] <- '1R1D')
HTRGR_BOTH <- within(HTRGR_BOTH, season[interval == 'rgrHT_yr2' ] <- '2R1D')
HTRGR_BOTH <- within(HTRGR_BOTH, season[interval == 'rgrHT_yrs1-2' ] <- '3R2D')
HTRGR_BOTH <- within(HTRGR_BOTH, season[interval == 'rgrHT_1yr' ] <- '2R1D')

HTRGR_BOTH <- within(HTRGR_BOTH, calendar_year[interval == 'rgrHT_wet1' ] <- '2008')
HTRGR_BOTH <- within(HTRGR_BOTH, calendar_year[interval == 'rgrHT_dry1' ] <- '2008')
HTRGR_BOTH <- within(HTRGR_BOTH, calendar_year[interval == 'rgrHT_wet2' ] <- '2009')
HTRGR_BOTH <- within(HTRGR_BOTH, calendar_year[interval == 'rgrHT_dry2' ] <- '2009')
HTRGR_BOTH <- within(HTRGR_BOTH, calendar_year[interval == 'rgrHT_wet3' ] <- '2010')
HTRGR_BOTH <- within(HTRGR_BOTH, calendar_year[interval == 'rgrHT_yr1'  & cohort == 1] <- '2008')
HTRGR_BOTH <- within(HTRGR_BOTH, calendar_year[interval == 'rgrHT_yr2'  & cohort == 1] <- '2009-2010')
HTRGR_BOTH <- within(HTRGR_BOTH, calendar_year[interval == 'rgrHT_yrs1-2'] <- '2008-2010')
HTRGR_BOTH <- within(HTRGR_BOTH, calendar_year[interval == 'rgrHT_yr1'  & cohort == 2] <- '2009-2010')

HTRGR_BOTH <- within(HTRGR_BOTH, duration[interval == 'rgrHT_wet1'] <- "seasonal")
HTRGR_BOTH <- within(HTRGR_BOTH, duration[interval == 'rgrHT_wet1'] <- "seasonal")
HTRGR_BOTH <- within(HTRGR_BOTH, duration[interval == 'rgrHT_wet2'] <- "seasonal")
HTRGR_BOTH <- within(HTRGR_BOTH, duration[interval == 'rgrHT_wet3'] <- "seasonal")
HTRGR_BOTH <- within(HTRGR_BOTH, duration[interval == 'rgrHT_dry1'] <- "seasonal")
HTRGR_BOTH <- within(HTRGR_BOTH, duration[interval == 'rgrHT_dry2'] <- "seasonal")
HTRGR_BOTH <- within(HTRGR_BOTH, duration[interval == 'rgrHT_yr1' & cohort == 1  ] <- 'annual')
HTRGR_BOTH <- within(HTRGR_BOTH, duration[interval == 'rgrHT_yr2' & cohort == 1  ] <- 'annual')
HTRGR_BOTH <- within(HTRGR_BOTH, duration[interval == 'rgrHT_yr1' & cohort == 2  ] <- 'annual')
HTRGR_BOTH <- within(HTRGR_BOTH, duration[interval == 'rgrHT_yrs1-2' & cohort == 1] <- 'multiyear')

HTRGR_BOTH <- within(HTRGR_BOTH, cohort_year[interval == 'rgrHT_wet1' & cohort == 1  ] <- 'cohort_yr_1')
HTRGR_BOTH <- within(HTRGR_BOTH, cohort_year[interval == 'rgrHT_dry1' & cohort == 1  ] <- 'cohort_yr_1')
HTRGR_BOTH <- within(HTRGR_BOTH, cohort_year[interval == 'rgrHT_wet2' & cohort == 1  ] <- 'cohort_yr_2')
HTRGR_BOTH <- within(HTRGR_BOTH, cohort_year[interval == 'rgrHT_dry2' & cohort == 1  ] <- 'cohort_yr_2')
HTRGR_BOTH <- within(HTRGR_BOTH, cohort_year[interval == 'rgrHT_wet3' & cohort == 1  ] <- 'cohort_yr_2')
HTRGR_BOTH <- within(HTRGR_BOTH, cohort_year[interval == 'rgrHT_wet2' & cohort == 2  ] <- 'cohort_yr_1')
HTRGR_BOTH <- within(HTRGR_BOTH, cohort_year[interval == 'rgrHT_dry2' & cohort == 2  ] <- 'cohort_yr_1')
HTRGR_BOTH <- within(HTRGR_BOTH, cohort_year[interval == 'rgrHT_wet3' & cohort == 2  ] <- 'cohort_yr_1')
HTRGR_BOTH <- within(HTRGR_BOTH, cohort_year[interval == 'rgrHT_yr1' & cohort == 1  ] <- 'cohort_yr_1')
HTRGR_BOTH <- within(HTRGR_BOTH, cohort_year[interval == 'rgrHT_yr2' & cohort == 1] <- 'cohort_yr_2')
HTRGR_BOTH <- within(HTRGR_BOTH, cohort_year[interval == 'rgrHT_yr1' & cohort == 2] <- 'cohort_yr_1')
HTRGR_BOTH <- within(HTRGR_BOTH, cohort_year[interval == 'rgrHT_yrs1-2' & cohort == 1] <- 'cohort_yr_1+2')

# Convert back to factors
HTRGR_BOTH$interval <- as.factor(HTRGR_BOTH$interval)
HTRGR_BOTH$season <- as.factor(HTRGR_BOTH$season)
HTRGR_BOTH$calendar_year <- as.factor(HTRGR_BOTH$calendar_year)
HTRGR_BOTH$cohort_year <- as.factor(HTRGR_BOTH$cohort_year)
HTRGR_BOTH$duration <- as.factor(HTRGR_BOTH$duration)
summary(HTRGR_BOTH)
str(HTRGR_BOTH)

HTRGR_BOTH<-HTRGR_BOTH[with(HTRGR_BOTH, order(block, trt, sdlg.type,seedling.id.no, decreasing = FALSE)), ]
LARGR_BOTH<-LARGR_BOTH[with(LARGR_BOTH, order(block, trt, sdlg.type,sdlg.id.no, decreasing = FALSE)), ]
str(HTRGR_BOTH)
# summary(ht.cohort2)
summary(LARGR_BOTH)

###################################################
######      ANALYSES  --- GROWTH
######################################################

# Are final height and length correlated?
COHORT1.wide<-COHORT1.wide[with(COHORT1.wide, order(sdlg.id.no, sdlg.type, block, trt)), ]
ht.cohort1<-ht.cohort1[with(ht.cohort1, order(seedling.id.no, sdlg.type, block, trt)), ]
sdlg_size_end<-as.data.frame(COHORT1.wide$LAt9)
sdlg_size_end<-cbind(sdlg_size_end,as.data.frame(ht.cohort1$ht.9),as.data.frame(ht.cohort1$trt),as.data.frame(ht.cohort1$block),as.data.frame(ht.cohort1$cohort))
names(sdlg_size_end)[1] <- "LA_final"
names(sdlg_size_end)[2] <- "HT_final"
names(sdlg_size_end)[3] <- "trt"
names(sdlg_size_end)[4] <- "block"
names(sdlg_size_end)[5] <- "cohort"
na.omit(sdlg_size_end)
plot(sdlg_size_end)
hist(sdlg_size_end$LA_final)
hist(sdlg_size_end$HT_final)
str(sdlg_size_end)
summary(sdlg_size_end)
cor(sdlg_size_end[1:2], method="spearman", use="complete.obs")  #Storng correlation between final leaf area and final height
# HTvLA<-lm(sdlg_size_end$LA_final ~ sdlg_size_end$HT_final)
# summary(HTvLA)

ggplot(sdlg_size_end, aes(x=HT_final, y=LA_final, color=trt)) + geom_point(shape=1)+
  geom_point(shape=1, position=position_jitter(width=1,height=.5))+ #shape 1 = hollow circles, jitter plot to seperate overlapping points
  geom_smooth(method=lm, se=FALSE)    #Add linear regression lines and Don't add shaded confidence region


hist(sdlg_size_end$LA_final)
hist(sdlg_size_end$HT_final)

GLM_HT_LA1<-glm(LA_final ~ 1, data = sdlg_size_end, family = gaussian)
GLM_HT_LA2<-glm(LA_final ~ HT_final, data = sdlg_size_end, family = gaussian)
GLM_HT_LA3<-glm(LA_final ~ trt, data = sdlg_size_end, family = gaussian)
GLM_HT_LA4<-glm(LA_final ~ HT_final+trt, data = sdlg_size_end, family = gaussian)
GLM_HT_LA5<-glm(LA_final ~ HT_final*trt, data = sdlg_size_end, family = gaussian)
GLM_HT_LA6<-glm(LA_final ~ HT_final*trt+block, data = sdlg_size_end, family = gaussian)

anova(GLM_HT_LA1,GLM_HT_LA2, test = "Chisq")  #imp fit by adding HT over just intercept
anova(GLM_HT_LA1,GLM_HT_LA3, test = "Chisq")  #no imp fit by adding trt over just intercept (but close)
anova(GLM_HT_LA2,GLM_HT_LA3, test = "Chisq")  #no imp fit by adding trt instead of ht
anova(GLM_HT_LA2,GLM_HT_LA4, test = "Chisq")  #adding trt and ht better than just ht 
anova(GLM_HT_LA4,GLM_HT_LA5, test = "Chisq")  #adding trt and trt*ht interaction best

AIC(GLM_HT_LA1,GLM_HT_LA2,GLM_HT_LA3,GLM_HT_LA4,GLM_HT_LA5,GLM_HT_LA6) #adding block doesn't result in lower AIC

# Anova version (can't use, unbalanced design)
# LAreg <- lm(LA_final ~ HT_final*trt, data =sdlg_size_end )
# anova(LAreg)
# summary(LAreg)

# 
# COHORT 2 (can combine later)
COHORT2.wide<-COHORT2.wide[with(COHORT2.wide, order(sdlg.id.no, sdlg.type, block, trt)), ]
ht.cohort2<-ht.cohort2[with(ht.cohort2, order(seedling.id.no, sdlg.type, block, trt)), ]
sdlg_size_C2end<-as.data.frame(COHORT2.wide$LAt4)
sdlg_size_C2end<-cbind(sdlg_size_C2end,as.data.frame(ht.cohort2$ht.4), as.data.frame(ht.cohort2$trt), as.data.frame(ht.cohort2$block),as.data.frame(ht.cohort2$cohort))
names(sdlg_size_C2end)[1] <- "LA_final"
names(sdlg_size_C2end)[2] <- "HT_final"
names(sdlg_size_C2end)[3] <- "trt"
names(sdlg_size_C2end)[4] <- "block"
names(sdlg_size_C2end)[5] <- "cohort"
na.omit(sdlg_size_C2end)
plot(sdlg_size_C2end)
hist(sdlg_size_C2end$LA_final)
hist(sdlg_size_C2end$HT_final)
str(sdlg_size_C2end)
summary(sdlg_size_C2end)
cor(sdlg_size_end[1:2], method="spearman", use="complete.obs")  #Storng correlation between final leaf area and final height

ggplot(sdlg_size_C2end, aes(x=HT_final, y=LA_final, color=trt)) + geom_point(shape=1)+
  geom_point(shape=1, position=position_jitter(width=1,height=.5))+ #shape 1 = hollow circles, jitter plot to seperate overlapping points
  geom_smooth(method=lm, se=FALSE)    #Add linear regression lines and Don't add shaded confidence region

hist(sdlg_size_C2end$LA_final)
hist(sdlg_size_C2end$HT_final)

GLM_HT_LA1<-glm(LA_final ~ 1, data = sdlg_size_C2end, family = gaussian)
GLM_HT_LA2<-glm(LA_final ~ HT_final, data = sdlg_size_C2end, family = gaussian)
GLM_HT_LA3<-glm(LA_final ~ trt, data = sdlg_size_C2end, family = gaussian)
GLM_HT_LA4<-glm(LA_final ~ HT_final+trt, data = sdlg_size_C2end, family = gaussian)
GLM_HT_LA5<-glm(LA_final ~ HT_final*trt, data = sdlg_size_C2end, family = gaussian)
GLM_HT_LA6<-glm(LA_final ~ HT_final*trt+block, data = sdlg_size_C2end, family = gaussian)

anova(GLM_HT_LA1,GLM_HT_LA2, test = "Chisq")  #imp fit by adding HT over just intercept
anova(GLM_HT_LA1,GLM_HT_LA3, test = "Chisq")  #no imp fit by adding trt over just intercept 
anova(GLM_HT_LA2,GLM_HT_LA3, test = "Chisq")  #no imp fit by adding trt instead of ht
anova(GLM_HT_LA2,GLM_HT_LA4, test = "Chisq")  #adding trt and ht better than just ht 
anova(GLM_HT_LA4,GLM_HT_LA5, test = "Chisq")  #adding trt and trt*ht interaction best

AIC(GLM_HT_LA1,GLM_HT_LA2,GLM_HT_LA3,GLM_HT_LA4,GLM_HT_LA5,GLM_HT_LA6) #for cohort 2 best model fit is just ht


###NEED TO COMPARE a) end of year 1 for cohorts b) cohort 1 year 1 vs. 3


# Both Together 
BOTH<-rbind(sdlg_size_end, sdlg_size_C2end)
summary(BOTH)
BOTH$cohort<-as.factor(BOTH$cohort)
BOTH$trt<-as.factor(BOTH$trt)
BOTH$block<-as.factor(BOTH$block)


GLM_HT_LA1<-glm(LA_final ~ 1, data = BOTH, family = gaussian)
GLM_HT_LA2<-glm(LA_final ~ HT_final, data = BOTH, family = gaussian)
GLM_HT_LA3<-glm(LA_final ~ trt, data = BOTH, family = gaussian)
GLM_HT_LA4<-glm(LA_final ~ HT_final+trt, data = BOTH, family = gaussian)
GLM_HT_LA5<-glm(LA_final ~ HT_final*trt, data = BOTH, family = gaussian)
GLM_HT_LA6<-glm(LA_final ~ HT_final*trt+block, data = BOTH, family = gaussian)
GLM_HT_LA7<-glm(LA_final ~ HT_final*trt*cohort, data = BOTH, family = gaussian)
GLM_HT_LA8<-glm(LA_final ~ HT_final*trt*cohort+block, data = BOTH, family = gaussian)

model.names <- c("1 Intercept", "2 Height", "3 Treatment", "4 Height + Treatment", "5 Height*Treatment", "6 Height*Treatment+Block", "7 Height*Treatment*Cohort", "8 Height*Treatment*Cohort+Block")

aov(GLM_HT_LA7)
summary(GLM_HT_LA7)

anova(GLM_HT_LA1,GLM_HT_LA2, test = "Chisq")  #imp fit by adding HT over just intercept
anova(GLM_HT_LA1,GLM_HT_LA3, test = "Chisq")  #no imp fit by adding trt over just intercept 
anova(GLM_HT_LA2,GLM_HT_LA3, test = "Chisq")  #no imp fit by adding trt instead of ht
anova(GLM_HT_LA2,GLM_HT_LA4, test = "Chisq")  #adding trt and ht better than just ht 
anova(GLM_HT_LA4,GLM_HT_LA5, test = "Chisq")  #adding trt and trt*ht interaction best
anova(GLM_HT_LA5,GLM_HT_LA7, test = "Chisq")  #adding trt*ht*cohort improvies
anova(GLM_HT_LA7,GLM_HT_LA8, test = "Chisq")  #adding block to 7 doesn't!


AIC(GLM_HT_LA1,GLM_HT_LA2,GLM_HT_LA3,GLM_HT_LA4,GLM_HT_LA5,GLM_HT_LA6,GLM_HT_LA7,GLM_HT_LA8) #for cohort 2 best model fit is just ht

# TO REPORT: http://www.ashander.info/posts/2015/10/model-selection-glms-aic-what-to-report/

# THE JAMIE WAY
summ.table <- do.call(rbind, lapply(list(GLM_HT_LA1, GLM_HT_LA2, GLM_HT_LA3,GLM_HT_LA4,GLM_HT_LA5,GLM_HT_LA6,GLM_HT_LA7,GLM_HT_LA8), broom::glance))
summ.table
table.cols <- c("df.residual", "deviance", "AIC")
reported.table <- summ.table[table.cols]
names(reported.table) <- c("Resid. Df", "Resid. Dev", "AIC")

reported.table[['dAIC']] <-  with(reported.table, AIC - min(AIC))
reported.table[['weight']] <- with(reported.table, exp(- 0.5 * dAIC) / sum(exp(- 0.5 * dAIC)))
reported.table$AIC <- NULL
reported.table$weight <- round(reported.table$weight, 2)
reported.table$dAIC <- round(reported.table$dAIC, 1)
row.names(reported.table) <- model.names
reported.table

# THE BOLKER WAY 
reported.table2 <- bbmle::AICtab(GLM_HT_LA1, GLM_HT_LA2, GLM_HT_LA3,GLM_HT_LA4,GLM_HT_LA5,GLM_HT_LA6,GLM_HT_LA7,GLM_HT_LA8, weights = TRUE, sort = FALSE, mnames = model.names)
reported.table2[["Resid. Dev"]]  <- summ.table[["deviance"]] # get the deviance from broom'd table
reported.table2





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

bmassC1<-filter(bmass,  cohort == "1") #cohort 1 after 2 years
bmassC1$cohort<-as.factor(bmassC1$cohort)

bmassC2<-filter(bmass,  cohort == "2") #cohort 2 after 1 year
bmassC2$cohort<-as.factor(bmassC2$cohort)

# VISUALIZATIONS



# HISTOGRAMS
# COHORT 1
hist(bmassC1$RSratio)
hist(bmassC1$Total.bmass)
# COHORT 2
hist(bmassC2$RSratio)
hist(bmassC2$Total.bmass)


# BOX PLOTS
# COHORT 1
boxplot(RSratio~trt,data=bmassC1) #Cohort 1 RS final
boxplot(Total.bmass~trt,data=bmassC1) #Cohort 1 Bmass final
# COHORT 2
boxplot(RSratio~trt,data=bmassC2) #Cohort 2 RS final 
boxplot(Total.bmass~trt,data=bmassC2) #Cohort bmass final

# PLOTS OF RS Ratio (Y) vs Total Biomass (X) BY TREATMENT
# COHORT 1

ggplot(bmassC1, aes(x=Total.bmass, y=RSratio, color=trt)) + geom_point(shape=1)+
  geom_point(shape=1, position=position_jitter(width=1,height=.5))+ #shape 1 = hollow circles, jitter plot to seperate overlapping points
  geom_smooth(method=lm, se=FALSE)    #Add linear regression lines and Don't add shaded confidence region

# COHORT 2
aov.bmC2<-aov(RSratio ~ trt+Total.bmass+block/trt, data = bmassC2)
summary(aov.bmC2)
ggplot(bmassC2, aes(x=Total.bmass, y=RSratio, color=trt)) + geom_point(shape=1)+
  geom_point(shape=1, position=position_jitter(width=1,height=.5))+ #shape 1 = hollow circles, jitter plot to seperate overlapping points
  geom_smooth(method=lm, se=FALSE)    #Add linear regression lines and Don't add shaded confidence region



##########
# ANALYSES (note: with focal only)
##########

##########
# COHORT 1
##########

# as ANCOVA (note there is an interaciton between trt and covariate, meaning final biomass could be influenced by density)
aov.bmC1<-aov(RSratio ~ trt * Total.bmass+Error(block:trt), data = bmassC1)
summary(aov.bmC1)



glm.1<-glm(RSratio ~ 1, family = gaussian, data = bmassC2)
summary(glm.1)
glm.2<-glm(RSratio ~ block, family = gaussian, data = bmassC2)
summary(glm.2)
glm.3<-glm(RSratio ~ trt, family = gaussian, data = bmassC2)
summary(glm.3)
glm.4<-glm(RSratio ~ Total.bmass, family = gaussian, data = bmassC2)
summary(glm.4)
glm.5<-glm(RSratio ~ trt+block, family = gaussian, data = bmassC2)
summary(glm.5)
glm.6<-glm(RSratio ~ trt+Total.bmass, family = gaussian, data = bmassC2)
summary(glm.6)
glm.7<-glm(RSratio ~ trt*Total.bmass+block, family = gaussian, data = bmassC2)
summary(glm.7)
# 
anova(glm.1,glm.2, test = "Chisq")  #no imp fit by adding block over just intercept
anova(glm.1,glm.3, test = "Chisq")  #imp fit by adding trt over just intercept
anova(glm.1,glm.4, test = "Chisq")  #no imp fit by adding total bmass over just intercept (though lose...0.06)
anova(glm.1,glm.5, test = "Chisq")  #no imp fit by adding total trt+block over just intercept
anova(glm.1,glm.6, test = "Chisq")  #no imp fit by adding trt+total bmass over just intercept (though lose...0.06)
anova(glm.1,glm.7, test = "Chisq")  #no imp fit by adding trt+total bmass + block over just intercept (though lose...0.06)

AIC(glm.1, glm.2, glm.3,glm.4, glm.5, glm.6, glm.7)

# TO REPORT: http://www.ashander.info/posts/2015/10/model-selection-glms-aic-what-to-report/

# THE JAMIE WAY
summ.table <- do.call(rbind, lapply(list(glm.1, glm.2, glm.3,glm.4, glm.5, glm.6, glm.7), broom::glance))
summ.table
table.cols <- c("df.residual", "deviance", "AIC")
reported.table <- summ.table[table.cols]
names(reported.table) <- c("Resid. Df", "Resid. Dev", "AIC")

reported.table[['dAIC']] <-  with(reported.table, AIC - min(AIC))
reported.table[['weight']] <- with(reported.table, exp(- 0.5 * dAIC) / sum(exp(- 0.5 * dAIC)))
reported.table$AIC <- NULL
reported.table$weight <- round(reported.table$weight, 2)
reported.table$dAIC <- round(reported.table$dAIC, 1)
#row.names(reported.table) <- model.names
reported.table

##########
# COHORT 2
##########


# as ANCOVA (note there is an interaciton between trt and covariate, meaning final biomass could be influenced by density)
aov.bmC2<-aov(RSratio ~ trt * Total.bmass+Error(block:trt), data = bmassC2)
summary(aov.bmC2)



glm.2.1<-glm(RSratio ~ 1, family = gaussian, data = bmassC2)
summary(glm.2.1)
glm.2.2<-glm(RSratio ~ block, family = gaussian, data = bmassC2)
summary(glm.2.2)
glm.2.3<-glm(RSratio ~ trt, family = gaussian, data = bmassC2)
summary(glm.2.3)
glm.2.4<-glm(RSratio ~ Total.bmass, family = gaussian, data = bmassC2)
summary(glm.2.4)
glm.2.5<-glm(RSratio ~ trt+block, family = gaussian, data = bmassC2)
summary(glm.2.5)
glm.2.6<-glm(RSratio ~ trt+Total.bmass, family = gaussian, data = bmassC2)
summary(glm.2.6)
glm.2.7<-glm(RSratio ~ trt*Total.bmass+block, family = gaussian, data = bmassC2)
summary(glm.2.7)
# 
anova(glm.2.1,glm.2.2, test = "Chisq")  #no imp fit by adding block over just intercept
anova(glm.2.1,glm.2.3, test = "Chisq")  #imp fit by adding trt over just intercept
anova(glm.2.1,glm.2.4, test = "Chisq")  #no imp fit by adding total bmass over just intercept (though lose...0.06)
anova(glm.2.1,glm.2.5, test = "Chisq")  #no imp fit by adding total trt+block over just intercept
anova(glm.2.1,glm.2.6, test = "Chisq")  #no imp fit by adding trt+total bmass over just intercept (though lose...0.06)
anova(glm.2.1,glm.2.7, test = "Chisq")  #no imp fit by adding trt+total bmass + block over just intercept (though lose...0.06)

AIC(glm.2.1, glm.2.2, glm.2.3,glm.2.4, glm.2.5, glm.2.6, glm.2.7)

# TO REPORT: http://www.ashander.info/posts/2015/10/model-selection-glms-aic-what-to-report/

# THE JAMIE WAY
summ.table2 <- do.call(rbind, lapply(list(glm.2.1, glm.2.2, glm.2.3,glm.2.4, glm.2.5, glm.2.6, glm.2.7), broom::glance))
summ.table2
table.cols2 <- c("df.residual", "deviance", "AIC")
reported.table2 <- summ.table2[table.cols2]
names(reported.table2) <- c("Resid. Df", "Resid. Dev", "AIC")

reported.table2[['dAIC']] <-  with(reported.table2, AIC - min(AIC))
reported.table2[['weight']] <- with(reported.table2, exp(- 0.5 * dAIC) / sum(exp(- 0.5 * dAIC)))
reported.table2$AIC <- NULL
reported.table2$weight <- round(reported.table2$weight, 2)
reported.table2$dAIC <- round(reported.table2$dAIC, 1)
#row.names2(reported.table2) <- model.names
reported.table2































































#########################################
# Analysis of LA - Experimental Seedlings
#########################################

# #choose seedlings focal + competitors
# add & sdlg.type == "focal" if you only want focal seedlings

rgr.test.C2<-filter(LARGR_BOTH,  duration == "annual" & cohort =="2") #cohort 2 after 1 year

rgr.test.C1<-filter(LARGR_BOTH,  duration == "multiyear" & cohort =="1") # cohort 1 after 2 years
rgr.test.C1.1<-filter(LARGR_BOTH,  cohort_year == "cohort_yr_1" & duration == "annual" & cohort =="1") # cohort 1 after 1 years
rgr.test.C1.2<-filter(LARGR_BOTH,  cohort_year == "cohort_yr_2" & duration == "annual" & cohort =="1") # cohort 1 YEAR 2 ONLY

rgr.testLA<-rbind(rgr.test.C2,rgr.test.C1.1) #cohort 1 and coghort 2 after 1 year

rgr.test.C1$cohort<-as.factor(rgr.test.C1$cohort)
rgr.test.C1$block<-as.factor(rgr.test.C1$block)
is.na(rgr.test.C1) <- do.call(cbind,lapply(rgr.test.C1, is.infinite)) #remove NAs
rgr.test.C1<-na.omit(rgr.test.C1)
rgr.test.C1<-droplevels(rgr.test.C1)

rgr.test.C2$cohort<-as.factor(rgr.test.C2$cohort)
rgr.test.C2$block<-as.factor(rgr.test.C2$block)
is.na(rgr.test.C2) <- do.call(cbind,lapply(rgr.test.C2, is.infinite)) #remove NAs
rgr.test.C2<-na.omit(rgr.test.C2)
rgr.test.C2<-droplevels(rgr.test.C2)

rgr.testLA$cohort<-as.factor(rgr.testLA$cohort)
rgr.testLA$block<-as.factor(rgr.testLA$block)
is.na(rgr.testLA) <- do.call(cbind,lapply(rgr.testLA, is.infinite)) #remove NAs
rgr.testLA<-na.omit(rgr.testLA)
rgr.testLA<-droplevels(rgr.testLA)
rgr.testLA<-& sdlg.type == "focal"

# http://conjugateprior.org/2013/01/formulae-in-r-anova/
# aov formula when A is random factor, B is fixed, and B is nested within A.
# A=block B=trt
# aov(Y ~ B + Error(A/B), data=d)

# foo<-lm(rgr.la ~ trt + block/trt, data=rgr.test.C1)
# anova(foo)
# summary(foo)

# using lme4
# lmer(Y ~ B + (1 | A), data=d)
# lmer(Y ~ 1 + B + (1 | A), data=d)
# lmer(rgr.la ~ trt + (1 | block), data=rgr.test.C1)
# # lmer(rgr.la ~ 1 + trt + (1 | block), data=rgr.test.C1)
# If A is random, B is fixed, and B is nested within A then
# lmer(Y ~ B + (1 | A:B), data=d)
# lmer(rgr.la ~ trt + (1 | block:trt), data=rgr.test.C1)


# Anovas
# Cohort 1 YEAR AFTER 2
C1<-lm(rgr.la ~ trt + block/trt, data=rgr.test.C1) 
anova(C1)
summary(C1)
with(rgr.test.C1, interaction.plot(trt, block, rgr.la),ylab = "RGR (leaf area)", xlab = "treatment", trace.label = "block")
boxplot(rgr.la~trt,data=rgr.test.C1) 
boxplot(rgr.la~block,data=rgr.test.C1) 
# Cohort 1 YER 1
# C1.1<-lm(rgr.la ~ trt + block/trt, data=rgr.test.C1) 
# anova(C1.1)
# summary(C1.1)
# with(rgr.test.C1.1, interaction.plot(trt, block, rgr.la),ylab = "RGR (leaf area)", xlab = "treatment", trace.label = "block")
# boxplot(rgr.la~trt,data=rgr.test.C1.1) 
# boxplot(rgr.la~block,data=rgr.test.C1.1) 
# 
# # Cohort 1 YER 2
# C1.2<-lm(rgr.la ~ trt + block/trt, data=rgr.test.C1) 
# anova(C1.2)
# summary(C1.2)
# with(rgr.test.C1.2, interaction.plot(trt, block, rgr.la),ylab = "RGR (leaf area)", xlab = "treatment", trace.label = "block")
# boxplot(rgr.la~trt,data=rgr.test.C1.2) 
# boxplot(rgr.la~block,data=rgr.test.C1.2) 



# cohort2 all
C2<-lm(rgr.la ~ trt + block/trt, data=rgr.test.C2)
anova(C2)
summary(C2)
with(rgr.test.C2, interaction.plot(trt, block, rgr.la),ylab = "RGR (leaf area)", xlab = "treatment", trace.label = "block")
# means.barplot <- qplot(x=trt, y=rgr.la, fill=trt,
#                        data=rgr.test.C2, geom="bar", stat="identity",
#                        position="dodge")
rgr.test.C2 %>% group_by(trt) %>% summarize(avg=mean(rgr.la))
rgr.test.C2 %>% group_by(trt) %>% summarize(sd=sd(rgr.la))
boxplot(rgr.la~trt,data=rgr.test.C2) 
boxplot(rgr.la~block,data=rgr.test.C2) 
# both - all at the end of the 1sy year
# Cboth<-lm(rgr.la ~ trt*cohort + block/trt/cohort, data=rgr.testLA)
# anova(Cboth)
# summary(Cboth)
# with(rgr.test.C1, interaction.plot(trt, block, rgr.la),ylab = "RGR (leaf area)", xlab = "treatment", trace.label = "block")
# boxplot(rgr.la~trt,data=rgr.testLA) 
# boxplot(rgr.la~block,data=rgr.testLA) 
# boxplot(rgr.la~cohort,data=rgr.testLA) 

# GLM LA
glm.0<-glm(rgr.la ~ 1, family = gaussian, data = rgr.testLA)
summary(glm.0)
glm.1<-glm(rgr.la ~ block, family = gaussian, data = rgr.testLA)
summary(glm.1)
glm.2<-glm(rgr.la ~ cohort, family = gaussian, data = rgr.testLA)
summary(glm.2)
glm.3<-glm(rgr.la ~ trt, family = gaussian, data = rgr.testLA)
summary(glm.3)
glm.4<-glm(rgr.la ~ trt+block, family = gaussian, data = rgr.testLA)
summary(glm.4)
glm.5<-glm(rgr.la ~ trt+cohort, family = gaussian, data = rgr.testLA)
summary(glm.5)
glm.6<-glm(rgr.la ~ cohort+block, family = gaussian, data = rgr.testLA)
summary(glm.6)
glm.7<-glm(rgr.la ~ trt*cohort+block, family = gaussian, data = rgr.testLA)
summary(glm.7)

anova(glm.0,glm.1, test = "Chisq")  #no imp fit by adding block over just intercept
anova(glm.0,glm.2, test = "Chisq")  #imp fit by adding cohort over just intercept
anova(glm.0,glm.3, test = "Chisq")  #no imp fit by adding trt over just intercept
anova(glm.3,glm.4, test = "Chisq")  #adding block to trt sig improves fit

AIC(glm.0,glm.1, glm.2, glm.3,glm.4, glm.5, glm.6, glm.7)


#SIMPLER AS ANOVA - LA (all seedlings, nested
aov.la<-aov(rgr.la ~ trt*cohort+trt/sdlg.id.no, data = rgr.testLA)
summary(aov.la)

boxplot(rgr.la~trt*cohort,data=rgr.testLA) 


#########################################
# Analysis of HT - Experimental Seedlings
#########################################

#Reduce dataset: only include "focal" seedlings.
# str(LARGR_BOTH)

###THIS ANALYSIS IS OF RGR from start of experiment to end - no seasonal effects
rgr.test.C2<-filter(LARGR_BOTH, sdlg.type == "focal" & duration == "annual" & cohort =="2")
rgr.test.C1<-filter(LARGR_BOTH, sdlg.type == "focal" & duration == "multiyear" & cohort =="1")
rgr.testLA<-rbind(rgr.test.C2,rgr.test.C1)

rgrHTC2.test<-filter(HTRGR_BOTH, sdlg.type == "focal" & duration == "annual" & cohort =="2")
rgrHTC1.test<-filter(HTRGR_BOTH, sdlg.type == "focal" & duration == "multiyear" & cohort =="1")
rgrHT.test<-rbind(rgrHTC2,rgrHTC1)

rgrHT.test$cohort<-as.factor(rgrHT.test$cohort)
rgrHT.test$block<-as.factor(rgrHT.test$block)
is.na(rgrHT.test) <- do.call(cbind,lapply(rgrHT.test, is.infinite)) #remove NAs
rgrHT.test<-na.omit(rgrHT.test)
rgrHT.test<-droplevels(rgrHT.test)

#SIMPLER AS ANOVA - LA 
aov.ht<-aov(rgr.ht ~ trt*cohort+block, data = rgrHT.test)
summary(aov.ht)
with(rgrHT.test, interaction.plot(trt, cohort, rgr.ht),ylab = "RGR (height)", xlab = "treatment", trace.label = "Cohort")




















#with all


# logit trasnform, but actually don't do this...
# bmass<-bmass[, "logitRS"] <-log10(bmass$RSratio+/(1-bmass$RSratio))









#Figures RGR

ggplot(rgr.testLA, aes(x=trt, y=rgr.la)) + geom_boxplot() 

ggplot(data=rgr.testLA, aes(x=trt, y=rgr.la, group=sdlg.id.no)) +
  geom_line()+
  geom_point()

p<-ggplot(rgr.testLA, aes(x=trt, y=rgr.la, group=cohort)) +
  geom_line(aes(color=cohort))+
  geom_point(aes(color=cohort))
p


plot(rgr.testLA$canopy.openess.percent,rgr.test$rgr.ht)


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



