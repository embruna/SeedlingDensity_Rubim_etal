
#R CODE FOR IMPORTING, MANIPULATING, AND ANALYZING THE DATASETS USED IN: Rubim et al. Seedling Density MS

##################
###Set WD and load packages you need.
##################

setwd("/Users/emiliobruna/Dropbox/Rubim Manuscripts/Seedling density[3]/EB Revision/Data-Rubim-CSV")
library(reshape2)



#################
###STEP 1: LOAD THE RAW DATA
##################
#CLear out everything from the environment 
rm(list=ls())
#load the  CSV files and save them as dataframes
Exp_Data_C2<-read.csv("EXP_DATA_COHORT2_18august2014.csv", dec=".", header = TRUE, sep = ",", check.names=FALSE )
#this makes sure treatment is an ordered factor (i.e., 1<2<4)
Exp_Data_C2$trt <- ordered(Exp_Data_C2$trt, levels = c("one", "two", "four"))
summary(Exp_Data_C2)


#####NOTE THAT STEP 2 and STEP 3 are essentially the same thing, should be made into a function!!!!

##################
###step 2: reshape data - wide to long - and lineup adjacent columns with 
###leaf lengths and percent of each leaf missing
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
cohort2.la[, "uncorrected.leaf.area"] <- 0.53+(0.831*cohort2.la$leaf.length)

#first removes NA from the column of percet of each leaf missing (in some cases NA because there was no measurment was taken
#then adds a column with the "corrected" leaf area (i.e., corrected leaf area - % missing) 
cohort2.la$leaf.percentage.missing[is.na(cohort2.la$leaf.percentage.missing)] <- 0
cohort2.la[, "corrected.leaf.area"] <-cohort2.la$uncorrected.leaf.area - (cohort2.la$uncorrected.leaf.area*(cohort2.la$leaf.percentage.missing/100))


##################
###DO EXACTLY THE SAME FOR COHORT 1
##################
#load the  CSV files and save them as dataframes
Exp_Data_C1<-read.csv("EXP_DATA_COHORT1_23august2014.csv", dec=".", header = TRUE, sep = ",", check.names=FALSE )
#this makes sure treatment is an ordered factor (i.e., 1<2<4)
Exp_Data_C1$trt <- ordered(Exp_Data_C1$trt, levels = c("one", "two", "four"))
summary(Exp_Data_C1)

##################
###step 2: reshape data - wide to long - and lineup adjacent columns with 
###leaf lengths and percent of each leaf missing
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
#Adds a column with the leaf area (uncorrected). LA is calculated using the formula in Bruna 2002 Oecologia 
cohort1.la[, "uncorrected.leaf.area"] <- 0.53+(0.831*cohort1.la$leaf.length)
#first removes NA from the column of percet of each leaf missing (in some cases NA because there was no measurment was taken
#then adds a column with the "corrected" leaf area (i.e., corrected leaf area - % missing) 
cohort1.la$leaf.percentage.missing[is.na(cohort1.la$leaf.percentage.missing)] <- 0
cohort1.la[, "corrected.leaf.area"] <-cohort1.la$uncorrected.leaf.area - (cohort1.la$uncorrected.leaf.area*(cohort1.la$leaf.percentage.missing/100))

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
COHORT1.long<-COHORT1.long[with(COHORT1.long, order(sdlg.type, block, trt, sdlg.id.no)), ]
COHORT2.long<-COHORT2.long[with(COHORT2.long, order(sdlg.type, block, trt, sdlg.id.no)), ]

#################
###CALCLULATIONS OF RGR BAsed on Leaf Area
##################
#Formula for RGR is rgr=[ln(W2)-ln(W1)]/t2-t1, where W is weegith at time t 

###Need to
COHORT1.long[, "rgr1.9"] <-(log(COHORT1.long$t9)-log(COHORT1.long$t1))/773

hist(Exp_Data$ht.final)
sort(Exp_Data$leafarea.initial)


#################
###CALCLULATIONS OF RGR BAsed on seedling height
##################
#rgr.ht=(log(Exp_Data$ht.final)-log(Exp_Data$ht.initial))/Exp_Data$days
