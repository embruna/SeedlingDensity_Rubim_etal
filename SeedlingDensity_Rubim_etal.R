
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
summary(Exp_Data_C2)
Exp_Data_C1<-read.csv("EXP_DATA_COHORT1_18august2014.csv", dec=".", header = TRUE, sep = ",", check.names=FALSE )
summary(Exp_Data_C1)

#####NOTE THAT STEP 2 and STEP 3 are essentially the same thing, should be made into a function!!!!

##################
###step 2: reshape data - wide to long - and lineup adjacent columns with 
###leaf lengths and percent of each leaf missing
##################

#select the columns for leaf length 
mdata.la<-Exp_Data_C2[1:33]
#Now melt them with reshape2 package to create a dataframe with 
#leaf-lengths in long form 
mdata.la <- melt(mdata.la, id.var=c("seedling.id.no","block","trt", "sdlg.no", "sdlg.type"))
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
mdata.miss<-Exp_Data_C2[,c(1,34:61)]
mdata.miss <- melt(mdata.miss, id.var=c("seedling.id.no"))
#then rename the column "value" for as "leaf.percentage.missing"
names(mdata.miss)[names(mdata.miss)=="value"] <- "leaf.percentage.missing"
#add that column to the dataframe of leaf lengths and rename that column
cohort2.la<-cbind(mdata.la,mdata.miss$leaf.percentage.missing)
names(cohort2.la)[names(cohort2.la)=="mdata.miss$leaf.percentage.missing"] <- "leaf.percentage.missing"
#the final steps are cleanup: remove NA, which leaves for each plant with records for only the leaves it actually had
cohort2.la<-na.omit(cohort2.la)
#and to add a column identifying which cohort this is (which you will need when you rbind cohort 1 and 2 for analyses)
cohort2.la[, "cohort"] <- "Cohort 2"


##################
###step 3: calclulate the area of each leaf, subtract the area missing from each leaf
###and sum to find the total leaf area of each plant in each time interval
##################
#Adds a column with the leaf area (uncorrected). LA is calculated using the formula in Bruna 2002 Oecologia 
cohort2.la[, "uncorrected.leaf.area"] <- 0.53+(0.831*cohort2.la$leaf.length)
#adds a column with the "corrected" leaf area (i.e., corrected leaf area - % missing) 
cohort2.la[, "corrected.leaf.area"] <-cohort2.la$uncorrected.leaf.area - (cohort2.la$uncorrected.leaf.area*(cohort2.la$leaf.percentage.missing/100))


##################
###step 4: DO EXACTLY THE SAME FOR COHORT 1
##################
#load the  CSV files and save them as dataframes
Exp_Data_C1<-read.csv("EXP_DATA_COHORT1_18august2014.csv", dec=".", header = TRUE, sep = ",", check.names=FALSE )
summary(Exp_Data_C1)

##################
###step 2: reshape data - wide to long - and lineup adjacent columns with 
###leaf lengths and percent of each leaf missing
##################

#select the columns for leaf length 
mdata.la.c1<-Exp_Data_C1[1:33]
#Now melt them with reshape2 package to create a dataframe with 
#leaf-lengths in long form 
mdata.la.c1 <- melt(mdata.la.c1, id.var=c("seedling.id.no","block","trt", "sdlg.no", "sdlg.type"))
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
mdata.miss.c1<-Exp_Data_C1[,c(1,34:61)]
mdata.miss.c1 <- melt(mdata.miss.c1, id.var=c("seedling.id.no"))
#then rename the column "value" for as "leaf.percentage.missing"
names(mdata.miss.c1)[names(mdata.miss.c1)=="value"] <- "leaf.percentage.missing"
#add that column to the dataframe of leaf lengths and rename that column
cohort1.la<-cbind(mdata.la.c1,mdata.miss.c1$leaf.percentage.missing)
names(cohort1.la)[names(cohort1.la)=="mdata.miss$leaf.percentage.missing"] <- "leaf.percentage.missing"
#the final steps are cleanup: remove NA, which leaves for each plant with records for only the leaves it actually had
cohort1.la<-na.omit(cohort1.la)
#and to add a column identifying which cohort this is (which you will need when you rbind cohort 1 and 2 for analyses)
cohort1.la[, "cohort"] <- "Cohort 1"


##################
###step 3: calclulate the area of each leaf, subtract the area missing from each leaf
###and sum to find the total leaf area of each plant in each time interval
##################
#Adds a column with the leaf area (uncorrected). LA is calculated using the formula in Bruna 2002 Oecologia 
cohort1.la[, "uncorrected.leaf.area"] <- 0.53+(0.831*cohort1.la$leaf.length)
#adds a column with the "corrected" leaf area (i.e., corrected leaf area - % missing) 
cohort1.la[, "corrected.leaf.area"] <-cohort1.la$uncorrected.leaf.area - (cohort1.la$uncorrected.leaf.area*(cohort1.la$leaf.percentage.missing/100))


##################
###step 4: ROW BIND THE TWO, then do some summarizing of leaf area for each plant to calclulate the RGR
##################






#################
###CALCLULATIONS OF RGR BAsed on Height & Leaf Area
##################
#Formula for RGR is rgr=[ln(W2)-ln(W1)]/t2-t1, where W is weegith at time t 
rgr.ht=(log(Exp_Data$ht.final)-log(Exp_Data$ht.initial))/Exp_Data$days
rgr.la=(log(Exp_Data$leafarea.final)-log(Exp_Data$leafarea.initial))/Exp_Data$days

hist(Exp_Data$ht.final)
sort(Exp_Data$leafarea.initial)



