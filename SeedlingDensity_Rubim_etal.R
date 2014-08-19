
#R CODE FOR IMPORTING, MANIPULATING, AND ANALYZING THE DATASETS USED IN: Rubim et al. Seedling Density MS

##################
###Set WD and load packages you need.
##################

setwd("/Users/emiliobruna/Dropbox/Rubim Manuscripts/Seedling density[3]/EB Revision/Data-Rubim-CSV")
library(gdata)
library(ggplot2)
library(reshape2)



#################
###STEP 1: LOAD THE RAW DATA
##################
#CLear out everything from the environment 
rm(list=ls())
#load the  CSV files and save them as dataframes
Exp_Data_C2<-read.csv("EXP_DATA_COHORT2_18august2014.csv", dec=".", header = TRUE, sep = ",", check.names=FALSE )
summary(Exp_Data_C2)

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
#now add a column with % of leaf area missing to that dataframe
#to do so create a dataframe with leaf-areas missing in long form
mdata.miss<-Exp_Data_C2[,c(1,34:61)]
mdata.miss <- melt(mdata.miss, id.var=c("seedling.id.no"))
#then rename the column "value" for as "leaf.percentage.missing"
names(mdata.miss)[names(mdata.miss)=="value"] <- "leaf.percentage.missing"
#add that column to the dataframe of leaf lengths 
cohort2.la<-cbind(mdata.la,mdata.miss$leaf.percentage.missing)
#the final step is to remove NA, which leaves each plant only with the leaves it actually had
cohort2.la<-na.omit(cohort2.la)

##################
###step 3: calclulate the area of each leaf, subtract the area missing from each leaf
###and sum to find the total leaf area of each plant in each time interval
##################



  




#Height data - wide to long
mdata.ht <- melt(Exp_Data_ht, id=c("block","trt", "sdlg.no", "sdlg.type"))
mdata2.ht <- mdata.ht[order(mdata.ht$block,mdata.ht$trt, mdata.ht$sdlg.no, mdata.ht$sdlg.type),]
mdata3.ht <- na.omit(mdata2.ht)


#################
###CALCLULATIONS OF RGR BAsed on Height & Leaf Area
##################
#Formula for RGR is rgr=[ln(W2)-ln(W1)]/t2-t1, where W is weegith at time t 
rgr.ht=(log(Exp_Data$ht.final)-log(Exp_Data$ht.initial))/Exp_Data$days
rgr.la=(log(Exp_Data$leafarea.final)-log(Exp_Data$leafarea.initial))/Exp_Data$days

hist(Exp_Data$ht.final)
sort(Exp_Data$leafarea.initial)



