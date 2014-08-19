
#R CODE FOR IMPORTING, MANIPULATING, AND ANALYZING THE DATASETS USED IN: Rubim et al. Seedling Density MS

##################
###Set WD and load packages you need.
##################

setwd("/Users/emiliobruna/Dropbox/Rubim Manuscripts/Seedling density[3]/EB Revision")
library(gdata)
library(ggplot2)
library(reshape)



#################
###DATA ENTRY AND CLEANUP
##################

#CLear out everything from the environment 
rm(list=ls())

#load the  CSV files and save them as dataframes
Exp_Data<-read.csv("cohort2-leaf_length.csv", dec=".", header = TRUE, sep = ",", check.names=FALSE )
Exp_Data_ht<-read.csv("cohort2-ht.csv", dec=".", header = TRUE, sep = ",", check.names=FALSE )

#look at the summary - all look ok?
summary(Exp_Data)


#step 2: reshape data - wide to long

mdata <- melt(Exp_Data, id=c("block","trt", "sdlg.no", "sdlg.type"))
mdata2 <- mdata[order(mdata$block,mdata$trt, mdata$sdlg.no, mdata$sdlg.type),]
mdata3 <- na.omit(mdata2)

mdata.ht <- melt(Exp_Data_ht, id=c("block","trt", "sdlg.no", "sdlg.type"))
mdata2.ht <- mdata.ht[order(mdata.ht$block,mdata.ht$trt, mdata.ht$sdlg.no, mdata.ht$sdlg.type),]
mdata3.ht <- na.omit(mdata2.ht)


write.csv(mdata3, "cohort2_long")

#step 3: calclulate leaf area for each leaf
#calclualte leaf area missing from each leaf
#subtract missing leaf area to get final leaf area
#sum to get total plant leaf area in each itme period
#calclulate realtive growth rate t0 vs tfinal for each cohort
#
#

#################
###CALCLULATIONS OF RGR BAsed on Height & Leaf Area
##################
#Formula for RGR is rgr=[ln(W2)-ln(W1)]/t2-t1, where W is weegith at time t 
rgr.ht=(log(Exp_Data$ht.final)-log(Exp_Data$ht.initial))/Exp_Data$days
rgr.la=(log(Exp_Data$leafarea.final)-log(Exp_Data$leafarea.initial))/Exp_Data$days

hist(Exp_Data$ht.final)
sort(Exp_Data$leafarea.initial)



