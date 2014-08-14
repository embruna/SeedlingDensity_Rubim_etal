
#R CODE FOR IMPORTING, MANIPULATING, AND ANALYZING THE DATASETS USED IN: Rubim et al. Seedling Density MS
##################
###Set WD and load packages you need.
##################

setwd("/Users/emiliobruna/Dropbox/Rubim Manuscripts/Seedling density[3]/EB Revision")
library(gdata)
library(ggplot2)



#################
###DATA ENTRY AND CLEANUP
##################

#CLear out everything from the environment 
rm(list=ls())

#load the  CSV files and save them as dataframes
Exp_Data<-read.csv("Seedling_Experimental_Data_11august2014.csv", dec=".", header = TRUE, sep = ",", check.names=FALSE )
#look at the summary - all look ok?
summary(Exp_Data)

###Will need to revise this to import the raw data and calclulate RGR from Raw data
###Forst do calclualtions, then join into one file, then do ANOVA
###this is a test to see if SSH worked
####Deleted this line
#################
###CALCLULATIONS OF RGR BAsed on Height & Leaf Area
##################
#Formula for RGR is rgr=[ln(W2)-ln(W1)]/t2-t1, where W is weegith at time t 
rgr.ht=(log(Exp_Data$ht.final)-log(Exp_Data$ht.initial))/Exp_Data$days
rgr.la=(log(Exp_Data$leafarea.final)-log(Exp_Data$leafarea.initial))/Exp_Data$days

hist(Exp_Data$ht.final)
sort(Exp_Data$leafarea.initial)



