# Introduction ------------------------------------------------------------

# Code for importing, cleaning, & analyzing data for Rubim et al. sdlg MS


# To Do List --------------------------------------------------------------


# These are the last data entry and manipulation" things that need to be dealt
# with prior to analyses

# Email of August 28

# 1 Me explica de novo o que querem dizer os favores nas colunas morre1, 
# morre2, etc.  Acredito que era que se uma planta competidora corria era 
# trocada por outra do viveiro, eh isso? E se era uma planta focal?  
# RESPOSTA: Foi as vezes que houve mortalidade no balde. quando era
# focal outra foi sorteada  e substitituida.

# In other words, 
# 1.a.4 5.a.9 means the first plant was alive peiods intervals 1-4, then 
# replaced with another from 6-9
# DECSION: key is focal plant - so should we do analyes only with 1-4? My gut 
# is YES because it meant that one died.Or can do for both, and average the RGR?
# NOTE: make sure days corresponds to dates of intervals; currently don't.
# DECSION 2: change ht and la from dead to zero?

# Momentos  Experimento_1	dias

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




# Load Packages -----------------------------------------------------------


# library(reshape2)
library(tidyverse)
# library(lme4)
library(nlme)
library(RColorBrewer)
library(lmerTest)
library(car) # to ensure type III SS
library(visreg) # visualizing regression models 
library(emmeans) # for estimating marginal means

# Import Raw Data - Canopy Cover & Biomass --------------------------------

canopy<-read_csv("data/canopy_cover.csv")
bmass<-read_csv("data/final_sdlg_biomass.csv")

# Make treatment an ordered factor (i.e., 1<2<4)
# Exclude blocks 17-19 from Cohort 1. These were added late in the experiment 
# and are not a valid comparison
bmass$trt <- ordered(bmass$trt, levels = c("one", "two", "four"))
bmass<-bmass[!(bmass$cohort == "1" & bmass$block>16), ]

bmass %>% group_by(cohort,trt) %>% summarize(n_distinct(block))

bmass<-bmass %>% 
  mutate(Total.bmass=root.biomass+lf.biomass+stem.biomass) %>% 
  mutate(RSratio=root.biomass/(lf.biomass+stem.biomass)) %>% 
  full_join(canopy, by = "block") %>% 
  select(-location,
         -site.openess.percent, 
         -mask.area.percent, 
         -sky.area.percent, 
         -LAI.4.ring, 
         -LAI.5.ring)
bmass$cohort<-as.factor(bmass$cohort)
bmass$sdlg.type<-"focal"
bmass$block<-as.integer(bmass$block)
bmass<-arrange(bmass,cohort,block,trt)
bmass<-as_tibble(bmass)



# Import Data - Leaf Area -------------------------------------------------

# COHORT 1
Exp_Data_C1<-read_csv("data/EXP_DATA_COHORT1_26nov2014.csv")
# Make treatment an ordered factor (i.e., 1<2<4)
Exp_Data_C1$trt <- ordered(Exp_Data_C1$trt, levels = c("one", "two", "four"))
# Exclude blocks 17-19 from Cohort 1. These were added late in the 
# experiment and are not a valid comparison
Exp_Data_C1<-Exp_Data_C1[Exp_Data_C1$block<17,]
summary(Exp_Data_C1)

# COHORT 2
Exp_Data_C2<-read_csv("data/EXP_DATA_COHORT2_20nov2014.csv")
# Make treatment an ordered factor (i.e., 1<2<4)
Exp_Data_C2$trt <- ordered(Exp_Data_C2$trt, levels = c("one", "two", "four"))
# summary(Exp_Data_C2)

colnames(Exp_Data_C1)
colnames(Exp_Data_C2)




# Calculate Total Leaf Area of Each Seedling ------------------------------

# NOTE THAT STEP 2 and STEP 3 are essentially the same thing, should 
# be made into a function!!!!

# COHORT 1: reshape data (wide to long) & lineup adjacent columns with 
# leaf lengths & % of each leaf missing
# select the columns for leaf length 
# convert to long form
mdata.lengths<-Exp_Data_C1 %>%
  select(cohort:lf14.t9) %>% 
  bind_cols(select(Exp_Data_C1, days)) %>% 
  gather(leaf.time,leaf.length,lf1.t1:lf14.t9) %>% 
  arrange(cohort,block,trt,sdlg.no)

# select the columns for missing leaf area  
# convert to long form
mdata.missing<-Exp_Data_C1 %>%
  select(cohort:sdlg.type) %>% 
  bind_cols(select(Exp_Data_C1, missing.lf3.t1:days)) %>% 
  gather(leaf.time,leaf.percentage.missing,missing.lf3.t1:missing.lf14.t9) %>% 
  arrange(cohort,block,trt,sdlg.no)

mdata.missing$leaf.time<-gsub("missing.","",mdata.missing$leaf.time)

# join the leaf leangths and misisng area
mdata.la.c1<- left_join(mdata.lengths,mdata.missing) %>% 
  separate(leaf.time, c("leaf", "interval")) %>%   # separate leaf no
  # and time interval into two columns
  drop_na(leaf.length) 

# Calculate the area of each leaf, subtract the area missing from each leaf,
# and sum to find the total leaf area of each plant in each time interval
cohort1.la<-mdata.la.c1 %>% 
  mutate(uncorrected.leaf.area=(0.53+(0.831*leaf.length))) %>%  
  # Adds a column with the leaf area (uncorrected). LA is calculated using the 
  # formula in Bruna 2002 Oecologia 
  replace_na(list(leaf.percentage.missing = 0)) %>% 
  mutate(corrected.leaf.area=uncorrected.leaf.area - 
           (uncorrected.leaf.area*leaf.percentage.missing/100))


# COHORT 2: reshape data (wide to long) & lineup adjacent columns
# with leaf lengths & % of each leaf missing
str(Exp_Data_C2)
Exp_Data_C2$lf6.t2<-as.numeric(Exp_Data_C2$lf6.t2)
Exp_Data_C2$lf7.t2<-as.numeric(Exp_Data_C2$lf7.t2)
Exp_Data_C2$lf8.t2<-as.numeric(Exp_Data_C2$lf8.t2)
Exp_Data_C2$lf8.t4<-as.numeric(Exp_Data_C2$lf8.t4)

# Exp_Data_C2$lf6.t2<-NULL
# Exp_Data_C2$lf7.t2<-NULL
# Exp_Data_C2$lf8.t2<-NULL
# Exp_Data_C2$lf8.t4<-NULL

# select the columns for leaf length 
# convert to long form
Exp_Data_C2$lf7.t4[Exp_Data_C2$lf7.t4=="x"]<-NA
Exp_Data_C2$lf7.t4<-as.numeric(as.character(Exp_Data_C2$lf7.t4))
mdata.lengths2<-Exp_Data_C2 %>%
  select(cohort:lf8.t4) %>% 
  bind_cols(select(Exp_Data_C2, days)) %>% 
  gather(leaf.time,leaf.length,lf1.t1:lf8.t4) %>% 
  arrange(cohort,block,trt,sdlg.no) 

# select the columns for missing leaf area  
# convert to long form
mdata.missing2<-Exp_Data_C2 %>%
  select(cohort:sdlg.type) %>% 
  bind_cols(select(Exp_Data_C2, missing.lf1.t1:days)) %>% 
  gather(leaf.time,leaf.percentage.missing,missing.lf1.t1:missing.lf8.t4) %>% 
  arrange(cohort,block,trt,sdlg.no)

mdata.missing2$leaf.time<-gsub("missing.","",mdata.missing2$leaf.time)

# join the leaf leangths and misisng area
mdata.la.c2<- left_join(mdata.lengths2,mdata.missing2) %>% 
  separate(leaf.time, c("leaf", "interval")) %>%   
  # separate leaf no  and time interval into two columns
  drop_na(leaf.length) 
mdata.la.c2$leaf.length<-as.numeric(mdata.la.c2$leaf.length)

# Calculate the area of each leaf, subtract the area missing from each leaf,
# and sum to find the total leaf area of each plant in each time interval
cohort2.la<-mdata.la.c2 %>% 
  mutate(uncorrected.leaf.area=(0.53+(0.831*leaf.length))) %>%  # Adds a column 
  # with the leaf area (uncorrected). LA is calculated using 
  # the formula in Bruna 2002 Oecologia 
  replace_na(list(leaf.percentage.missing = 0)) %>% 
  mutate(corrected.leaf.area=uncorrected.leaf.area -
           (uncorrected.leaf.area*leaf.percentage.missing/100))


# STEP 4: AGGREGATE DATA FOR ANALYSIS  


# Formula for RGR is rgr=[ln(W2)-ln(W1)]/t2-t1, where W is wight at time t 

# TODO: # NEED TO FIGURE OUT WHIHC OF THESE IS TRUE - 
# ARE PLANTS WITH ZERO IN LAST COLUMN ALIVE OR DEAD???
# Some of the plants had a leaf area of "0" in the last measuremnt
# If this means the plant is dead, then use this to convert those values to NA 
# so that they don't mess up the calclulation of RGR
# COHORT1.wide$LAt9[COHORT1.wide$LAt9== 0] <- NA
# COHORT2.wide$LAt4[COHORT2.wide$LAt4== 0] <- NA


COHORT1<-cohort1.la %>% 
  group_by(cohort,block,trt,sdlg.type,sdlg.no,interval,days,seedling.id.no) %>% 
  summarize(leaf.area=sum(corrected.leaf.area)) %>% 
  select(cohort, 
         block, 
         trt, 
         sdlg.type, 
         sdlg.no, 
         seedling.id.no, 
         interval, 
         days, 
         leaf.area)

COHORT1.wide<-COHORT1 %>% spread(interval,leaf.area) %>% 
  mutate(rgrLA_wet1=(log(t3)-log(t1))/84) %>%  #t1-t3
  mutate(rgrLA_dry1=(log(t5)-log(t3))/(282-84)) %>% #t3-t5
  mutate(rgrLA_wet2=(log(t7)-log(t5))/(495-282)) %>% #t5-t7
  mutate(rgrLA_dry2=(log(t8)-log(t7))/(684-495)) %>%  #t7-t8
  mutate(rgrLA_wet3=(log(t9)-log(t8))/(days-684)) %>%  #t7-t8
  mutate(rgrLA_yr1=(log(t6)-log(t1))/369) %>%     #5 08 march 08 - 29  march 09 
  mutate(rgrLA_yr2=(log(t9)-log(t6))/(days-396)) %>% #March 8 08 
  mutate(rgrLA_yrs1.2=(log(t9)-log(t1))/days) 


COHORT2<-cohort2.la %>% 
  group_by(cohort,block,trt,sdlg.type,sdlg.no,interval,days,seedling.id.no) %>% 
  summarize(leaf.area=sum(corrected.leaf.area)) %>% 
  select(cohort, 
         block, 
         trt, 
         sdlg.type, 
         sdlg.no, 
         seedling.id.no, 
         interval, 
         days, 
         leaf.area)

COHORT2.wide<-COHORT2 %>% spread(interval,leaf.area) %>% 
  mutate(rgrLA_wet2=(log(t2)-log(t1))/91) %>% 
  mutate(rgrLA_dry2=(log(t3)-log(t1))/(280-91)) %>% 
  mutate(rgrLA_wet3=(log(t4)-log(t3))/(days-280)) %>% 
  mutate(rgrLA_yr1=(log(t4)-log(t1))/days) 

all.plants<-bind_rows(COHORT1.wide,COHORT2.wide) 
all.plants$trt <- ordered(all.plants$trt, levels = c("one", "two", "four"))  
all.plants$cohort<-as.factor(all.plants$cohort)

all.plants<-left_join(all.plants,bmass)
all.plants$index <- seq.int(nrow(all.plants))

# Add an indicator of whether the seedlings planted in the experiment 
# survived to the end of the study

surv1<-all.plants %>% filter(cohort=="1") %>% ungroup() %>% select(index,t1:t9)
surv1<-surv1 %>% select(index) %>% 
  mutate(data.frame(ad=rowSums(is.na(surv1) > 0)))

surv2<-all.plants %>% filter(cohort=="2") %>% ungroup() %>% select(index,t1:t4)
surv2<-surv2 %>% select(index) %>% 
  mutate(data.frame(ad=rowSums(is.na(surv2) > 0)))

surv3<-bind_rows(surv1,surv2) %>%
  mutate(survive=ifelse(ad>0,"dead","alive")) 

all.plants<-full_join(all.plants,surv3)
rm(surv3,surv2,surv1)
all.plants$survive<-as.factor(all.plants$survive)

# Exclude blocks 17-19 from Cohort 1. These were added late in the experiment 
# and are not a valid comparison
all.plants <- all.plants %>% 
  filter(block<17)


all.plants$block <- as.factor(all.plants$block)

write.csv(all.plants, "./data_clean/all_plants.csv")





hist(all.plants$canopy.openess.percent)
summary(all.plants$canopy.openess.percent)



# Analyses - Initial Leaf Area ---------------------------------------------

# Testing for differences in Initial LA

# https://www.zoology.ubc.ca/~schluter/R/Model.html#Two_fixed_and_a_random%E2%80%93split_plot
# INIT_m1 <- glmer(t1~trt+cohort+(1|block), data = all.plants, family=Gamma)
INIT_m1 <- lmer(t1~trt*cohort+(1|block), data = all.plants)
# Plot the residuals against the fitted values in a partial check of assumptions.
plot(INIT_m1)
# To visualize the model fit, visreg() seems to work in this case, 
# since there is only one random factor.
visreg(INIT_m1, xvar = "trt", by = "cohort", scales=list(rot = 90))
# Use emmeans to obtain model-based estimates of the treatment (A) means.
# means of each treatment combination
emmeans(INIT_m1, c("trt", "cohort"), data = all.plants) 

# means of "A" treatment levels averaged over levels of "B"
emmeans(INIT_m1, c("trt"), data = all.plants)      

summary(INIT_m1)
VarCorr(INIT_m1)
# Use anova() to test the fixed effects (the grand mean and the treatment A
anova(INIT_m1, type = 3)

# Summary Table
init_summary<-all.plants %>% 
  group_by(cohort,trt) %>% 
  summarise(avg_LA_init=mean(t1),sd_LA_init=sd(t1))
init_summary


# Figure - Initial Leaf Area ----------------------------------------------

plot_init_la<-ggplot(init_summary, aes(x=trt, y=avg_LA_init))+
  geom_bar(stat="summary") +
  ylab("Initial leaf area (mean cm2 +/- 1 SD)") +
  xlab("Seedling density") +
  ggtitle("Initial Leaf Area") + 
  facet_grid(cols=vars(cohort), scales="fixed")+
  scale_y_continuous(expand = c(0, 0))+
  geom_errorbar(aes(ymin=avg_LA_init-sd_LA_init, ymax=avg_LA_init+sd_LA_init),
                width=.2)

# scale_colour_manual(values=c("gray54", "orangered2","darkblue"))

plot_init_la<-plot_init_la+ theme_classic()+theme(legend.direction = 'horizontal', 
                                                  legend.position = c(0.1,0.85),
                                                  plot.title = element_text(color="black", size=18, hjust=0.05, vjust=-.2),
                                                  legend.background = element_rect(size=0.5, linetype="solid", color="black"),
                                                  axis.title.x=element_text(colour="black", size = 18, vjust=-0.5),            #sets x axis title size, style, distance from axis #add , face = "bold" if you want bold
                                                  axis.title.y=element_text(colour="black", size = 18, vjust=2),            #sets y axis title size, style, distance from axis #add , face = "bold" if you want bold
                                                  axis.text=element_text(colour="black", size = 14),                             #sets size and style of labels on axes
                                                  legend.title = element_blank(), #remove title of legend
                                                  legend.text = element_text(color="black", size=14))   #plot margin - top, right, bottom, left
plot_init_la






# Analyses - Survivorship -------------------------------------------------
 
# Cohort 1: 4 focal from trt=4 died and need to be replaced for the analyses (block 2,6,8,9)
# Cohort 2: 1 focal from trt=2 died and need to be replaced for the analyses (block 19)
# Cohort 2: 3 focal from trt=4 died and need to be replaced for the analyses (block 9,15,18)

# Orgabize the data to analyze survival
# Will get Number initial, number final, if a pot had deaths (true/false), and
# How manys eeddlings died in a pot

fisher_data<-all.plants %>% 
  group_by(cohort,block,trt,survive) %>% 
  summarize(Nfinal=n()) %>% 
  filter(survive=="alive") %>% 
  mutate(Ninit = case_when(trt == "one" ~ 1, 
                           trt == "two" ~ 2,
                           trt == "four" ~ 4)) %>% 
  mutate(death_in_pots=Ninit-Nfinal) %>% 
  mutate(pot_had_death=death_in_pots>0) %>% 
  ungroup()
fisher_data
###########

# Overall survival by cohort

surv_overall<-fisher_data %>% 
  group_by(cohort) %>% 
  summarize(alive = sum(Nfinal), 
            dead = sum(death_in_pots),
            blocks = (n_distinct(block))) %>% 
  mutate(perc_mort = dead / (alive + dead)*100) %>% 
  mutate(total_sdlgs = alive + dead) %>% 
  select(cohort, 
         blocks, 
         total_sdlgs, 
         alive, 
         dead, 
         perc_mort)

surv_overall
###########

# Mortality by treatment and cohort
surv_trt_cohort<-fisher_data %>% 
  group_by(cohort,trt) %>% 
  summarize(alive=sum(Nfinal), dead=sum(death_in_pots)) %>% 
  mutate(total_sdlgs = alive + dead) %>% 
  mutate(perc_mort=dead/(alive+dead)*100)
surv_trt_cohort


# figure - mortality % -------------------------------------------------


plot_mort_perc<-ggplot(surv_trt_cohort, aes(x=trt, y=perc_mort))+
  geom_bar(stat="summary") +
  ylab("Seedlings dying (%)") +
  xlab("Seedling density") +
  ggtitle("Percent of Seedlings in each Treatmnent Dying (pooled across pots)") + 
  facet_grid(cols=vars(cohort), scales="fixed")+
  scale_y_continuous(expand = c(0, 0))

# scale_colour_manual(values=c("gray54", "orangered2","darkblue"))

plot_mort_perc <- plot_mort_perc+ theme_classic()+theme(legend.direction = 'horizontal', 
                                                  legend.position = c(0.1,0.85),
                                                  plot.title = element_text(color="black", size=18, hjust=0.05, vjust=-.2),
                                                  legend.background = element_rect(size=0.5, linetype="solid", color="black"),
                                                  axis.title.x=element_text(colour="black", size = 18, vjust=-0.5),            #sets x axis title size, style, distance from axis #add , face = "bold" if you want bold
                                                  axis.title.y=element_text(colour="black", size = 18, vjust=2),            #sets y axis title size, style, distance from axis #add , face = "bold" if you want bold
                                                  axis.text=element_text(colour="black", size = 14),                             #sets size and style of labels on axes
                                                  legend.title = element_blank(), #remove title of legend
                                                  legend.text = element_text(color="black", size=14))   #plot margin - top, right, bottom, left
plot_mort_perc


# Analysis - Seedling survival (by pot) -----------------------------------

# No of Pots that had at least 1 death in them by cohort and treatment

fisher_data_by_pot<-fisher_data %>% 
  group_by(cohort,block,trt,pot_had_death) %>% 
  summarize(pots_with_death=sum(pot_had_death==TRUE),pots_no_death=sum(pot_had_death==FALSE)) %>% 
  group_by(cohort,trt) %>% 
  summarize(pots_no_death=sum(pots_no_death), pots_with_death=sum(pots_with_death)) %>% 
  mutate(perc_mort=pots_with_death/(pots_no_death+pots_with_death)*100)
ungroup()
fisher_data_by_pot

##########

fisher_data_summary1<-fisher_data_by_pot %>% 
  filter(cohort=="1") %>% 
  select(pots_with_death, pots_no_death) %>% 
  ungroup()
byPotC1<-select(fisher_data_summary1,c(pots_with_death, pots_no_death))
byPotC1<-fisher.test(byPotC1)
p_val_C1surv<-round(byPotC1$p.value,4)

fisher_data_summary2<-fisher_data_by_pot %>% 
  filter(cohort=="2") %>% 
  select(pots_with_death, pots_no_death) %>% 
  ungroup()
byPotC2<-select(fisher_data_summary2,c(pots_with_death, pots_no_death))
byPotC2<-fisher.test(byPotC2)
p_val_C2surv<-round(byPotC2$p.value,3)



plot_mort_pot<-ggplot(fisher_data_by_pot, aes(x=trt, y=perc_mort))+
  geom_bar(stat="summary") +
  ylab("Pots with 1+ dead seedlings (%)") +
  xlab("Seedling density") +
  ggtitle("Percent of Pots in each Treatmnent with at least 1 dead seedling") + 
  facet_grid(cols=vars(cohort), scales="fixed")+
  scale_y_continuous(expand = c(0, 0))

plot_mort_pot <- plot_mort_pot+ theme_classic()+theme(legend.direction = 'horizontal', 
                                                        legend.position = c(0.1,0.85),
                                                        plot.title = element_text(color="black", size=18, hjust=0.05, vjust=-.2),
                                                        legend.background = element_rect(size=0.5, linetype="solid", color="black"),
                                                        axis.title.x=element_text(colour="black", size = 18, vjust=-0.5),            #sets x axis title size, style, distance from axis #add , face = "bold" if you want bold
                                                        axis.title.y=element_text(colour="black", size = 18, vjust=2),            #sets y axis title size, style, distance from axis #add , face = "bold" if you want bold
                                                        axis.text=element_text(colour="black", size = 14),                             #sets size and style of labels on axes
                                                        legend.title = element_blank(), #remove title of legend
                                                        legend.text = element_text(color="black", size=14))   #plot margin - top, right, bottom, left
plot_mort_pot


ftable(all.plants$survive,all.plants$cohort)
ftable(all.plants$trt,all.plants$survive,all.plants$cohort)
ftable(all.plants$cohort,all.plants$survive,all.plants$block)
# Let’s construct a simp
# GLM with binomial (survival = TRUE /FALSE
# Logistic regression
# https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#pronunciation-of-lmerglmeretc.
m1 <- glmer(survive ~ 1 + (1|block),
            family = binomial,
            data = all.plants)
summary(m1)
AIC(m1)

m2 <- glmer(survive ~ 1 + (1|cohort) + (1|block),
            family = binomial,
            data = all.plants)

summary(m2)
AIC(m2)

m3 <- glmer(survive ~ cohort + (1|block),
            family = binomial,
            data = all.plants)

summary(m3)
AIC(m3)



m3b <- glmer(survive ~ trt + (1|cohort),
            family = binomial,
            data = all.plants)

summary(m3b)
AIC(m3b)
Anova(m3b,type="III")

m4 <- glmer(survive ~ trt + (1|cohort) + (1|block),
                               family = binomial,
                               data = all.plants)
summary(m4)
AIC(m4)


m5 <- glmer(survive ~ trt * cohort,
            family = binomial,
            data = all.plants)
summary(m5)
AIC(m5)

m5 <- aov(survive ~ trt:cohort,
            data = all.plants)
summary(m5)
AIC(m5)

car::Anova(m4)

anova(m1,m2)
anova(m2,m3)
anova(m3,m4)
anova(m4,m5)
summary(m1)
Anova(m5,type="III")





# Analyses - Final Leaf Area -----------------------------------------------



# LEAF AREA, BIOMASS, RS ALLOCATION
# library(xtable)
# https://stackoverflow.com/questions/37497948/aov-error-term-in-r-whats-the-difference-bw-errorid-and-errorid-timevar
# aov(Y ~ Error(A), data=d)               # Lone random effect
# aov(Y ~ B + Error(A/B), data=d)         # A random, B fixed, B nested within A
# aov(Y ~ (B*X) + Error(A/(B*X)), data=d) # B and X interact within levels of A

# A=cohort
# B=block
# X=cohort

# first Is there a relationship between final leaf area and biomass
# c1.1<-all.plants %>%
#   filter(cohort=="1") %>% 
#   filter(sdlg.type=="focal") %>% 
#   select(cohort,block,trt,sdlg.type,sdlg.no,days,t1,t6,lf.biomass,stem.biomass,
#          root.biomass,Total.bmass,RSratio,canopy.openess.percent) %>% 
#   rename("la"="t6")
# 
# c2.1<-all.plants %>% 
#   filter(sdlg.type=="focal") %>% 
#   filter(cohort=="2") %>% 
#   select(cohort,block,trt,sdlg.type,sdlg.no,days,t1,t4,lf.biomass,stem.biomass,
#          root.biomass,Total.bmass,RSratio,canopy.openess.percent) %>% 
#   rename("la"="t4")
# 
# 
# 
# after1yr<-bind_rows(c1.1,c2.1)
# 
# bmass.means<-after1yr %>% 
#   group_by(cohort,trt) %>% 
#   filter(sdlg.type=="focal") %>% 
#   drop_na(Total.bmass) %>%  
#   summarize(mean.tbmass=mean(Total.bmass),sd.tbmass=sd(Total.bmass))
# 
# after1yr$block<-as.factor(after1yr$block)
# 
# p <- ggplot(bmass.means, aes(x=trt, y=mean.tbmass, fill=cohort)) + 
#   geom_bar(stat="identity", position=position_dodge()) +
#   geom_errorbar(aes(ymin=mean.tbmass-sd.tbmass, ymax=mean.tbmass+sd.tbmass), width=.2,
#                 position=position_dodge(.9))
# p
# https://stat.ethz.ch/~meier/teaching/anova/random-and-mixed-effects-models.html
# (see also https://www.zoology.ubc.ca/~schluter/R/Model.html#fit_a_linear_mixed-effects_model
# http://bbolker.github.io/mixedmodels-misc/glmmFAQ.html
# https://www.muscardinus.be/2017/07/lme4-random-effects/#:~:text=Implicit%20nesting&text=That%20is%20each%20level%20of,class2%20variable%20in%20our%20data.
# and https://rcompanion.org/handbook/G_03.html)

# Are there correlations among response variables
data=filter(all.plants,cohort=="1")
A<-data$Total.bmass
B<-data$t4
C<-data$RSratio
D<-data$t1
cor.test(C,D)  # YES
cor.test(A,C)  # YES
cor.test(B,D)  # YES
cor.test(A,D)


# analyses - lead area (all surviving plants in a pot) ------------------------------

# Need to figure out if there are "replacement" plants in a pot

# replacements <- all.plants %>% 
#   group_by(cohort,trt,block) %>% 
#   summarize(seedlings_in_pot=n())


# LA SUMMARY TABLES
summ_table_la_c1<-all.plants %>% 
  filter(cohort=="1") %>%
  filter(survive=="alive") %>% 
  group_by(cohort,trt) %>% 
  drop_na(t4) %>% 
  drop_na(t9) %>% 
  summarize(mean_y0=mean(t1),sd_y0=sd(t1),mean_y1=mean(t4),sd_y1=sd(t4),mean_y2=mean(t9),sd_y2=sd(t9))
summ_table_la_c1

summ_table_la_c2<-all.plants %>% 
  filter(cohort=="2") %>% 
  group_by(cohort,trt) %>% 
  filter(survive=="alive") %>% 
  drop_na(t4) %>% 
  summarize(mean_y0=mean(t1),sd_y0=sd(t1),mean_y1=mean(t4),sd_y1=sd(t4))
summ_table_la_c2

la_summary<-bind_rows(summ_table_la_c1,summ_table_la_c2) %>% ungroup()

# Summary Table - LA initial, yr 1, yr2
la_summary

# ratio of b, trt 4 vs trt 1 and 2
percLA4v1<-(la_summary %>% select(cohort,trt,mean_y1) %>% filter(trt=="four"))[3]/
  (la_summary %>% select(cohort,trt,mean_y1) %>%  filter(trt=="one"))[3]
percLA4v1$cohort<-c(1,2)
percLA4v1

percLA2v1<-(la_summary %>% select(cohort,trt,mean_y1) %>% filter(trt=="two"))[3]/
  (la_summary %>% select(cohort,trt,mean_y1) %>%  filter(trt=="one"))[3]
percLA2v1$cohort<-c(1,2)
percLA2v1




######### LA AFTER 1 year of GROWTH (both cohorts)
library(lmerTest)

# LEAF AREA after 1 year, both cohorts, only those alive whole experiment
# (main effect: trt and cohort, covariate = initial leaf area, random effect of block)
mod_yr1_LA <- lmer(t4 ~ trt+cohort+t1+(1|block), data=filter(all.plants,survive=="alive"))
anova(mod_yr1_LA)
rand(mod_yr1_LA)

# LEAF AREA after 2 years, only those alive whole experiment
# (main effect: trt, covariate = initial leaf area, random effect of block)
mod_yr2_LA <- lmer(t9 ~ trt+t1+(1|block), data=filter(all.plants,survive=="alive"))
anova(mod_yr2_LA)
rand(mod_yr2_LA)

# analyses final biomass - focal sdlgs only -------------------------------

BM_summary<-all.plants %>% 
  filter(sdlg.type=="focal") %>% 
  group_by(cohort,trt) %>% 
  drop_na(Total.bmass) %>% 
  summarize(mean_bmass=mean(Total.bmass),sd_bmass=sd(Total.bmass))
BM_summary

# ratio of bm, trt 4 vs trt 1 and 2

percBM4v1<-(BM_summary %>% select(cohort,trt,mean_bmass) %>% filter(trt=="four"))[3]/
  (BM_summary %>% select(cohort,trt,mean_bmass) %>%  filter(trt=="one"))[3]
percBM4v1$cohort<-c(1,2)
percBM4v1

percBM2v1<-(BM_summary %>% select(cohort,trt,mean_bmass) %>% filter(trt=="two"))[3]/
  (BM_summary %>% select(cohort,trt,mean_bmass) %>%  filter(trt=="one"))[3]
percBM2v1$cohort<-c(1,2)
percBM2v1




# TOTAL BIOMASS COHORT 1 (main effect: trt, covariate= final leaf area, random effect of block and block*cohort)
BM.mod1 <- lmer(Total.bmass ~ trt+t9+(1|block), data=filter(all.plants,sdlg.type=="focal" & cohort=="1"))
anova(BM.mod1)
rand(BM.mod1)

# TOTAL BIOMASS COHORT 2 (main effect: trt, covariate= final leaf area, random effect of block and block*cohort)
BM.mod2 <- lmer(Total.bmass ~ trt+t4+(1|block), data=filter(all.plants,sdlg.type=="focal" & cohort=="2"))
anova(BM.mod2)
rand(BM.mod2)

#####################
# RS Ratio
#####################
# TOTAL BIOMASS COHORT 1 (main effect: trt, covariate= fTotal.bmass random effect of block and block*cohort)
RS.mod1 <- lmer(RSratio ~ trt+Total.bmass+(1|block), data=filter(all.plants,sdlg.type=="focal" & cohort=="1"))
anova(RS.mod1)
rand(RS.mod1)

# TOTAL BIOMASS COHORT 2 (main effect: trt, covariate= Total.bmass, random effect of block and block*cohort)
RS.mod2 <- lmer(RSratio ~ trt+Total.bmass+(1|block), data=filter(all.plants,sdlg.type=="focal" & cohort=="2"))
anova(RS.mod2)
rand(RS.mod2)


# Figures -----------------------------------------------------------------


# Cohort 2: Final Leaf Area, Biomass, RS Ratio by trt

# % change in LA vs. Canopy Openness by TRT (A: Cohort 1, B: Cohort 2)
# Will need to calibrate by days because seedlimngs ewere there different times
# Maybe isolate just the ones that were alive the whole time?

# % of seedlings surviving in each treatment (A: Cohort 1, B: Cohort 2)

# Tables  -----------------------------------------------------------------

# Initial Leaf Area for each Cohort







# FIGURES
# 1) Final LA cohort, RS, Total Biomass x trt facet by cohort
# 2) la vs 1|canopy.openess.percent
# 2) la figure vs other variables as in original


  
figdata<-all.plants %>% 
  select(cohort,trt,RSratio,Total.bmass,t9,t4) %>% 
  gather(metric,value,RSratio:t4) 
  # filter(cohort=="2" | (cohort==1 & metric!="t4")) %>% 
  # filter(cohort=="1")
figdata$metric<-gsub("t4", "Final\nLeaf Area\n(sq. cm)",figdata$metric)
figdata$metric<-gsub("t9", "Final\nLeaf Area\n(sq. cm)",figdata$metric)
figdata$metric <- gsub("RSratio","R:S\nRatio",figdata$metric)
figdata$metric <- gsub("Total.bmass","Total\nBiomass\n(g)",figdata$metric)

figdata <- figdata %>% 
  mutate(cohort = ifelse(cohort == "1", 
                         "Cohort 1\n(24 mos)", 
                         "Cohort 2\n(12 mos)")) %>%
  rename("Density Treatment"="trt")


figdata$metric<-as.factor(figdata$metric)
# create a dataset
# Most basic violin chart
library(RColorBrewer)
library(lemon)
library(egg)

# Cohort 1: Final Leaf Area, Biomass, RS Ratio by trt


p <- ggplot(figdata, aes(x=`Density Treatment`, y=value,fill=`Density Treatment`)) + # fill=name allow to automatically dedicate a color for each group
  # geom_violin() +
  geom_boxplot(width=0.9,outlier.size = 0.7) + 
  # stat_summary(fun=mean, geom="point", size=2,color="red")+
  # geom_bar(stat="summary") +
  facet_grid(cols = vars(cohort), rows=vars(metric), scales="free_y")+
  geom_jitter(color="black", size=.7, alpha=0.8,width = 0.15, height = 0)+
  scale_fill_brewer(palette = "Blues")
  
p <- p + theme_classic() + 
  theme(legend.direction = 'vertical', 
        legend.position='none',
        # legend.position = c(0.1,0.85),
         # plot.title = element_text(color="black", size=18, hjust=0.05, vjust=-.2),
         # legend.background = element_rect(size=0.5, linetype="solid", color="black"),
         axis.title.x=element_text(colour="black", size = 18, vjust=-0.5),            #sets x axis title size, style, distance from axis #add , face = "bold" if you want bold
         axis.title.y=element_blank(),            #sets y axis title size, style, distance from axis #add , face = "bold" if you want bold
         axis.text=element_text(colour="black", size = 14),                             #sets size and style of labels on axes
        legend.title = element_blank(), #remove title of legend
        legend.text = element_text(color="black", size=14),
        panel.spacing.x =unit(1.0, "lines"), 
        panel.spacing.y=unit(1,"lines"),
        strip.text.x = element_text(size = 14,margin = margin(0,0,1,0, "lines")),
        strip.text.y = element_text(size = 14, angle=0,margin = margin(0,0,0,0, "lines")),
        strip.background.y = element_rect(fill = NA, colour = NA),
        strip.background.x = element_rect(fill = NA, colour = NA),
        
        )+
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=1.5)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=1.5)
tag_facet(p)

facet_labels<-c("A","B","C","D","E","F")
p<-tag_facet(p,open="", close="", tag_pool=facet_labels,hjust=-1.5, size=6)
p



  



























plot(all.plants$canopy.openess.percent,(all.plants$t9/all.plants$t1*100))

















# p<-ggplot(data=la_b1, aes(x=Total.bmass, y=trt)) +
#   geom_bar(stat="identity")
# p









foo<-filter(la_b2,cohort=="2" & sdlg.type=="focal")  
cor(foo$Total.bmass,foo$final.la,method = "spearman", use= "pairwise.complete.obs")
la_b2 %>% 
  ggplot( aes(x=canopy.openess.percent, y=final.la, color=sdlg.type,fill=sdlg.type)) +
  geom_point()+
  geom_smooth(method='lm')
#   
# bmass %>% filter(cohort=="1") %>%
#   ggplot( aes(x=trt, y=Total.bmass, fill=trt)) +
#   geom_boxplot() +
#   # scale_fill_viridis(discrete = TRUE, alpha=0.6) +
#   geom_jitter(color="black", size=0.4, alpha=0.9) +
#   theme_classic() +
#   theme(
#     legend.position="none",
#     plot.title = element_text(size=11)
#   ) +
#   ggtitle("A boxplot with jitter") +
#   xlab("")


# Conditional case: given survival, was there a difference in biomass?

la_b1<-all.plants %>%
  filter(cohort=="1") %>% 
  select(cohort,block,trt,sdlg.type,sdlg.no,days,t4,lf.biomass,stem.biomass,
         root.biomass,Total.bmass,RSratio,canopy.openess.percent,survive) %>% 
  rename("final.la"="t4")

la_b2<-all.plants %>%
  filter(cohort=="2") %>% 
  select(cohort,block,trt,sdlg.type,sdlg.no,days,t9,lf.biomass,stem.biomass,
         root.biomass,Total.bmass,RSratio,canopy.openess.percent,survive) %>% 
  rename("final.la"="t9")

bmass_data<-bind_rows(la_b1,la_b2) %>% 
    filter(survive=="alive") 
  
  
boxplot(Total.bmass~trt*cohort,data=bmass_data) 
boxplot(RSratio~trt*cohort,data=bmass_data) 






# ######################################################
# Analysis: total leaf area after 12 months (Cohorts 1 and 2)
# ######################################################

ONE_YR_LA_C1<-dplyr::select(COHORT1.wide,cohort, block, trt, sdlg.type, t1,t6, days)
names(ONE_YR_LA_C1)[5] <- "LA_initial"
names(ONE_YR_LA_C1)[6] <- "LA_final"

ONE_YR_LA_C2<-dplyr::select(COHORT2.wide,cohort, block, trt, sdlg.type, t1,t4, days)
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

# Calclute means and SD of initial and final by trt
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

# anova(GLM_LA1.1,GLM_LA1.2, test = "Chisq")  #imp fit by adding LA t1  over just intercept
# anova(GLM_LA1.2,GLM_LA1.3, test = "Chisq")  # imp fit by adding trt over just intercept 
# anova(GLM_LA1.3,GLM_LA1.4, test = "Chisq")  # imp fit by adding cohort instead of intercept
# anova(GLM_LA1.4,GLM_LA1.5, test = "Chisq")  #imp fit by adding block instead of intercept
# anova(GLM_LA1.4,GLM_LA1.6, test = "Chisq")  
# anova(GLM_LA1.7,GLM_LA1.2, test = "Chisq")  
# anova(GLM_LA1.2,GLM_LA1.7, test = "Chisq")  
# anova(GLM_LA1.2,GLM_LA1.9, test = "Chisq")  
# 
# AIC(GLM_LA1.1,GLM_LA1.2,GLM_LA1.3,GLM_LA1.4,GLM_LA1.5,GLM_LA1.6,GLM_LA1.7,GLM_LA1.8,GLM_LA1.9,GLM_LA1.10,
#     GLM_LA1.11,GLM_LA1.12,GLM_LA1.13,GLM_LA1.14,GLM_LA1.15,GLM_LA1.16) #for cohort 2 best model fit is just ht
# 
# # TO REPORT: http://www.ashander.info/posts/2015/10/model-selection-glms-aic-what-to-report/
# summ.tableLA12 <- do.call(rbind, lapply(list(GLM_LA1.1,GLM_LA1.2,GLM_LA1.3,GLM_LA1.4,GLM_LA1.5,GLM_LA1.6,
#                                              GLM_LA1.7,GLM_LA1.8,GLM_LA1.9,GLM_LA1.10,GLM_LA1.11,GLM_LA1.12,
#                                              GLM_LA1.13,GLM_LA1.14,GLM_LA1.15,GLM_LA1.16), broom::glance))
# summ.tableLA12
# table.cols12 <- c("df.residual", "deviance", "AIC")
# reported.tableLA12 <- summ.tableLA12[table.cols]
# names(reported.tableLA12) <- c("Resid. Df", "Resid. Dev", "AIC")
# 
# reported.tableLA12[['dAIC']] <-  with(reported.tableLA12, AIC - min(AIC))
# reported.tableLA12[['weight']] <- with(reported.tableLA12, exp(- 0.5 * dAIC) / sum(exp(- 0.5 * dAIC)))
# reported.tableLA12$AIC <- NULL
# reported.tableLA12$weight <- round(reported.tableLA12$weight, 2)
# reported.tableLA12$dAIC <- round(reported.tableLA12$dAIC, 1)
# row.names(reported.tableLA12) <- model.namesLA1
# reported.tableLA12
# 
# write.table(reported.tableLA12, file = "LA_1yr_GLM.csv", , sep = ",", na = "NA", row.names = TRUE)


# ######################################################
# total leaf area after 24 months (Cohort 1)
# ######################################################
# 

# Are final height and length correted?
# 
# plot(sdlg_size_end)
# hist(sdlg_size_end$LA_final)
# hist(sdlg_size_end$HT_final)
# str(sdlg_size_end)
# summary(sdlg_size_end)
# cor(sdlg_size_end[2:3], method="spearman", use="complete.obs")  #Storng corretion between final leaf area and final height
# 
# # Calclute LA and SD by trt
# na.omit(sdlg_size_end)%>% group_by(trt) %>%summarise(avg=mean(LA_final)) #final seedling la by trt
# na.omit(sdlg_size_end)%>% group_by(trt) %>%summarise(sd=sd(LA_final))
# 
# # do same but only seedlings aive the who experiment to see if the smaller size is due to density or transplanting in smaller ones when some died.
# all.exp.sdlgs<-filter(sdlg_size_end, days >=773)
# na.omit(all.exp.sdlgs)%>% group_by(trt) %>%summarise(avg=mean(LA_final)) #final seedling la by trt
# na.omit(all.exp.sdlgs)%>% group_by(trt) %>%summarise(sd=sd(LA_final))
# boxplot(LA_final~trt,data=all.exp.sdlgs) 

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
# 
# GLM_LA2.1<-glm(LA_final ~ 1, data = sdlg_size_end, family = gaussian) #INTERCEPT
#   summary(GLM_LA2.1)
# GLM_LA2.2<-glm(LA_final ~ LA_initial, data = sdlg_size_end, family = gaussian) #INITIAL LEAF AREA 
#   summary(GLM_LA2.2)
# GLM_LA2.3<-glm(LA_final ~ block, data = sdlg_size_end, family = gaussian) #BLOCK
#   summary(GLM_LA2.3)
# GLM_LA2.4<-glm(LA_final ~ trt, data = sdlg_size_end, family = gaussian) #TREATMENT
#   summary(GLM_LA2.4)
# GLM_LA2.5<-glm(LA_final ~ days, data = sdlg_size_end, family = gaussian) #days
#   summary(GLM_LA2.5)
# GLM_LA2.6<-glm(LA_final ~ block+LA_initial, data = (sdlg_size_end), family = gaussian) #INITIAL LA + BLOCK
#   summary(GLM_LA2.6)
# GLM_LA2.7<-glm(LA_final ~ trt+LA_initial, data = sdlg_size_end, family = gaussian) # INITIAL LA + TRT
#   summary(GLM_LA2.7)
# GLM_LA2.8<-glm(LA_final ~ trt*days, data = sdlg_size_end, family = gaussian) #TRT*DAYS
#   summary(GLM_LA2.8)
# GLM_LA2.9<-glm(LA_final ~ trt+LA_initial+block, data = sdlg_size_end, family = gaussian) #INITIAL LA+TRT+BLOCK
#   summary(GLM_LA2.9)
# GLM_LA2.10<-glm(LA_final ~ trt+LA_initial+days, data = sdlg_size_end, family = gaussian) #INITIAL LA+TRT+DAYS
#   summary(GLM_LA2.10)
# GLM_LA2.11<-glm(LA_final ~ trt*days+LA_initial, data = sdlg_size_end, family = gaussian) #INITIAL LA+TRT*DAYS
#   summary(GLM_LA2.11)
# GLM_LA2.12<-glm(LA_final ~ trt+LA_initial+days+block, data = sdlg_size_end, family = gaussian) #INITIAL LA+TRT+BLOCK+DAYS
#   summary(GLM_LA2.12)
# GLM_LA2.13<-glm(LA_final ~ trt*days+block, data = sdlg_size_end, family = gaussian) #INITIAL LA+TRT*DAYS
#   summary(GLM_LA2.13)
# GLM_LA2.14<-glm(LA_final ~ trt*days+LA_initial+block, data = sdlg_size_end, family = gaussian) #INITIAL LA+TRT*DAYS
#   summary(GLM_LA2.14)
# GLM_LA2.15<-glm(LA_final ~ trt+days, data = sdlg_size_end, family = gaussian) #TRT*DAYS
#   summary(GLM_LA2.15)
# 
# 
# model.names <- c("1 Intercept", "2 Initial Leaf Area", "3 Block", "4 Density","5 Days Since Transplant",
#                  "6 Initial Leaf Area+Block",
#                  "7 Density+Initial Leaf Area",
#                  "8 Density*Days Since Transplant",
#                  "9 Density+Initial Leaf Area+Block",
#                  "10 Density+Initial Leaf Area+Days Since Transplant",
#                  "11 Density*Days Since Transplant+Initial Leaf Area",
#                  "12 Density+Days Since Transplant+Initial Leaf Area+Block",
#                  "13 Density*Days Since Transplant+Block",
#                  "14 Density*Days Since Transplant+Initial Leaf Area+Block",
#                  "15 Density+Days Since Transplanting")
# 
# 
# AIC(GLM_LA2.1,GLM_LA2.2,GLM_LA2.3,GLM_LA2.4,GLM_LA2.5,GLM_LA2.6,GLM_LA2.7,GLM_LA2.8, GLM_LA2.9, 
#     GLM_LA2.10, GLM_LA2.11, GLM_LA2.12, GLM_LA2.13, GLM_LA2.14, GLM_LA2.15) #adding block doesn't result in lower AIC
# 
# 
# # TO REPORT: http://www.ashander.info/posts/2015/10/model-selection-glms-aic-what-to-report/
# summ.tableLA24 <- do.call(rbind, lapply(list(GLM_LA2.1,GLM_LA2.2,GLM_LA2.3,GLM_LA2.4,GLM_LA2.5,GLM_LA2.6,GLM_LA2.7,
#                                              GLM_LA2.8, GLM_LA2.9, GLM_LA2.10, GLM_LA2.11, GLM_LA2.12, GLM_LA2.13, GLM_LA2.14, GLM_LA2.15), broom::glance))
# summ.tableLA24
# table.colsLA24 <- c("df.residual", "deviance", "AIC")
# reported.tableLA24 <- summ.tableLA24[table.colsLA24]
# names(reported.tableLA24) <- c("Resid. Df", "Resid. Dev", "AIC")
# 
# reported.tableLA24[['dAIC']] <-  with(reported.tableLA24, AIC - min(AIC))
# reported.tableLA24[['weight']] <- with(reported.tableLA24, exp(- 0.5 * dAIC) / sum(exp(- 0.5 * dAIC)))
# reported.tableLA24$AIC <- NULL
# reported.tableLA24$weight <- round(reported.tableLA24$weight, 2)
# reported.tableLA24$dAIC <- round(reported.tableLA24$dAIC, 1)
# row.names(reported.tableLA24) <- model.names
# reported.tableLA24

# write.table(reported.tableLA24, file = "LA_2yr_GLM.csv", , sep = ",", na = "NA", row.names = TRUE)
# 
# anova(GLM_LA2.1,GLM_LA2.2, test = "Chisq")  #do not imp fit by adding INITIAL LA over just intercept
# anova(GLM_LA2.1,GLM_LA2.3, test = "Chisq")  #imp fit by adding BLOCK over just intercept 
# anova(GLM_LA2.1,GLM_LA2.4, test = "Chisq")  #imp fit by adding TRT P=0.056
# anova(GLM_LA2.1,GLM_LA2.5, test = "Chisq")  #adding DAYS SINCE TPLANT imporives fit
# anova(GLM_LA2.3,GLM_LA2.6, test = "Chisq")  #adding trt and trt*ht interaction best
# anova(GLM_LA2.10,GLM_LA2.12, test = "Chisq")  #adding trt and trt*ht interaction and days
# ######################################################
# ######################################################
# BIOMASS ANALYSES - SET UP
# ######################################################
# ####################################################### ######################################################
# 
# bmass[, "RSratio"] <-(bmass$root.biomass/(bmass$lf.biomass+bmass$stem.biomass))
# bmass[, "Total.bmass"] <-bmass$root.biomass+bmass$lf.biomass+bmass$stem.biomass
# bmass<-full_join(bmass, canopy, by = "block")
# bmass$location<-NULL
# bmass$site.openess.percent<-NULL
# bmass$mask.openess.percent<-NULL
# bmass$sky.area.percent<-NULL
# bmass$LAI.4.ring<-NULL 
# bmass$LAI.5.ring<-NULL
# bmass$cohort<-as.factor(bmass$cohort)
# bmass$block<-as.integer(bmass$block)
# 
# bmass<-arrange(bmass,cohort,block,trt)
# 
# LA_for_Bmass<-dplyr::filter(ONE_YEAR_LA,sdlg.type=="focal")
# LA_for_Bmass<-droplevels(LA_for_Bmass)
# LA_for_Bmass<-arrange(LA_for_Bmass,cohort,block,trt)
# 
# str(LA_for_Bmass)
# str(bmass)
# 
# bmass<-cbind(bmass,LA_for_Bmass$LA_initial,LA_for_Bmass$LA_final)
# 
# names(bmass)[11] <- "LA_initial"
# names(bmass)[12] <- "LA_final"
# 
# bmassC1<-filter(bmass,  cohort == "1") #cohort 1 after 2 years
# bmassC1$cohort<-as.factor(bmassC1$cohort)
# 
# bmassC2<-filter(bmass,  cohort == "2") #cohort 2 after 1 year
# bmassC2$cohort<-as.factor(bmassC2$cohort)
# 
# # Corr of leaf area and biomass
# corr1<-bmass$Total.bmass
# corr2<-bmass$LA_final
# CORR<-cbind(corr1,corr2)
# cor(CORR, method="spearman", use="complete.obs")  #Storng corretion between final leaf area and final height

# Visualizations
# 
# ggplot(bmass, aes(x=LA_initial, y=Total.bmass, color=trt)) + geom_point(shape=1)+
#   geom_point(shape=1, position=position_jitter(width=1,height=.5))+ #shape 1 = hollow circles, jitter plot to seperate overlapping points
#   geom_smooth(method=lm, se=FALSE)    #Add linear regression lines and Don't add shaded confidence region


# ######################################################
# total biomass after 12 months (Cohort 2)
# ######################################################
# # 
# # Histogram
# hist(bmassC2$Total.bmass)
# # BOX PLOTS
# boxplot(Total.bmass~trt,data=bmassC2) #Cohort bmass final
# 
# # GLMs
# GLM_BM_1<-glm(Total.bmass ~ 1, data = bmassC2, family = gaussian)
# summary(GLM_BM_1)
# GLM_BM_2<-glm(Total.bmass ~ LA_initial, data = bmassC2, family = gaussian)
# summary(GLM_BM_2)
# GLM_BM_3<-glm(Total.bmass ~ trt, data = bmassC2, family = gaussian)
# summary(GLM_BM_3)
# GLM_BM_4<-glm(Total.bmass ~ block, data = bmassC2, family = gaussian)
# summary(GLM_BM_4)
# GLM_BM_5<-glm(Total.bmass ~ trt+block, data = bmassC2, family = gaussian)
# summary(GLM_BM_5)
# GLM_BM_6<-glm(Total.bmass ~ trt+LA_initial+block, data = bmassC2, family = gaussian)
# summary(GLM_BM_6)
# 
# anova(GLM_BM_1,GLM_BM_1, test = "Chisq")  #imp fit by adding HT over just intercept
# anova(GLM_BM_1,GLM_BM_1, test = "Chisq")  #no imp fit by adding trt over just intercept (but close)
# anova(GLM_BM_1,GLM_BM_1, test = "Chisq")  #no imp fit by adding trt instead of ht
# anova(GLM_BM_1,GLM_BM_1, test = "Chisq")  #adding trt and ht better than just ht 
# anova(GLM_BM_1,GLM_BM_1, test = "Chisq")  #adding trt and trt*ht interaction best
# anova(GLM_BM_1,GLM_BM_1, test = "Chisq")  #adding trt and trt*ht interaction and days
# 
# AIC(GLM_BM_1,GLM_BM_2,GLM_BM_3,GLM_BM_4,GLM_BM_5,GLM_BM_6) #adding block doesn't result in lower AIC
# 
# # TO REPORT: http://www.ashander.info/posts/2015/10/model-selection-glms-aic-what-to-report/
# 

# 
# bm1yr_plot<-ggplot(bmassC2, aes(x=LA_final, y=Total.bmass, color=trt))+
#   geom_point(shape=16, position=position_jitter(width=1,height=.5), size=3)+ #shape 1 = hollow circles, jitter plot to seperate overlapping points
#   geom_smooth(method=lm, se=FALSE)+    #Add linear regression lines and Don't add shaded confidence region
#   ylab("Total biomass (g)") +
#   xlab("Final leaf Area (cm2)")+
#   ggtitle("A) 12 months")+
#   scale_colour_manual(values=c("gray54", "orangered2","darkblue"))
# 
# bm1yr_plot<-bm1yr_plot+ theme_classic()+theme(legend.direction = 'vertical', 
#                                               legend.position = c(0.1,0.85),
#                                               plot.title = element_text(color="black", size=18, hjust=0.05, vjust=-.2),
#                                               legend.background = element_rect(size=0.5, linetype="solid", color="black"),
#                                               axis.title.x=element_text(colour="black", size = 18, vjust=-0.5),            #sets x axis title size, style, distance from axis #add , face = "bold" if you want bold
#                                               axis.title.y=element_text(colour="black", size = 18, vjust=2),            #sets y axis title size, style, distance from axis #add , face = "bold" if you want bold
#                                               axis.text=element_text(colour="black", size = 14),                             #sets size and style of labels on axes
#                                               legend.title = element_blank(), #remove title of legend
#                                               legend.text = element_text(color="black", size=14))   #plot margin - top, right, bottom, left
# 
# 
# print(bm1yr_plot)



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


plot.sdlgs<-read.csv("data/demog_plot_sdlgs.csv", dec=".", header = TRUE, sep = ",", check.names=FALSE )
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



