setwd("/Users/Wenwen/UW/Global Health - MobileWAChX")
#setwd("/Users/keshet/UW/Global Health - MobileWAChX")

library(foreign); library(dplyr); library(lubridate);library(haven);
library(stringr); library(reshape2);
library(survival); library(survminer); library(ggplot2); library(ggthemes);
library(ggpubr); library(magrittr); library(readxl);
library(lme4); library(MASS); library(gmodels); library(gee); library(geepack);
library(sandwich); library(arsenal); library(lmtest); library(gridExtra)

##############################################################################
# read abstraction data for CD4 variables at enrollment
##############################################################################
full <- read_dta('3. DataMerges/AbstractionMerges/Abstraction_CRFs_020420.dta') 

# obtain maternal cd4 data
cd4 <- full %>% 
  dplyr::select(ptidno, Enrolled, eventdate, cd4count, cd4pct) %>% 
  subset(cd4count!='' & cd4count!='NA') %>% unique() %>%
  arrange(ptidno, eventdate) %>% group_by(ptidno) %>% filter(row_number()==1)

cd4 <- cd4 %>% dplyr::select(ptidno, cd4count)

# obtain maternal height and weight data
hw <- full %>% 
  dplyr::select(ptidno, eventdate, 
                hivcare7ab_cccheight, hivvisit7ab_cccheight,
                hivcare7ab_cccweight, hivvisit7ab_cccweight) %>% 
  # combine two weight variables
  mutate(height = ifelse(!is.na(hivcare7ab_cccheight), hivcare7ab_cccheight, hivvisit7ab_cccheight),
         weight = ifelse(!is.na(hivcare7ab_cccweight), hivcare7ab_cccweight, hivvisit7ab_cccweight)) %>%
  # clean this variable by changing- -2 to NA
  mutate(height = ifelse(height==-2, NA, height),
         weight = ifelse(weight==-2, NA, weight)) %>%
  # subset with non-missing data
  subset(!is.na(height) | !is.na(weight)) %>% unique() %>%
  dplyr::select(ptidno, eventdate, height, weight) %>% 
  arrange(ptidno, eventdate)

summary(hw$height)
summary(hw$weight)


############################################################################################
# read merged interview data for enrollment data
############################################################################################
interview <- read.dta('3. DataMerges/InterviewMerges/AllMerged_Interview_CRFs_results.dta')
interview$ptidno <- str_sub(interview$ptidno, 1, 8)


# enroll
enroll <- interview %>% subset(visit=='Enrollment') %>% unique()
enroll$enroldate <- as.Date(enroll$enroldate, '%m/%d/%y')

# exclude 1 maternal death before delivery
enroll <- enroll %>% subset(ptidno!='25420408')


############################################################################################
# read delivery data
############################################################################################
delivery <- read.csv('3. DataMerges/delivery/Delivery dates_verified_12Nov2019_WJ.csv')
delivery$deldate <- as.Date(delivery$deldate, '%d-%b-%y')



############################################################################################
# read SAE log
############################################################################################
sae <- read_excel('11. SAE logs/SAE_verified_28july2020.xlsx', sheet = 1)
sae$saetype <- tolower(sae$saetype)
inf <- sae %>% mutate(saetype = ifelse(str_detect(saetype, 'misca')==T, 'miscarriage',
                                       ifelse(str_detect(saetype, 'still')==T, 'stillbirth',
                                              ifelse(str_detect(saetype, 'infant')==T & str_detect(saetype, 'death')==T, 'infdeath', NA)))) %>%
  subset(!is.na(saetype))

############################################################################################
# merge with enroll
############################################################################################
mydataall <- merge(enroll, inf, by='ptidno', all.x=T)
mydataall <- merge(mydataall, delivery, all.x = T)
mydataall$saedate <- as.Date(mydataall$saedate)
mydataall$deldate <- as.Date(mydataall$deldate)
mydataall$en_lmp <- ifelse(mydataall$en_lmp=='2030-01-01', NA, mydataall$en_lmp)
mydataall$en_lmp <- as.Date(mydataall$en_lmp, origins='1970-01-01')

########################################################################
############### NO NEED TO RUN, THESE ARE ONLY LMP CHECK ###############
########################################################################

d <- mydataall %>% 
  dplyr::select(ptidno, enroldate, deldate,
                en_lmp, en_lmpsure, en_lmpsuremo) %>%
  mutate(gestageatenr = round((enroldate - en_lmp)/7,1),
         gestageatdel = round((deldate - en_lmp)/7,1),
         lmpsure = ifelse(en_lmpsure=='Yes', 'Yes with date', 'No'),
         lmpsure2 = ifelse(en_lmpsure=='Yes', 'Yes', ifelse(en_lmpsuremo=='Yes', 'Yes', 'No')),
         lmpsure3 = ifelse(en_lmpsure=='Yes', 'Da', ifelse(en_lmpsuremo=='Yes', 'Mo', 'No')))

table(d$lmpsure)
table(d$lmpsure2)
table(d$lmpsure3)


# distribution function
gestageplot <- function(m, lmp, title) {
  # density plot
  a <- ggplot(d, aes(x=m))+geom_density(aes(fill=lmp), alpha=.5)+
    scale_x_continuous('') + 
    theme(panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.position = 'top', axis.text.x = element_text(hjust = 1))
  # boxplot
  b <- ggplot(d, aes(x=lmp, y=m, fill=lmp))+
    geom_boxplot(alpha=0.5, width=0.5)+
    scale_y_continuous(title) + 
    theme(panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.position = 'none', axis.text.x = element_text(hjust = 1))+
    coord_flip()
   # combine two graphs
  ggarrange(a,b, nrow=2, heights = 2:1)
}

# gestage at enrollment 
## by lmpsure
gestageplot(d$gestageatenr, d$lmpsure, 'gestational age at enrollment (week) by 2-level lmp (date Y/N)')
descrip(d$gestageatenr, strata=d$lmpsure)
t.test(d$gestageatenr ~ d$lmpsure, var.equal = TRUE)

## by lmpsure2
gestageplot(d$gestageatenr, d$lmpsure2, "gestational age at enrollment (week) by 2-level lmp (month Y/N)")
descrip(d$gestageatenr, strata=d$lmpsure2)
t.test(d$gestageatenr ~ d$lmpsure2, var.equal = TRUE)

## by lmpsure3
gestageplot(d$gestageatenr, d$lmpsure3, "gestational age at enrollment (week) by 3-level lmp (date/month/N)")
descrip(d$gestageatenr, strata=d$lmpsure3)
pairwise.t.test(d$gestageatenr, d$lmpsure3, pool.sd = T)


# gestage at delivery 
## by lmpsure
gestageplot(d$gestageatdel, d$en_lmpsure, "gestational age at delivery (week) by 2-level lmp (date Y/N)")
descrip(d$gestageatdel, strata=d$en_lmpsure)
t.test(d$gestageatdel ~ d$lmpsure, var.equal = TRUE)

## by lmpsure2
gestageplot(d$gestageatdel, d$lmpsure2, "gestational age at delivery (week) by 2-level lmp (month Y/N)")
descrip(d$gestageatdel, strata=d$lmpsure2)
pairwise.t.test(d$gestageatdel, d$lmpsure2, pool.sd = T)

## by lmpsure3
gestageplot(d$gestageatdel, d$lmpsure3, "gestational age at delivery (week) by 3-level lmp (date/month/N)")
descrip(d$gestageatdel, strata=d$lmpsure3)
pairwise.t.test(d$gestageatdel, d$lmpsure3, pool.sd = T)

## calculate median and IQR
stat <- d %>% arrange(en_lmpsure3) %>% group_by(en_lmpsure3) %>% 
  summarise(med=median(gestageatdel, na.rm=T),
            low=quantile(gestageatdel, 0.25, na.rm=T),
            upp=quantile(gestageatdel, 0.75, na.rm=T))

########################################################################
############### LMP CHECK FINISHED ###############
########################################################################


#------------------------------ data cleaning desicion ------------------------------#

# exclude women missing LMP, or unsuring about any LMP
mydata <- mydataall %>% subset(!is.na(en_lmp) & (en_lmpsuremo=='Yes' | is.na(en_lmpsuremo)))

#------------------------------ data cleaning decision ------------------------------#


############################################################
# data check of ms / sb
############################################################


# generate gestage at enr and gestage at del
mydata <- mydata %>% 
  mutate(gestageatenr = enroldate - en_lmp,
         gestageatdel = deldate - en_lmp) %>%
  # assign 0 to non-sae
  mutate(saetype = ifelse(is.na(saetype), 0, saetype))

summary(as.numeric(mydata$gestageatenr)/7)
summary(as.numeric(mydata$gestageatdel)/7)

addmargins(table(mydata$saetype))

# re-calculate miscarriage and stillbirth based on lmp
mydata <- mydata %>%
  mutate(nonalive = ifelse(saetype!=0 & saetype!='infdeath', 1, 0),
         adv_sb = ifelse(nonalive==1 & gestageatdel >=140, 1, 0),
         adv_ptb = ifelse(gestageatenr>=259, NA,
                          ifelse(adv_sb==1, NA,
                                 ifelse(gestageatdel<259, 1, 0))),
         adv_nnd = ifelse(adv_sb==1, NA, ifelse(saetype=='infdeath' & saedate-deldate<=28, 1, 0)))

addmargins(table(mydata$adv_ptb, mydata$nonalive, useNA='ifany'))

# exclude miscarriage < 140 days gestage at delivery
mydata <- mydata %>% subset(!(nonalive==1 & gestageatdel<140))

# calculate early sb among enrolled < 28wk = 196 days
# calculate late sb among enrolled < 32wk = 224 days
# calculate term sb among enrolled >= 32wk = 224 days

mydata <- mydata %>% mutate(adv_esb = ifelse(gestageatenr>=196, NA,
                                             ifelse(adv_sb==1 & gestageatdel<196, 1, 
                                                    ifelse(adv_sb==1 & gestageatdel>=196, 0, 
                                                           ifelse(adv_sb==0, 0, NA)))))
mydata <- mydata %>% mutate(adv_lsb = ifelse(gestageatenr>=252, NA,
                                             ifelse(adv_sb==1 & gestageatdel<252 & gestageatdel>=196, 1, 
                                                    ifelse(adv_sb==1 & gestageatdel>=252, 0, 
                                                           ifelse(adv_sb==1 & gestageatdel<196, NA, 
                                                                  ifelse(adv_sb==0, 0, NA))))))
mydata <- mydata %>% mutate(adv_tsb = ifelse(adv_sb==1 & gestageatdel>=252, 1, 
                                                    ifelse(adv_sb==0, 0, NA)))

addmargins(table(mydata$adv_esb, mydata$nonalive, useNA='ifany'))
addmargins(table(mydata$adv_lsb, mydata$nonalive, useNA='ifany'))
addmargins(table(mydata$adv_tsb, mydata$nonalive, useNA='ifany'))

# calculate severe adv = very ptb among enrolled 28-32wk = 196-224 days
mydata <- mydata %>% mutate(adv_vptb = ifelse(gestageatenr>=224, NA,
                                              ifelse(nonalive==0 & (gestageatdel>196 & gestageatdel<224), 1, 
                                                     ifelse(nonalive==0 & (gestageatdel<=196 | gestageatdel>=224), 0, 
                                                            ifelse(adv_ptb==0, 0, NA)))))
addmargins(table(mydata$adv_vptb, mydata$nonalive, useNA='ifany'))

# calculate any adv
mydata <- mydata %>% mutate(adv_any = ifelse(adv_sb==1, 1, 
                                             ifelse(adv_ptb==1, 1, 
                                                    ifelse(adv_nnd==1, 1, 0))))
table(mydata$adv_any)
addmargins(table(mydata$adv_esb, mydata$adv_sb, useNA='ifany'))
addmargins(table(mydata$adv_vptb, mydata$adv_ptb, useNA='ifany'))

#------------------------------ here are the cleaned outcomes ------------------------------#

# gestational age at delivery by stillbirth - barchart
mydata <- mydata %>% group_by(adv_sb) %>% mutate(mean_gestageatdel=mean(gestageatdel))

ggplot(mydata, aes(x=gestageatdel/7))+
  geom_density(aes(fill=as.character(adv_sb)), alpha=0.5)+
  geom_vline(aes(xintercept=mean_gestageatdel/7, color=as.character(adv_sb)),linetype="dashed")+
  scale_x_continuous('gestational age at delivery (weeks)', breaks = seq(0,48,4), labels = seq(0,48,4), expand = c(0, 0)) + 
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = 'top', axis.text.x = element_text(hjust = 1))
  



###################################################################################
# data cleaning of demographics
###################################################################################
# merge with cd4count
mydata <- merge(cd4, mydata, by='ptidno', all.y=T)
mydata$cd4count <- as.numeric(mydata$cd4count)

# adolescent vs adult, mother
mydata$aya <- ifelse(mydata$en_age<=24, 'aya', 'adult')

# site hand
table(mydata$site_hand)
mydata$site <- ifelse(mydata$site_hand==11 | mydata$site_hand==12, 'Nairobi', 'Western')

# education
table(mydata$en_educlevel)
mydata$primarycom <- ifelse(mydata$en_educlevel=='Less than primary' | mydata$en_educlevel==' No formal education', 'No', 'Yes')
mydata$secondarycom <- ifelse(mydata$en_educlevel=='Secondary - completed' | mydata$en_educlevel=='Above secondary', 'Yes', 'No')

# currently in school
table(mydata$en_inschool)
mydata$en_inschool <- ifelse(mydata$en_inschool=='No answer', NA, as.character(mydata$en_inschool))

# marriage
table(mydata$en_marital)
mydata$curmarried <- ifelse(mydata$en_marital == 'Currently married' | mydata$en_marital == 'Come we stay', 'Yes', 'No')

# partner tested for HIV
table(mydata$di_parttest)
mydata$parttested <- ifelse(mydata$di_parttest==1, 'Yes', ifelse(mydata$di_parttest==0, 'No', NA))

# employment
table(mydata$en_employ)
mydata$employ <- ifelse(mydata$en_employ == 'no answer', NA, ifelse(mydata$en_employ == 'unemployed', 'No', 'Yes'))

# monthly income
mydata$income <- ifelse(mydata$en_monthincome==-2, NA, mydata$en_monthincome)

# abuse in the last year and recent
table(mydata$ip_abuselastyr)
mydata$abuselastyr <- ifelse(mydata$ip_abuselastyr=='No answer', NA, as.character(mydata$ip_abuselastyr))
table(mydata$ip_abuserecent)
mydata$abuserecent <- ifelse(is.na(mydata$ip_abuserecent), 'No', ifelse(mydata$ip_abuserecent=='Yes', 'Yes', 'No'))
table(mydata$abuserecent)

# number of ANC visits during current preganancy >1
table(mydata$en_ancvisits)
mydata$ancvisitgt1 <- ifelse(mydata$en_ancvisits >1, 'Yes', 'No')

# age of being sex active
mydata$sexactive <- ifelse(mydata$en_sexdebut<=0, NA, mydata$en_sexdebut)
mydata$sexactive_bin <- ifelse(mydata$sexactive<15, '<15', '>=15')

# Have ever used family planning
table(mydata$fp_fpever)
mydata$fpever <- ifelse(mydata$fp_fpever=='Yes', 'Yes', ifelse(mydata$fp_fpever=='No', 'No', NA))

# pregnancy history
mydata <- mydata %>%
  mutate(primi=ifelse(en_pregnum==1, 'Yes', 'No'),
         fetaldeathhis=ifelse(primi=='Yes', 'na-primi', ifelse(en_sponabortion+en_stillbirths>0, 'yes', 'no')),
         infdeathhis=ifelse(primi=='Yes', 'na-primi', ifelse(en_childrendied>0, 'yes', 'no')))

# history of STI
table(mydata$en_sti)
tab <- (tableby(en_sti ~ en_ulcers+en_warts+en_gon+en_syph+en_trich+en_chlam+en_othersti+en_sti_sp,
                data=mydata, control=desc))
summary(tab)# ulcers=2, warts=2, gon=4, syph=15, chlam=1, cand=1, cer=1, hiv=1
## one woman said Candidiasis, it is a yeast infection, exclude
mydata$sti <- ifelse(mydata$en_sti_sp=='Candidiasis', NA, 
                     ifelse(mydata$en_sti=='Yes', 'Yes', ifelse(mydata$en_sti=='No', 'No', NA)))
table(mydata$sti)
## clean syph
mydata$syph <- ifelse(is.na(mydata$en_syph), 'No', ifelse(mydata$en_syph=='No', 'No', 'Yes'))
table(mydata$sti, mydata$syph)

# When you got pregnant with this pregnancy, did you want to have a / another baby
table(mydata$fp_pregwant)
mydata$pregwant <- ifelse(mydata$fp_pregwant=='Yes', 'Yes', ifelse(mydata$fp_pregwant=='No', 'No', NA))

# ARV status
table(mydata$ad_currentarv)
mydata$currentarv <- ifelse(str_sub(mydata$ad_currentarv, 1, 3)=='Yes', 'Yes, on', ifelse(str_sub(mydata$ad_currentarv, 1, 3)=='No,', 'No, new', 'Not'))
table(mydata$currentarv)

# ARV regimen
addmargins(table(mydata$ad_currentregimen))
addmargins(table(mydata$ad_currentregimen_sp))

## create a new ART regimen variable
mydata$adreg <- ifelse(mydata$ad_currentregimen!='Other', as.character(mydata$ad_currentregimen), mydata$ad_currentregimen_sp)
addmargins(table(mydata$adreg))

### 3-drug ART first regimen (TDF/AZT/...)
mydata$adreg_1st <- ifelse(mydata$adreg=='No answer', NA, 
                                   ifelse(str_detect(mydata$adreg, 'TDF')==T, 'TDF', 
                                          ifelse(str_detect(mydata$adreg, 'AZT')==T, 'ZDV', 'unknown')))

table(mydata$adreg_1st); prop.table(table(mydata$adreg_1st))*100

mydata$adreg_1st_bi <- ifelse(mydata$adreg_1st=='unknown', NA, mydata$adreg_1st)
table(mydata$adreg_1st, mydata$adreg_1st_bi)

### 3-drug ART third regimen (EFV/NVP/LPV/r/...)
mydata$adreg_3rd <- ifelse(mydata$adreg=='No answer', NA, 
                        ifelse(str_detect(mydata$adreg, 'NVP')==T, 'NVP', 
                               ifelse(str_detect(mydata$adreg, 'EFV')==T, 'EFV', 
                                      ifelse(str_detect(mydata$adreg, 'LPV')==T, 'LPV/r', 'unknown'))))

table(mydata$adreg_3rd); prop.table(table(mydata$adreg_3rd))*100

mydata$adreg_3rd_bi <- ifelse(mydata$adreg_3rd=='unknown' | mydata$adreg_3rd=='LPV/r', NA, mydata$adreg_3rd)
table(mydata$adreg_3rd, mydata$adreg_3rd_bi)

# knowledge score
table(mydata$ad_motivation_findout)
mydata <- mydata %>% mutate_at(vars(starts_with('ad_')), as.numeric)
mydata$knowledge_howtake <- mydata$ad_knowledge_howtake-1
mydata$knowledge_skip_rev <- ifelse(mydata$ad_knowledge_skip==1, 0, 7-mydata$ad_knowledge_skip)
mydata$knowledge_sideeffec <- mydata$ad_knowledge_sideeffec-1
mydata$knowledge_livelong <- mydata$ad_knowledge_livelong-1

mydata$motivation_findout <- mydata$ad_motivation_findout-1
mydata$motivation_support <- mydata$ad_motivation_support-1
mydata$motivation_forever <- mydata$ad_motivation_forever-1
mydata$motivation_sideeffe <- mydata$ad_motivation_sideeffe-1

mydata$skills_support <- mydata$ad_skills_support-1
mydata$skills_refills <- mydata$ad_skills_refills-1
mydata$skills_sideeffects <- mydata$ad_skills_sideeffects-1
mydata$skills_remember <- mydata$ad_skills_remember-1
mydata$skills_emotion <- mydata$ad_skills_emotion-1
mydata$skills_physgood <- mydata$ad_skills_physgood-1
mydata$skills_physbad <- mydata$ad_skills_physbad-1

mydata$imb_i <- mydata$knowledge_howtake + mydata$knowledge_skip_rev + mydata$knowledge_sideeffec + mydata$knowledge_livelong
mydata$imb_m <- mydata$motivation_findout + mydata$motivation_support + mydata$motivation_forever + mydata$motivation_sideeffe
mydata$imb_b <- mydata$skills_support + mydata$skills_refills + mydata$skills_sideeffects + mydata$skills_remember + mydata$skills_emotion + mydata$skills_physgood + mydata$skills_physbad

mydata$imbscore <- mydata$imb_i + mydata$imb_m + mydata$imb_b

# standardize as % of max
mydata$imbscore_pct = mydata$imbscore/75*100
summary(mydata$imbscore_pct)

# know you were positive when pregnant
table(mydata$en_hivpospreg)
mydata$know <- ifelse(mydata$en_hivpospreg=='Yes', 'Yes', ifelse(mydata$en_hivpospreg=='No', 'No', NA))

# disclosure to anyone
mydata$disc <- ifelse(mydata$di_disctopart=='Yes' | mydata$di_disctoother1!='None', 'Yes', 'No') 

# depression score >5, >10
table(mydata$dp_tired)
mydata <- mydata %>% mutate_at(vars(starts_with('dp_')), as.numeric)
mydata$dep_score <- -9+mydata$dp_tired + mydata$dp_slow + mydata$dp_sleep + mydata$dp_selfharm + mydata$dp_littleinterest + mydata$dp_feeldown + mydata$dp_failure + mydata$dp_concentrate + mydata$dp_appetite
mydata$dep_scoregt5 <- ifelse(mydata$dep_score>5, 'Yes', 'No')
mydata$dep_scoregt10 <- ifelse(mydata$dep_score>10, 'Yes', 'No')


# social support score
table(mydata$ss_supportdoc)
mydata$ss_score <- mydata$ss_supportbed + mydata$ss_supportlisten + mydata$ss_supportcrisis + 
  mydata$ss_supportdoc + mydata$ss_supportlove + mydata$ss_supportfun + mydata$ss_supportinfo + 
  mydata$ss_supportprob + mydata$ss_supportrelax + mydata$ss_supportmeal + mydata$ss_supportadvice +
  mydata$ss_supportmindoff + mydata$ss_supportchores + mydata$ss_supportworries + mydata$ss_supportsuggest + 
  mydata$ss_supportenjoy + mydata$ss_supportunderstand + mydata$ss_supportwanted


# food security score
mydata <- mydata %>% mutate_at(vars(starts_with('fa_')), as.numeric)

mydata$faworry1 <- ifelse(mydata$fa_worryfood==3, mydata$fa_worryfoodfreq, 0)
mydata$faprefer2 <- ifelse(mydata$fa_foodnotprefer==3, mydata$fa_foodnotpreferfreq, 0)
mydata$falimited3 <- ifelse(mydata$fa_foodlimited==3, mydata$fa_foodlimitedfreq, 0)
mydata$falike4 <- ifelse(mydata$fa_foodnotlike==3, mydata$fa_foodnotlikefreq, 0)
mydata$fasmall5 <- ifelse(mydata$fa_smallmeal==3, mydata$fa_smallmealrfreq, 0)
mydata$fafew6 <- ifelse(mydata$fa_fewmeal==3, mydata$fa_fewmealfreq, 0)
mydata$falack7 <- ifelse(mydata$fa_nofood==3, mydata$fa_nofoodfreq, 0)
mydata$fasleep8 <- ifelse(mydata$fa_sleephungry==3, mydata$fa_sleephungryfreq, 0)
mydata$faday9 <- ifelse(mydata$fa_wholedaynofood==3, mydata$fa_wholedaynofoodfreq, 0)

mydata$hfia_cat <- ifelse(mydata$faworry1<=1 & mydata$faprefer2==0 & mydata$falimited3==0 & mydata$falike4==0 & mydata$fasmall5==0 & mydata$fafew6==0 & mydata$falack7==0 & mydata$fasleep8==0 & mydata$faday9==0, 'level_1', 
                          ifelse((mydata$faworry1==2 | mydata$faworry1==3 | mydata$faworry1==3 | mydata$faprefer2==1 | mydata$faprefer2==2 | mydata$faprefer2==3 | mydata$falimited3==1 | mydata$falike4==1) & mydata$fasmall5==0 & mydata$fafew6==0 & mydata$falack7==0 & mydata$fasleep8==0 & mydata$faday9==0, 'level_2', 
                                 ifelse((mydata$falimited3==2 | mydata$falimited3==3 | mydata$falike4==2 | mydata$falike4==3 | mydata$fasmall5==1 | mydata$fasmall5==2 | mydata$fafew6==1 | mydata$fafew6==2) & mydata$falack7==0 & mydata$fasleep8==0 & mydata$faday9==0, 'level_3', 
                                        ifelse(mydata$fasmall5==3 | mydata$fafew6==3 | mydata$falack7==1 | mydata$falack7==2 | mydata$falack7==3 | mydata$fasleep8==1 | mydata$fasleep8==2 | mydata$fasleep8==3 | mydata$faday9==1 | mydata$faday9==2 | mydata$faday9==3, 'level_4', 0))))

mydata$hfia_bin <- ifelse(mydata$hfia_cat=='level_3' | mydata$hfia_cat=='level_4', 'insecure', 'foodsecure')
table(mydata$hfia_bin, mydata$hfia_cat)

# 1 "Food secure" 2 "Mildly food insecure" 3 "Moderately food insecure" 4 "Severely food insecure"


# travel time to clinic
mydata$timetoclinic <- ifelse(mydata$en_scale_clinictravel=='Hours (round to nearest)', 60*mydata$en_clinictravelm, mydata$en_clinictravelm)
mydata$timetoclinic <- ifelse(mydata$timetoclinic<0, NA, mydata$timetoclinic)
  

#######################################################
# enrollment VL and art date
########################################################
vl <- read.csv('7. VL results/0. allVL_merged/VL_23June2020_long.csv') %>% subset(from=='Aenrollment')
art <- read.csv('3. DataMerges/ARTStartDate/artdate.csv')
art$artstartdate <- as.Date(art$artstartdate, '%m/%d/%Y')

# art time
mydata <- merge(mydata, art, all.x = T)
mydata$enroldate <- as.Date(mydata$enroldate, '%m/%d/%Y')
mydata$en_lmp <- as.Date(mydata$en_lmp)
mydata$artb4preg <- ifelse(mydata$artstartdate - mydata$en_lmp <0, 'Yes', 'No')
table(mydata$artb4preg)


# merge mydata with enroll VL
mydata <- merge(mydata, vl, by='ptidno', all.x = T)
mydata$enrolvl <- as.numeric(as.character(mydata$Result))


#######################################################
# read infant sex data
#######################################################
infsex <- read.csv('16. Infant sex/infant sex_20210602.csv') %>% dplyr::select(-X)
table(infsex$sex)

# merge mydata with infant sex
mydata <- merge(mydata, infsex, by='ptidno', all.x = T)

# change sae to NA
mydata$sex <- ifelse(mydata$sex=='sae', NA, mydata$sex)

addmargins(table(mydata$sex, mydata$nonalive, useNA = 'ifany'))


# check data availibility
View(as.data.frame(summary(tableby(sex ~ as.character(adv_any)+
                                     as.character(adv_sb)+
                                     as.character(adv_esb)+
                                     as.character(adv_ptb)+
                                     as.character(adv_vptb)+
                                     as.character(adv_nnd),
                                   data=mydata, control=desc))))

tab <- (tableby(sex ~ as.character(adv_any)+
                  as.character(adv_sb)+
                  as.character(adv_esb)+
                  as.character(adv_ptb)+
                  as.character(adv_vptb)+
                  as.character(adv_nnd),
                data=mydata, control=desc))

dat <- as.data.frame(summary(tab))

###################################################################
# Table 1. Descriptive demographics, total N=774
###################################################################

## use the formula in the tableby function to output all statistics
desc  <- tableby.control(test=TRUE, total=TRUE,
                         numeric.stats=c("N", "medianq1q3"), digits = 0,
                         cat.stats = c("N", "countpct"), digits.count = 0, digits.pct = 1,
                         stats.labels = list(N='Count', medianq1q3 = "Median (Q1, Q3)", countpct = "Count (Pct)"))

# change outcome here
tab <- (tableby( ~ site+en_age+aya+(en_age<20)+
                  primarycom+secondarycom+en_inschool+
                  curmarried+hfia_cat+hfia_bin+(hfia_cat=='level_4')+
                  employ+income+(income>10000)+
                  dep_scoregt5+dep_scoregt10+
                  ss_score+(ss_score<64)+
                  abuselastyr+abuserecent+
                  en_ancvisits+ancvisitgt1+timetoclinic+(timetoclinic>30)+(timetoclinic>60)+
                  gestageatenr+(gestageatenr<196)+
                  fpever+primi+sti+syph+fetaldeathhis+infdeathhis+pregwant+
                  know+disc+(enrolvl>1000)+cd4count+(cd4count<350)+currentarv+
                  adreg_1st_bi+adreg_3rd_bi+artb4preg+imbscore_pct+parttested,
                data=mydata, control=desc))

dat <- as.data.frame(summary(tab))


###################################################################
# Table 2. Descriptive of outcomes
###################################################################
tab <- (tableby(site ~ as.character(adv_any)+
                   as.character(adv_sb)+
                   as.character(adv_esb)+
                   as.character(adv_ptb)+
                   as.character(adv_vptb)+
                   as.character(adv_nnd),
                 data=mydata, control=desc))

tab <- (tableby( ~ as.character(adv_any)+
                  as.character(adv_sb)+
                  as.character(adv_esb)+
                  as.character(adv_lsb)+
                  as.character(adv_tsb)+
                  as.character(adv_ptb)+
                  as.character(adv_vptb)+
                  as.character(adv_nnd),
                data=subset(mydata, mydata$gestageatenr<196), control=desc))

dat <- as.data.frame(summary(tab))


###################################################################
# Figure 2. Prevalence of adv by ART use (stat)
###################################################################
tab <- (tableby(adreg_1st_bi ~ adv_any+adv_sb+adv_ptb+adv_nnd+adv_esb+adv_vptb,
                data=mydata, control=desc))
dat <- as.data.frame(summary(tab))


tab <- (tableby(adreg_3rd ~ adv_any+adv_sb+adv_ptb+adv_nnd+adv_esb+adv_vptb,
                data=mydata, control=desc))
dat <- as.data.frame(summary(tab))

tab <- (tableby(adreg_3rd_bi ~ adv_any+adv_sb+adv_ptb+adv_nnd+adv_esb+adv_vptb,
                data=mydata, control=desc))
dat <- as.data.frame(summary(tab))


mydata$adreg_3drug <- paste(mydata$adreg_1st_bi,mydata$adreg_3rd_bi)
mydata$adreg_3drug <- ifelse(str_detect(mydata$adreg_3drug, 'NA')==T, 'unknown', mydata$adreg_3drug)

tab <- (tableby(adreg_3drug ~ adv_any+adv_sb+adv_ptb+adv_nnd+adv_esb+adv_vptb,
                data=mydata, control=desc))
dat <- as.data.frame(summary(tab))


tab <- (tableby(artb4preg ~ adv_any+adv_sb+adv_ptb+adv_nnd+adv_esb+adv_vptb,
                data=mydata, control=desc))
dat <- as.data.frame(summary(tab))


mydata$adreg_1st_bi_time <- paste(mydata$adreg_1st_bi,mydata$artb4preg)
mydata$adreg_1st_bi_time <- ifelse(str_detect(mydata$adreg_1st_bi_time, 'NA')==T, NA, mydata$adreg_1st_bi_time)

tab <- (tableby(adreg_1st_bi_time ~ adv_any+adv_sb+adv_ptb+adv_nnd+adv_esb+adv_vptb,
                data=mydata, control=desc))
dat <- as.data.frame(summary(tab))


mydata$adreg_3rd_bi_time <- paste(mydata$adreg_3rd_bi,mydata$artb4preg)
mydata$adreg_3rd_bi_time <- ifelse(str_detect(mydata$adreg_3rd_bi_time, 'NA')==T, NA, mydata$adreg_3rd_bi_time)

tab <- (tableby(adreg_3rd_bi_time ~ adv_any+adv_sb+adv_ptb+adv_nnd+adv_esb+adv_vptb,
                data=mydata, control=desc))

dat <- as.data.frame(summary(tab))

###################################################################
# Figure 3. time to stillbirth among women enrolled <20wk
# by  dep_scoregt5, VLenroll>1000, time-to-clinic>60min, parttested
###################################################################
fit <- survfit(Surv(gestageatenr/7, gestageatdel/7, adv_sb) ~ dep_scoregt5, data = subset(mydata, mydata$gestageatenr<140)); fit
ir <- pyears(Surv(gestageatenr/365, gestageatdel/365, adv_sb) ~ dep_scoregt5, data = subset(mydata, mydata$gestageatenr<140), scale = 1); summary(ir)
survdiff(Surv((gestageatdel-gestageatenr)/365, adv_sb) ~ dep_scoregt5, data = subset(mydata, mydata$gestageatenr<140))

# beautiful palette 
# palette = c("#E7B800", "#2E9FDF")
# palette = c("#FF9E29", "#86AA00")
# palette = "Dark2"

dp <- ggsurvplot(survfit(Surv(gestageatenr/7, gestageatdel/7, adv_sb) ~ dep_scoregt5, 
                         data=subset(mydata, mydata$gestageatenr<140)),
                   risk.table=T, conf.int=F, linetype = "strata",
           xlim=c(20,40), ylim=c(0.75,1), break.time.by=4, fontsize=4,
           legend.title = "Depression",
           legend.labs = c("No", "Yes"),
           xlab='Gestational age at delivery (week)', size=0.5, 
           ggtheme = theme_bw()+theme(panel.grid = element_blank(),
                                      plot.title = element_text(size = 10)))+guides(colour=guide_legend(nrow=2))


fit <- survfit(Surv(gestageatenr/7, gestageatdel/7, adv_sb) ~ parttested, data = subset(mydata, mydata$gestageatenr<140)); fit
ir <- pyears(Surv(gestageatenr/365, gestageatdel/365, adv_sb) ~ parttested, data = subset(mydata, mydata$gestageatenr<140), scale = 1); summary(ir)
survdiff(Surv((gestageatdel-gestageatenr)/365, adv_sb) ~ parttested, data = subset(mydata, mydata$gestageatenr<140))

pt <- ggsurvplot(survfit(Surv(gestageatenr/7, gestageatdel/7, adv_sb) ~ parttested, 
                         data = subset(mydata, mydata$gestageatenr<140)), 
                 risk.table=TRUE, conf.int=FALSE,  linetype = "strata",
           xlim=c(20,40), ylim=c(0.75,1), break.time.by=4, fontsize=4, 
           legend.title = "Partner HIV testing",
           legend.labs=c('Not tested', 'Tested'),        
           xlab='Gestational age at delivery (week)', size=0.5, 
           palette = c("#E7B800", "#2E9FDF"),
           ggtheme = theme_bw()+theme(panel.grid = element_blank(),
                                      plot.title = element_text(size = 10)))+guides(colour=guide_legend(nrow=2))

ss <- ggsurvplot(survfit(Surv(gestageatenr/7, gestageatdel/7, adv_lsb) ~ (ss_score<63), 
                         data = subset(mydata, mydata$gestageatenr<196)), 
                 risk.table=TRUE, conf.int=FALSE,  linetype = "strata",
                 xlim=c(20,40), ylim=c(0.75,1), break.time.by=4, fontsize=4, 
                 legend.title = "Social support score",
                 legend.labs=c('≥median', '<median'),        
                 xlab='Gestational age at delivery (week)', size=0.5, 
                 palette = c("#FF9E29", "#86AA00"),
                 ggtheme = theme_bw()+theme(panel.grid = element_blank(),
                                            plot.title = element_text(size = 10)))+guides(colour=guide_legend(nrow=2))

grid.arrange(dp$plot, pt$plot, ss$plot, ncol=1, nrow=3)

###################################################################
# Figure 4. time to stillbirth 
# by ART regimen and timing
###################################################################
fit <- survfit(Surv(gestageatenr/7, gestageatdel/7, adv_sb) ~ artb4preg, 
               data = subset(mydata, mydata$gestageatenr<140)); fit

fit <- survfit(Surv(gestageatenr/7, gestageatdel/7, adv_lsb) ~ artb4preg, 
               data = subset(mydata, mydata$gestageatenr<196)); fit

pyears(Surv((gestageatdel-gestageatenr)/365, adv_sb) ~ 1, data = subset(mydata, mydata$gestageatenr<140), scale = 1)
pyears(Surv((gestageatdel-gestageatenr)/365, adv_lsb) ~ 1, data = subset(mydata, mydata$gestageatenr<196), scale = 1)

summary(pyears(Surv(gestageatenr, gestageatdel, adv_sb) ~ artb4preg, data = subset(mydata, mydata$gestageatenr<140)))


survdiff(Surv((gestageatdel-gestageatenr)/365, adv_sb) ~ adreg_1st_bi, data = subset(mydata, mydata$gestageatenr<140))
survdiff(Surv((gestageatdel-gestageatenr)/365, adv_sb) ~ adreg_3rd_bi, data = subset(mydata, mydata$gestageatenr<140))
survdiff(Surv((gestageatdel-gestageatenr)/365, adv_sb) ~ artb4preg, data = subset(mydata, mydata$gestageatenr<140))


ggsurvplot(fit, data=subset(mydata, mydata$gestageatenr<140), risk.table=TRUE, conf.int=FALSE, linetype = "strata",
           xlim=c(20,48), ylim=c(0.7,1), break.time.by=4, fontsize=4, 
           xlab='Gestational age at delivery (week)', size=0.5, 
           legend.labs=c('On ART during pregnancy', 'On ART pre-conception'),
           ggtheme = theme_bw()+theme(panel.grid = element_blank()))


###################################################################
# Table 2. any stillbirth among women enrolled <20 wk
###################################################################

############## prevalence ##############
tab <- (tableby(adv_sb ~ site+en_age+aya+(en_age<20)+
                   primarycom+secondarycom+en_inschool+
                   curmarried+employ+(income>10000)+
                  dep_scoregt5+dep_scoregt10+
                  hfia_bin+(hfia_cat=='level_4')+
                   (ss_score<64)+
                   abuselastyr+abuserecent+
                   ancvisitgt1+(timetoclinic>60)+
                   gestageatenr+(gestageatenr<196)+
                   primi+sti+syph+pregwant+
                   know+disc+(enrolvl>1000)+(cd4count<400)+(cd4count<350)+(cd4count<200)+
                   adreg_1st_bi+adreg_3rd_bi+artb4preg+(imbscore_pct<75)+parttested+
                   sex,
                 data=subset(mydata,mydata$gestageatenr<140), control=desc))
dat <- as.data.frame(summary(tab))


############## PR ##############
demo <- list('site', 'en_age', "aya", "(en_age<20)", 
             'primarycom', 'secondarycom', 'en_inschool', "curmarried", 
             "employ","(income>10000)",
             "dep_scoregt5","dep_scoregt10","(hfia_cat=='level_3' | hfia_cat=='level_4')","(hfia_cat=='level_4')",
             '(ss_score<64)','abuselastyr','abuserecent',
             "ancvisitgt1",'(timetoclinic>60)', 
             'primi', 'sti', 'syph', "pregwant",
             "know","disc",
             '(enrolvl>1000)', '(cd4count<400)','(cd4count<350)','(cd4count<200)',
             "adreg_1st_bi","adreg_3rd_bi","artb4preg",
             '(imbscore_pct<75)','parttested', 'sex')

result <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(result) <- c('PR (95%CI)', 'p-value')

for (i in demo){
  reg <- glm(paste("adv_sb ~ gestageatenr+site+", i[[1]]), 
             data=subset(mydata, mydata$gestageatenr<140), family="binomial"(link='log'))
  rob <- coeftest(reg, vcov. = vcovHAC)
  PR <- round(exp(rob[,1]), 2)
  lower <- round(exp(rob[,1]+qnorm(0.025)*rob[,2]), 2)
  upper <- round(exp(rob[,1]+qnorm(0.975)*rob[,2]), 2)
  CI <- paste0('(', lower, '-', upper, ')')
  p <- round(rob[,4], 3)
  PRCIp <- paste0(PR, ' ', CI,'; ',p)
  row <- as.data.frame(cbind(i, PRCIp))[4, ]
  result <- rbind(result, row)
}

reg <- glm(adv_sb ~ site+gestageatenr+adreg_3rd_bi, data=subset(mydata, mydata$gestageatenr<140), family="binomial"(link='log'))
summary(reg)
rob <- coeftest(reg, vcov. = vcovHAC)
PR <- round(exp(rob[,1]), 2)
lower <- round(exp(rob[,1]+qnorm(0.025)*rob[,2]), 2)
upper <- round(exp(rob[,1]+qnorm(0.975)*rob[,2]), 2)
CI <- paste0('(', lower, '-', upper, ')')
p <- round(rob[,4], 3)
print(cbind(PR, CI, p))



############## incidence rate ##############
demo <- list('(site=="Western")', "(en_age<25)", "(en_age<20)", 
             'primarycom', 'secondarycom', 'en_inschool', "curmarried", 
             "employ","(income>10000)",
             "dep_scoregt5","dep_scoregt10","(hfia_cat=='level_3' | hfia_cat=='level_4')","(hfia_cat=='level_4')",
             '(ss_score<64)','abuselastyr','abuserecent',
             "ancvisitgt1",'(timetoclinic>60)', 
             'primi', 'sti', 'syph', "pregwant",
             "know","disc",'(enrolvl>1000)', '(cd4count<400)', '(cd4count<350)', '(cd4count<200)', 
             "(adreg_1st_bi=='TDF')","(adreg_3rd_bi=='NVP')","artb4preg",
             '(imbscore_pct<75)','parttested', '(sex=="male")')


result_ir <- data.frame(matrix(nrow = 0, ncol = 4))
colnames(result_ir) <- c('var', 'N', 'eventpy', 'ir100py')
output_ir <- function(x) as.formula(paste("Surv((gestageatdel-gestageatenr), adv_sb) ~", x))

for (i in demo){
  eq <- pyears(output_ir(i), data = subset(mydata, mydata$gestageatenr<140))
  N <- eq$observations
  ir <- round(eq$event/eq$pyears*100, 2)
  eventpy <- paste0(eq$event, '/', round(eq$pyears, 1))
  row1 <- as.data.frame(cbind(i, N, eventpy, ir))[1,]
  row2 <- as.data.frame(cbind(i, N, eventpy, ir))[2,]
  result_ir <- rbind(result_ir, row1, row2)
}

result_ir$group <- rownames(result_ir)

############## cox PH model ##############
output_ph_unadjusted <- function(x) as.formula(paste("Surv(gestageatenr, gestageatdel, adv_sb) ~", x))
output_ph_adjusted <- function(x) as.formula(paste("Surv(gestageatenr, gestageatdel, adv_sb) ~site+", x))

result_ph <- data.frame(matrix(nrow = 0, ncol = 3))
colnames(result_ph) <- c('var', 'crude', 'adjusted')

for (i in demo){
  # index 1 for site-unadjusted 
  reg1 <- coxph(output_ph_unadjusted(i), 
               data = subset(mydata, mydata$gestageatenr<140), robust=T, id=ptidno)
  beta1 <- coef(summary(reg1))[,1]
  HR1 <- round(exp(beta1), 2)
  err1 <- sqrt(diag(reg1$var))
  lower1 <- round(exp(beta1 - 1.96*err1), 2)
  upper1 <- round(exp(beta1 + 1.96*err1), 2)
  CI1 <- paste0('(', lower1, '-', upper1, ')')
  p1 <- round(coef(summary(reg1))[,6], 3)
  HRCIp1 <- paste0(HR1, ' ', CI1,'; ', p1)

  # index 2 for site-adjusted 
  reg2 <- coxph(output_ph_adjusted(i), 
                data = subset(mydata, mydata$gestageatenr<140), robust=T, id=ptidno)
  beta2 <- coef(summary(reg2))[-1,1]
  HR2 <- round(exp(beta2), 2)
  err2 <- sqrt(diag(reg2$var))
  lower2 <- round(exp(beta2 - 1.96*err2), 2)
  upper2 <- round(exp(beta2 + 1.96*err2), 2)
  CI2 <- paste0('(', lower2, '-', upper2, ')')
  p2 <- round(coef(summary(reg2))[,6], 3)
  HRCIp2 <- paste0(HR2, ' ', CI2,'; ', p2)
  
  # cbind crude HR and adjusted HR
  row <- as.data.frame(cbind(i, HRCIp1, HRCIp2))[2,]
  result_ph <- rbind(result_ph, row)
}
# add index to keep the order of variables
result_ph$x <- as.numeric(rownames(result_ph))


# merge IR data and HR data by group
result_merge_by <- merge(result_ir, result_ph, by='i', all.y=T) %>% 
  # change HR of ref group to "Ref"
  mutate(HRCIp1 = ifelse(str_sub(group, 1, 3)=='Yes' | str_sub(group, 1, 4)=='TRUE', HRCIp1, 'Ref'),
         HRCIp2 = ifelse(str_sub(group, 1, 3)=='Yes' | str_sub(group, 1, 4)=='TRUE', HRCIp2, 'Ref')) 
# generate a total group showing N
result_merge_all <- result_merge_by %>% arrange(x) %>% group_by(i) %>%
  summarise(N = last(N),
            x = last(x)) %>% ungroup() %>% mutate(group = i, HRCIp1='', HRCIp2='') 
# bind with by-group data
result_merge <- bind_rows(result_merge_by, result_merge_all) %>% 
  # arrange by variable order - x
  arrange(x, HRCIp1) %>%
  mutate(N=ifelse(!is.na(ir), '', N), # remove N in by-group data
    group=ifelse(!is.na(ir), paste0('  ', group), group)) %>% # add space to by-group name
  dplyr::select(group, N, eventpy, ir, HRCIp1, HRCIp2)
result_merge[is.na(result_merge)] <- ""

# recalculate median ss_score and imbscore_pct
summary(subset(mydata, gestageatenr<196)$ss_score)
summary(subset(mydata, gestageatenr<140)$imbscore_pct)

reg <- coxph(Surv(gestageatenr, gestageatdel, adv_sb) ~site+sti, 
      data = subset(mydata, mydata$gestageatenr<140), robust=T, id=ptidno)

reg <- coxph(Surv(gestageatenr, gestageatdel, adv_lsb) ~site+sti, 
             data = subset(mydata, mydata$gestageatenr<196), robust=T, id=ptidno)

beta <- coef(summary(reg))[,1]
HR <- round(exp(beta), 2)
err <- sqrt(diag(reg$var))
lower <- round(exp(beta - 1.96*err), 2)
upper <- round(exp(beta + 1.96*err), 2)
CI <- paste0('(', lower, '-', upper, ')')
p <- round(coef(summary(reg))[,6], 3)
print(cbind(HR, CI, p))


###################################################################
# Table 3. late stillbirth among women enrolled <28 wk
###################################################################
fit <- survfit(Surv(gestageatenr/7, gestageatdel/7, adv_lsb) ~ hfia_bin, data = subset(mydata, mydata$gestageatenr<196)); fit
survdiff(Surv((gestageatdel-gestageatenr)/365, adv_lsb) ~ (ip_sexabuse=='Yes'), data = subset(mydata, mydata$gestageatenr<196))

ss <- ggsurvplot(survfit(Surv(gestageatenr/7, gestageatdel/7, adv_lsb) ~ (ss_score>63), 
                         data = subset(mydata, mydata$gestageatenr<196)), 
                 risk.table=TRUE, conf.int=FALSE,  linetype = "strata",
                 xlim=c(28,48), ylim=c(0.9,1), break.time.by=4, fontsize=4, 
                 legend.title = "Social support score",
                 legend.labs=c('Below 63', 'Above 63'),        
                 xlab='Gestational age at delivery (week)', size=0.5, 
                 ggtheme = theme_bw()+theme(panel.grid = element_blank(),
                                            plot.title = element_text(size = 10)))

pt <- ggsurvplot(survfit(Surv(gestageatenr/7, gestageatdel/7, adv_lsb) ~ parttested, 
                         data = subset(mydata, mydata$gestageatenr<196)), 
                 risk.table=TRUE, conf.int=FALSE,  linetype = "strata",
                 xlim=c(28,48), ylim=c(0.9,1), break.time.by=4, fontsize=4, 
                 legend.title = "Partner tested for HIV",
                 legend.labs=c('No', 'Yes'),        
                 xlab='Gestational age at delivery (week)', size=0.5, 
                 ggtheme = theme_bw()+theme(panel.grid = element_blank(),
                                            plot.title = element_text(size = 10)))

grid.arrange(ss$plot, pt$plot, ss$table, pt$table, ncol=2, nrow=2, heights=c(2.8,1))

############## prevalence ##############
tab <- (tableby(adv_lsb ~ site+en_age+aya+(en_age<20)+
                  primarycom+secondarycom+en_inschool+
                  curmarried+employ+(income>10000)+
                  dep_scoregt5+dep_scoregt10+
                  hfia_bin+(hfia_cat=='level_4')+
                  (ss_score<63)+
                  abuselastyr+abuserecent+
                  ancvisitgt1+(timetoclinic>60)+
                  gestageatenr+(gestageatenr<196)+
                  primi+sti+syph+pregwant+
                  know+disc+(enrolvl>1000)+(cd4count<400)+(cd4count<350)+(cd4count<200)+
                  adreg_1st_bi+adreg_3rd_bi+artb4preg+(imbscore_pct<75)+parttested+
                  sex,
                data=subset(mydata,mydata$gestageatenr<196), control=desc))
dat <- as.data.frame(summary(tab))


result <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(result) <- c('PR (95%CI)', 'p-value')

for (i in demo){
  reg <- glm(paste("adv_lsb ~ gestageatenr+site+", i[[1]]), 
             data=subset(mydata, mydata$gestageatenr<196), family="binomial"(link='log'))
  rob <- coeftest(reg, vcov. = vcovHAC)
  PR <- round(exp(rob[,1]), 2)
  lower <- round(exp(rob[,1]+qnorm(0.025)*rob[,2]), 2)
  upper <- round(exp(rob[,1]+qnorm(0.975)*rob[,2]), 2)
  CI <- paste0('(', lower, '-', upper, ')')
  p <- round(rob[,4], 3)
  PRCIp <- paste0(PR, ' ', CI,'; ',p)
  row <- as.data.frame(cbind(i, PRCIp))[4, ]
  result <- rbind(result, row)
}
############## HR ##############
demo <- list('(site=="Western")', "(en_age<25)", "(en_age<20)", 
             'primarycom', 'secondarycom', 'en_inschool', "curmarried", 
             "employ","(income>10000)",
             "dep_scoregt5","dep_scoregt10","(hfia_cat=='level_3' | hfia_cat=='level_4')","(hfia_cat=='level_4')",
             '(ss_score<63)','abuselastyr','abuserecent',
             "ancvisitgt1",'(timetoclinic>60)', 
             'primi', 'sti', 'syph', "pregwant",
             "know","disc",'(enrolvl>1000)', '(cd4count<400)', '(cd4count<350)', '(cd4count<200)', 
             "(adreg_1st_bi=='TDF')","(adreg_3rd_bi=='NVP')","artb4preg",
             '(imbscore_pct<75)','parttested', '(sex=="male")')

result_ir <- data.frame(matrix(nrow = 0, ncol = 4))
colnames(result_ir) <- c('var', 'N', 'eventpy', 'ir100py')
output_ir <- function(x) as.formula(paste("Surv((gestageatdel-gestageatenr), adv_lsb) ~", x))

for (i in demo){
  eq <- pyears(output_ir(i), data = subset(mydata, mydata$gestageatenr<196))
  N <- eq$observations
  ir <- round(eq$event/eq$pyears*100, 2)
  eventpy <- paste0(eq$event, '/', round(eq$pyears, 1))
  row1 <- as.data.frame(cbind(i, N, eventpy, ir))[1,]
  row2 <- as.data.frame(cbind(i, N, eventpy, ir))[2,]
  result_ir <- rbind(result_ir, row1, row2)
}

result_ir$group <- rownames(result_ir)


############## cox PH model ##############
output_ph_unadjusted <- function(x) as.formula(paste("Surv(gestageatenr, gestageatdel, adv_lsb) ~", x))
output_ph_adjusted <- function(x) as.formula(paste("Surv(gestageatenr, gestageatdel, adv_lsb) ~site+", x))

result_ph <- data.frame(matrix(nrow = 0, ncol = 3))
colnames(result_ph) <- c('var', 'crude', 'adjusted')

for (i in demo){
  # index 1 for site-unadjusted 
  reg1 <- coxph(output_ph_unadjusted(i), 
                data = subset(mydata, mydata$gestageatenr<196), robust=T, id=ptidno)
  beta1 <- coef(summary(reg1))[,1]
  HR1 <- round(exp(beta1), 2)
  err1 <- sqrt(diag(reg1$var))
  lower1 <- round(exp(beta1 - 1.96*err1), 2)
  upper1 <- round(exp(beta1 + 1.96*err1), 2)
  CI1 <- paste0('(', lower1, '-', upper1, ')')
  p1 <- round(coef(summary(reg1))[,6], 3)
  HRCIp1 <- paste0(HR1, ' ', CI1,'; ', p1)
  
  # index 2 for site-adjusted 
  reg2 <- coxph(output_ph_adjusted(i), 
                data = subset(mydata, mydata$gestageatenr<196), robust=T, id=ptidno)
  beta2 <- coef(summary(reg2))[-1,1]
  HR2 <- round(exp(beta2), 2)
  err2 <- sqrt(diag(reg2$var))
  lower2 <- round(exp(beta2 - 1.96*err2), 2)
  upper2 <- round(exp(beta2 + 1.96*err2), 2)
  CI2 <- paste0('(', lower2, '-', upper2, ')')
  p2 <- round(coef(summary(reg2))[,6], 3)
  HRCIp2 <- paste0(HR2, ' ', CI2,'; ', p2)
  
  # cbind crude HR and adjusted HR
  row <- as.data.frame(cbind(i, HRCIp1, HRCIp2))[2,]
  result_ph <- rbind(result_ph, row)
}
# add index to keep the order of variables
result_ph$x <- as.numeric(rownames(result_ph))


# merge IR data and HR data by group
result_merge_by <- merge(result_ir, result_ph, by='i', all.y=T) %>% 
  # change HR of ref group to "Ref"
  mutate(HRCIp1 = ifelse(str_sub(group, 1, 3)=='Yes' | str_sub(group, 1, 4)=='TRUE', HRCIp1, 'Ref'),
         HRCIp2 = ifelse(str_sub(group, 1, 3)=='Yes' | str_sub(group, 1, 4)=='TRUE', HRCIp2, 'Ref')) 
# generate a total group showing N
result_merge_all <- result_merge_by %>% arrange(x) %>% group_by(i) %>%
  summarise(N = last(N),
            x = last(x)) %>% ungroup() %>% mutate(group = i, HRCIp1='', HRCIp2='') 
# bind with by-group data
result_merge <- bind_rows(result_merge_by, result_merge_all) %>% 
  # arrange by variable order - x
  arrange(x, HRCIp1) %>%
  mutate(N=ifelse(!is.na(ir), '', N), # remove N in by-group data
         group=ifelse(!is.na(ir), paste0('  ', group), group)) %>% # add space to by-group name
  dplyr::select(group, N, eventpy, ir, HRCIp1, HRCIp2)
result_merge[is.na(result_merge)] <- ""


###################################################################
# Table 4. PR of PTB among women enrolled <37 week
###################################################################
table(mydata$adv_ptb)
addmargins(table(mydata$adv_ptb, (mydata$gestageatenr<259)))

# recalculate median ss_score and imbscore_pct
summary(subset(mydata,gestageatenr<259)$ss_score)
summary(subset(mydata,gestageatenr<259)$imbscore_pct)

tab <- (tableby(adv_ptb ~ site+en_age+aya+(en_age<20)+
                  primarycom+secondarycom+en_inschool+
                  curmarried+employ+(income>10000)+
                  dep_scoregt5+dep_scoregt10+
                  hfia_bin+(hfia_cat=='level_4')+
                  (ss_score<64)+
                  abuselastyr+abuserecent+
                  ancvisitgt1+(timetoclinic>60)+
                  gestageatenr+(gestageatenr<196)+
                  primi+sti+syph+pregwant+
                  know+disc+(enrolvl>1000)+(cd4count<400)+(cd4count<350)+(cd4count<200)+
                  adreg_1st_bi+adreg_3rd_bi+artb4preg+(imbscore_pct<75)+parttested+
                  sex,
                data=subset(mydata,mydata$gestageatenr<259), control=desc))

dat <- as.data.frame(summary(tab))

# log-binomial model
demo <- list('(site=="Western")', "(en_age<25)", "(en_age<20)", 
             'primarycom', 'secondarycom', 'en_inschool', "curmarried", 
             "employ","(income>10000)",
             "dep_scoregt5","dep_scoregt10","(hfia_cat=='level_3' | hfia_cat=='level_4')","(hfia_cat=='level_4')",
             '(ss_score<64)','abuselastyr','abuserecent',
             "ancvisitgt1",'(timetoclinic>60)', 
             'primi','sti', 'syph', "pregwant",
             "know","disc",'(enrolvl>1000)', '(cd4count<400)', '(cd4count<350)', '(cd4count<200)', 
             "(adreg_1st_bi=='TDF')","(adreg_3rd_bi=='NVP')","artb4preg",
             '(imbscore_pct<75)','parttested', '(sex=="male")')

result <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(result) <- c('PR (95% CI)', 'p-value')

for (i in demo){
  reg <- glm(paste("adv_ptb ~ site+gestageatenr+", i[[1]]), 
             data=subset(mydata, mydata$gestageatenr<259), family="binomial"(link='log'))
  rob <- coeftest(reg, vcov. = vcovHAC)
  beta <- rob[,1]
  err <- rob[,2]
  PR <- round(exp(beta), 2)
  lower <- round(exp(beta+qnorm(0.025)*err), 2)
  upper <- round(exp(beta+qnorm(0.975)*err), 2)
  CI <- paste0('(', lower, '-', upper, ')')
  p <- round(rob[,4], 3)
  PRCIp <- paste0(PR, ' ', CI,'; ',p)
  row <- as.data.frame(cbind(i, PRCIp))[4, ]
  result <- rbind(result, row)
}

reg <- glm(adv_ptb ~ site+gestageatenr+sti, 
           data=subset(mydata, mydata$gestageatenr<259), family="binomial"(link='log'))


reg <- glm(adv_ptb ~ site+secondarycom+abuserecent+syph+(enrolvl>1000), 
           data=subset(mydata, mydata$gestageatenr<259), family="binomial"(link='log'))

reg <- glm(adv_vptb ~ site+secondarycom+primi+(enrolvl>1000), 
           data=subset(mydata, mydata$gestageatenr<224), family="binomial"(link='log'))

reg <- glm(adv_nnd ~ site+syph+gestageatdel, 
           data=mydata, family="binomial"(link='log'))

reg <- glm(adv_nnd ~ site+syph+gestageatdel, 
           data=mydata, family="binomial"(link='log'))

summary(reg)
rob <- coeftest(reg, vcov. = vcovHAC)
PR <- round(exp(rob[,1]), 2)
lower <- round(exp(rob[,1]+qnorm(0.025)*rob[,2]), 2)
upper <- round(exp(rob[,1]+qnorm(0.975)*rob[,2]), 2)
CI <- paste0('(', lower, '-', upper, ')')
p <- round(rob[,4], 3)
print(cbind(PR, CI, p))


###################################################################
# Table 6. PR of very PTB among women enrolled <28 week
###################################################################
table(mydata$adv_vptb, (mydata$gestageatenr<224))

# recalculate median ss_score and imbscore_pct
summary(subset(mydata,gestageatenr<224)$ss_score)
summary(subset(mydata,gestageatenr<224)$imbscore_pct)


tab <- (tableby(adv_vptb ~ site+en_age+aya+(en_age<20)+
                  primarycom+secondarycom+en_inschool+
                  curmarried+employ+(income>10000)+
                  dep_scoregt5+dep_scoregt10+
                  hfia_bin+(hfia_cat=='level_4')+
                  (ss_score<64)+
                  abuselastyr+abuserecent+
                  ancvisitgt1+(timetoclinic>60)+
                  gestageatenr+(gestageatenr<196)+
                  primi+sti+syph+pregwant+
                  know+disc+(enrolvl>1000)+(cd4count<400)+(cd4count<350)+(cd4count<200)+
                  adreg_1st_bi+adreg_3rd_bi+artb4preg+(imbscore_pct<75)+parttested+
                  sex,
                data=subset(mydata,mydata$gestageatenr<224), control=desc))

dat <- as.data.frame(summary(tab))

# log-binomial model - adjust for site
result <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(result) <- c('PR (95% CI)', 'p-value')

for (i in demo){
  reg <- glm(paste("adv_vptb ~ site+", i[[1]]), 
             data=subset(mydata, mydata$gestageatenr<224), family="binomial"(link='log'))
  rob <- coeftest(reg, vcov. = vcovHAC)
  beta <- rob[,1]
  err <- rob[,2]
  PR <- round(exp(beta), 2)
  lower <- round(exp(beta+qnorm(0.025)*err), 2)
  upper <- round(exp(beta+qnorm(0.975)*err), 2)
  CI <- paste0('(', lower, '-', upper, ')')
  p <- round(rob[,4], 3)
  PRCIp <- paste0(PR, ' ', CI,'; ',p)
  row <- as.data.frame(cbind(i, PRCIp))[3, ]
  result <- rbind(result, row)
}



reg <- glm(adv_vptb ~ (adreg_3rd=='EFV'), data=subset(mydata, mydata$gestageatenr<224), family="binomial"(link='log'))
summary(reg)
rob <- coeftest(reg, vcov. = vcovHAC)
PR <- round(exp(rob[,1]), 2)
lower <- round(exp(rob[,1]+qnorm(0.025)*rob[,2]), 2)
upper <- round(exp(rob[,1]+qnorm(0.975)*rob[,2]), 2)
CI <- paste0('(', lower, '-', upper, ')')
p <- round(rob[,4], 3)
print(cbind(PR, CI, p))
print(paste0(PR, ' ', CI,'; ',p))

########################################################
# Neonatal death log-binomial model
########################################################
tab <- (tableby(adv_nnd ~ site+en_age+aya+(en_age<20)+
                  primarycom+secondarycom+en_inschool+
                  curmarried+employ+(income>10000)+
                  dep_scoregt5+dep_scoregt10+
                  hfia_bin+(hfia_cat=='level_4')+
                  (ss_score<64)+
                  abuselastyr+abuserecent+
                  ancvisitgt1+(timetoclinic>60)+
                  gestageatenr+(gestageatenr<196)+
                  primi+sti+syph+pregwant+
                  know+disc+(enrolvl>1000)+(cd4count<400)+(cd4count<350)+(cd4count<200)+
                  adreg_1st_bi+adreg_3rd_bi+artb4preg+(imbscore_pct<75)+parttested+
                  sex+adv_ptb,
                data=mydata, control=desc))

dat <- as.data.frame(summary(tab))


result <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(result) <- c('PR (95% CI)', 'p-value')

for (i in demo){
  reg <- glm(paste("adv_nnd ~", i[[1]]), 
             data=mydata, family="binomial"(link='log'))
  rob <- coeftest(reg, vcov. = vcovHAC)
  beta <- rob[,1]
  err <- rob[,2]
  PR <- round(exp(beta), 2)
  lower <- round(exp(beta+qnorm(0.025)*err), 2)
  upper <- round(exp(beta+qnorm(0.975)*err), 2)
  CI <- paste0('(', lower, '-', upper, ')')
  p <- round(rob[,4], 3)
  PRCIp <- paste0(PR, ' ', CI,'; ',p)
  row <- as.data.frame(cbind(i, PRCIp))[2, ]
  result <- rbind(result, row)
}

reg <- glm(adv_nnd ~ gestageatdel, data=mydata, family="binomial"(link='log'))
reg <- glm(adv_nnd ~ site+syph, data=mydata, family="binomial"(link='log'))

reg <- glm(adv_nnd ~ site+sti, data=subset(mydata, gestageatdel>=259), family="binomial"(link='log'))
reg <- glm(adv_nnd ~ site+syph, data=subset(mydata, gestageatdel>=259), family="binomial"(link='log'))


summary(reg)
rob <- coeftest(reg, vcov. = vcovHAC)
PR <- round(exp(rob[,1]), 2)
lower <- round(exp(rob[,1]+qnorm(0.025)*rob[,2]), 2)
upper <- round(exp(rob[,1]+qnorm(0.975)*rob[,2]), 2)
CI <- paste0('(', lower, '-', upper, ')')
p <- round(rob[,4], 3)
print(cbind(PR, CI, p))
print(paste0(PR, ' ', CI,'; ',p))



tab <- (tableby(adv_ptb ~ ip_abuserecentwho1+ip_sexabuse, data=subset(mydata, gestageatenr<249), control=desc))
as.data.frame(summary(tab))

addmargins(table(mydata$ip_abuserecentwho1, mydata$abuserecent))
chisq.test(mydata$dep_scoregt5, mydata$abuserecent)

# first: 1-84, second: 85-182
prop.table(table(mydata$ancvisitgt1, mydata$gestageatenr<182), 2)




tab <- (tableby(adv_ptb ~ (en_sponabortion>0)+(en_elecabortion>0)+
                  (en_stillbirths>0)+(en_childrendied>0),
                data=subset(mydata,mydata$gestageatenr<249), control=desc))

dat <- as.data.frame(summary(tab))

reg <- glm(adv_ptb ~ site+(en_childrendied>0), data=subset(mydata, gestageatenr<249), family="binomial"(link='log'))


summary(reg)
rob <- coeftest(reg, vcov. = vcovHAC)
PR <- round(exp(rob[,1]), 2)
lower <- round(exp(rob[,1]+qnorm(0.025)*rob[,2]), 2)
upper <- round(exp(rob[,1]+qnorm(0.975)*rob[,2]), 2)
CI <- paste0('(', lower, '-', upper, ')')
p <- round(rob[,4], 3)
print(cbind(PR, CI, p))






# check maternal weight data availibility
lmpdel <- mydata %>% dplyr::select(ptidno, en_lmp, deldate)

# women ever with wt data
hwmerge <- merge(hw, lmpdel, by='ptidno') 

hwmerge2 <- hwmerge %>% 
  # women ever with wt data at least 2 measurements between lmp and del
  subset(eventdate>=en_lmp & eventdate<deldate) %>%
  group_by(ptidno) %>% count() %>% subset(n>1)

hwmerge2 <- merge(hwmerge2, hwmerge, all.x=T) %>% subset(eventdate>=en_lmp & eventdate<deldate)

hwmerge2 %>% group_by(ptidno) %>% count() %>% subset(n>1)


hwmerge3 <- hwmerge2 %>% group_by(ptidno) %>% 
  summarise(wt_lmp_1mo=sum(eventdate<=en_lmp+30),
            wt_lmp_2mo=sum(eventdate<=en_lmp+60),
            wt_lmp_3mo=sum(eventdate<=en_lmp+90),
            wt_lmp_4mo=sum(eventdate<=en_lmp+120),
            wt_lmp_6mo=sum(eventdate<=en_lmp+180),
            
            wt_del_1wk=sum(eventdate>=deldate-7 & eventdate<=deldate),
            wt_del_2wk=sum(eventdate>=deldate-14 & eventdate<=deldate),
            wt_del_1mo=sum(eventdate>=deldate-30 & eventdate<=deldate),
            wt_del_2mo=sum(eventdate>=deldate-60 & eventdate<=deldate))

tab <- (tableby((wt_del_1wk>0) ~ (wt_lmp_1mo>0)+(wt_lmp_2mo>0)+(wt_lmp_3mo>0)+(wt_lmp_4mo>0)+(wt_lmp_6mo>0),
                data=wtmerge3, control=desc))
dat <- as.data.frame(summary(tab))


# calculate wt change between first & last
hwmerge3 <- hwmerge2 %>% arrange(ptidno, eventdate) %>% group_by(ptidno) %>% 
  mutate(interval = last(eventdate) - first(eventdate),
         wtchange = last(weight) - first(weight)) %>%
  # subset for positive change
  subset(interval>=30 & wtchange>0) %>%
  dplyr::select(ptidno, interval, wtchange) %>% unique()

# caculate monthly wt change
hwmerge3 <- hwmerge3 %>% mutate(monthlywt = wtchange / as.numeric(interval) * 30,
                                weeklywt = wtchange / as.numeric(interval) * 7)
  
mydata_new <- merge(mydata, hwmerge3, all.y=T)
summary(mydata_new$weeklywt)

# wt change with sb 0.42, 0.27
summary(pyears(Surv(gestageatenr, gestageatdel, adv_sb) ~ (weeklywt<0.27), data = subset(mydata_new, mydata_new$gestageatenr<140)))
survdiff(Surv((gestageatdel-gestageatenr)/365, adv_sb) ~ (weeklywt<0.27), data = subset(mydata_new, mydata_new$gestageatenr<140))

reg <- coxph(Surv(gestageatenr, gestageatdel, adv_sb) ~(weeklywt<0.27), 
             data = subset(mydata_new, mydata$gestageatenr<140), robust=T, id=ptidno)

summary(pyears(Surv(gestageatenr, gestageatdel, adv_lsb) ~ (weeklywt<0.27), data = subset(mydata_new, mydata_new$gestageatenr<196)))
survdiff(Surv((gestageatdel-gestageatenr)/365, adv_lsb) ~ (weeklywt<0.42), data = subset(mydata_new, mydata_new$gestageatenr<196))

reg <- coxph(Surv(gestageatenr, gestageatdel, adv_lsb) ~(weeklywt<0.27), 
             data = subset(mydata_new, mydata$gestageatenr<196), robust=T, id=ptidno)

beta <- coef(summary(reg))[,1]
HR <- round(exp(beta), 2)
err <- sqrt(diag(reg$var))
lower <- round(exp(beta - 1.96*err), 2)
upper <- round(exp(beta + 1.96*err), 2)
CI <- paste0('(', lower, '-', upper, ')')
p <- round(coef(summary(reg))[,6], 3)
print(cbind(HR, CI))

tab <- (tableby(adv_ptb ~ (weeklywt<0.27),
                data=subset(mydata_new,mydata_new$gestageatenr<259), control=desc))

tab <- (tableby(adv_vptb ~ (weeklywt<0.27),
                data=subset(mydata_new,mydata_new$gestageatenr<224), control=desc))

tab <- (tableby(adv_nnd ~ (weeklywt<0.27),
                data=mydata_new, control=desc))

dat <- as.data.frame(summary(tab))


reg <- glm(adv_vptb ~ site+(weeklywt<0.27), data=subset(mydata_new,mydata_new$gestageatenr<224), family="binomial"(link='log'))
reg <- glm(adv_nnd ~ site+(weeklywt<0.27), data=mydata_new, family="binomial"(link='log'))

summary(reg)
rob <- coeftest(reg, vcov. = vcovHAC)
PR <- round(exp(rob[,1]), 2)
lower <- round(exp(rob[,1]+qnorm(0.025)*rob[,2]), 2)
upper <- round(exp(rob[,1]+qnorm(0.975)*rob[,2]), 2)
CI <- paste0('(', lower, '-', upper, ')')
p <- round(rob[,4], 3)
print(cbind(PR, CI, p))

