# script for MI pooling studies
# report results at 24 & 48 months
# include only studies that can report at least 24 months data

# studies included - 
# LEADER, Harmony Outcomes, REWIND

# LEADER - liraglutide arm 

library(tidyverse)
library(metaRMST)
library(ggthemes)
library(survminer)
library(IPDfromKM)

li <- read.table("F:/GLP1_agonists/data/leader/mi_liraglutide.txt")

summary(li)

# convert V2 to survival - 

li$V2 <- 100 - li$V2

summary(li)

# now to extract information 

leader_l <- preprocess(dat = li, totalpts = 4668,
                       maxy = 100)


leader_ipd_li <- getIPD(prep = leader_l,
                        armID = 1)


plot(leader_ipd_li)

mi_li <- leader_ipd_li$IPD


# control arm 


control <- 
  read.table("F:/GLP1_agonists/data/leader/mi_control.txt")

summary(control)

# convert V2 to survival - 

control$V2 <- 100 - control$V2

summary(control)

# now to extract information 

leader_c <- preprocess(dat = control,totalpts = 4672,
maxy = 100)


leader_ipd_c <- getIPD(prep = leader_c,
                       armID = 0)


plot(leader_ipd_c)

mi_control <- leader_ipd_c$IPD


# combine data and then plot 

df_leader_mi <- rbind(mi_li, mi_control)

mi_leader_s <- 
  survfit(Surv(time, status) ~ treat, data = df_leader_mi)


ggsurvplot(mi_leader_s, 
           fun = "event",
           ylim = c(0,0.20),
           xlim = c(0,54),
           break.x.y = 6,
           break.y.by = 0.05,
           censor.size = 0,
           palette = c("gray","blue2"))

# time already in months.

write_csv(df_leader_mi,
          'F:\\GLP1_agonists\\data\\pooled_data\\mi_leader.csv'
)


# Harmony Outcomes 


abi <- read.table("F:/GLP1_agonists/data/Harmony/mi_abiglutide.txt")

summary(abi)

# convert V2 to survival - 

abi$V2 <- 100 - abi$V2

summary(abi)

# now to extract information 

har_a <- preprocess(dat = abi, totalpts = 4731,
                       maxy = 100)


harmony_ipd_abi <- getIPD(prep = har_a,
                        armID = 1)


plot(harmony_ipd_abi)

mi_abi <- harmony_ipd_abi$IPD


# Harmony outcomes - control arm 


control <- 
  read.table("F:/GLP1_agonists/data/Harmony/mi_control.txt")

summary(control)

# convert V2 to survival - 

control$V2 <- 100 - control$V2

summary(control)

# now to extract information 

harmony_c <- preprocess(dat = control,totalpts = 4732,
maxy = 100)


harmony_ipd_c <- getIPD(prep = harmony_c,
                       armID = 0)

plot(harmony_ipd_c)

harmony_mi_control <- harmony_ipd_c$IPD


# combine data and then plot 

df_harmony_mi <- rbind(mi_abi, harmony_mi_control)

# plot the graph now 

harmony_mi_s <- survfit(Surv(time, status) ~ treat, 
data = df_harmony_mi)

ggsurvplot(harmony_mi_s, 
           fun = "event",
           ylim = c(0,0.16),
           xlim = c(0,28),
           break.x.y = 4,
           break.y.by = 0.02,
           censor.size = 0,
           palette = c("lightblue", "red"))


# save the dataset
# time is in months, so no change 

write_csv(df_harmony_mi,
	'F:\\GLP1_agonists\\data\\pooled_data\\mi_harmony.csv')



# REWIND

# dulaglutide


dula <- 
  read.table("F:/GLP1_agonists/data/rewind/mi_dulaglutide.txt")

summary(dula)

# convert V2 to survival - 

dula$V2 <- 100 - dula$V2

summary(dula)

# now to extract information 

rewind_dula <- preprocess(dat = dula,
                          totalpts = 4949,
maxy = 100)


rewind_ipd_d <- getIPD(prep = rewind_dula,
                       armID = 1)

plot(rewind_ipd_d)

rewind_mi_dula <- rewind_ipd_d$IPD


# rewind, control arm 


control <- 
  read.table("F:/GLP1_agonists/data/rewind/mi_control.txt")

summary(control)

# convert V2 to survival - 

control$V2 <- 100 - control$V2

summary(control)

# now to extract information 

rewind_c <- preprocess(dat = control,
                       totalpts = 4952,
                        maxy = 100)


rewind_ipd_c <- getIPD(prep = rewind_c,
                        armID = 0)

plot(rewind_ipd_c)

rewind_mi_control <- rewind_ipd_c$IPD


# combine datasets.

df_rewind_mi <- rbind(rewind_mi_control, rewind_mi_dula)

# survobj and then plot 

rewind_mi_s <- survfit(Surv(time,status) ~
                       treat, data = df_rewind_mi)


ggsurvplot(rewind_mi_s,
           data = df_rewind_mi,
           xlim = c(0,6),
           ylim = c(0,0.18),
           fun = "event",
           break.y.by = 0.03,
           palette = c("blue","red"))

# looks good, now save dataset.

write_csv( df_rewind_mi,
'F:\\GLP1_agonists\\data\\pooled_data\\mi_rewind.csv')


# sustain6 

# semaglutide arm 


sema <- 
  read.table("F:/GLP1_agonists/data/sustain6/mi_semaglutide.txt")

summary(sema)

# convert V2 to survival - 

sema$V2 <- 100 - sema$V2

summary(sema)

# now to extract information 

sustain_sema <- preprocess(dat = sema,
                          totalpts = 1648,
                          maxy = 100)


sustain_ipd_sema <- getIPD(prep = sustain_sema,
                       armID = 1)

plot(sustain_ipd_sema)

df_sustain_mi_sema <- sustain_ipd_sema$IPD

# sustain6 control arm 


control <- 
  read.table("F:/GLP1_agonists/data/sustain6/mi_control.txt")

summary(control)

# convert V2 to survival - 

control$V2 <- 100 - control$V2

summary(control)

# now to extract information 

sustain6_c <- preprocess(dat = control,
                       totalpts = 1649,
                       maxy = 100)


sustain6_ipd_c <- getIPD(prep = sustain6_c,
                       armID = 0)

plot(sustain6_ipd_c)

sustain6_mi_control <- sustain6_ipd_c$IPD

# combine datasets

mi_sustain6 <- rbind(sustain6_mi_control,
                     df_sustain_mi_sema)

# survobj and then plot

mi_sus_s <- survfit(Surv(time, status) ~ treat,
                    data = mi_sustain6)


ggsurvplot(mi_sus_s,
           data = mi_sustain6,
           fun = "event",
           ylim = c(0, 0.05),
           break.y.by = 0.01,
           xlim = c(0, 104),
           break.x.by = 8,
           palette = c("gray","blue"))

# save the dataset.

write_csv( mi_sustain6,
           'F:\\GLP1_agonists\\data\\pooled_data\\mi_sustain6.csv')

