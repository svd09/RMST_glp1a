# script for cv mortality 
# extract the data from the information 
# save the datasets after combining the data.

library("tidyverse")
library("survival")
library("IPDfromKM")
library("metaRMST")
library("broom")
library("survminer")

# get the data 1 trial at a time and then format it to extract information.

# leader 
# liraglutide arm - 

li <- read.table("F:/GLP1_agonists/data/leader/cvmort_liraglutide.txt")

summary(li)

# convert V2 to survival - 

li$V2 <- with(li, ifelse(V2 == -0.0479, 0, V2)) 

li$V2 <- 100 - li$V2

li$V2 <- with(li, ifelse(V2 > 100, 100, V2))

summary(li)

# now to extract information 

nrisk_l <- c(4668, 4641, 4599, 4558, 4505, 4445, 4382, 4322, 1723, 484)

trisk_l <- c(0,6, 12, 18, 24, 32, 36, 42, 48, 54)

leader_l <- preprocess(dat = li,
                             trisk = trisk_l,
                             nrisk = nrisk_l,
                             maxy = 100)




leader_ipd_li <- getIPD(prep = leader_l,
                             armID = 1)


plot(leader_ipd_li)

cv_mort_li <- leader_ipd_li$IPD


# control arm 


control <- 
read.table("F:/GLP1_agonists/data/leader/cvmort_control.txt")

summary(control)

# convert V2 to survival - 

control$V2 <- 100 - control$V2

summary(control)

# now to extract information 

nrisk_c <- c(4672, 4648, 4601, 4479, 4407, 4338, 4267, 1709, 465)

trisk_c <- c(0,6, 12, 18, 24, 32, 36, 42, 48, 54)

leader_c <- preprocess(dat = control,
                       trisk = trisk_c,
                       nrisk = nrisk_c,
                       maxy = 100)




leader_ipd_c <- getIPD(prep = leader_c,
                        armID = 0)


plot(leader_ipd_c)

cvmort_control <- leader_ipd_c$IPD


# combine data and then plot 

df_leader_cvmort <- rbind(cv_mort_li, cvmort_control)

cvmort_leader_s <- 
  survfit(Surv(time, status) ~ treat, data = df_leader_cvmort)


ggsurvplot(cvmort_leader_s, 
           fun = "event",
           ylim = c(0,0.20),
           xlim = c(0,54),
           break.x.y = 6,
           break.y.by = 0.05,
           censor.size = 0,
           palette = c("gray","blue2"))

# time already in months.

write_csv(df_leader_cvmort,
  'F:\\GLP1_agonists\\data\\pooled_data\\cvmort_leader.csv'
)

# sustain 6 

sema <- read.table("F:/GLP1_agonists/data/sustain6/cvmort_semaglutide.txt")

summary(sema)

# convert V2 to survival - 

sema$V2 <- 100 - sema$V2

summary(sema)

# now to extract information 

nrisk_s <- c(1648, 1634, 1627, 1617, 1607, 1589, 1579)

trisk_s <- c(0,16, 32, 48, 64, 80, 96)

sustain6_s <- preprocess(dat = sema,
                       trisk = trisk_s,
                       nrisk = nrisk_s,
                       maxy = 100)


sustain6_ipd_sema <- getIPD(prep = sustain6_s,
                        armID = 1)


plot(sustain6_ipd_sema)

cv_mort_sema <- sustain6_ipd_sema$IPD


# now to get control arm data 


control <- 
read.table("F:/GLP1_agonists/data/sustain6/cvmort_control.txt")

summary(control)

# convert V2 to survival - 

control$V2 <- 100 - control$V2

summary(control)

# now to extract information 

nrisk_c <- c(1649, 1637, 1623, 1617, 1600, 1584, 1566)

trisk_c <- c(0,16, 32, 48, 64, 80, 96)

sustain6_c <- preprocess(dat = control,
                         trisk = trisk_c,
                         nrisk = nrisk_c,
                         maxy = 100)


sustain6_ipd_control <- getIPD(prep = sustain6_c,
                            armID = 0)


plot(sustain6_ipd_control)

cv_mort_control <- sustain6_ipd_control$IPD


# combine the dataset and then plot again 

cvmort_sustain6 <- rbind(cv_mort_sema, cv_mort_control)

surv_sustain6 <- survfit(Surv(time, status) ~ treat, 
                         data = cvmort_sustain6)

# plot graph 


ggsurvplot(surv_sustain6,
           data = cvmort_sustain6,
           fun = "event",
           ylim = c(0,0.05),
           break.y.by = 0.01,
           xim = c(0, 104),
           break.x.by = 8,
           palette = c("gray","blue2"),
           censor.size = 0)

# now to convert weeks to months
# save with time in months 

cvmort_sustain6$time <- cvmort_sustain6$time/4

# save dataaset

write_csv(cvmort_sustain6,
  'F:\\GLP1_agonists\\data\\pooled_data\\cvmort_sustain6.csv')


# EXCSEL 


exe <- 
read.table("F:/GLP1_agonists/data/excsel/cvmort_exenatide.txt")

summary(exe)

# convert V2 to survival - 

exe$V2 <- 100 - exe$V2

summary(exe)

# now to extract information 

nrisk_e <- c(7356, 7234, 6433, 4095, 2698, 892)

trisk_e <- c(0,1,2,3,4,5)

excsel_exe <- preprocess(dat = exe,
                         trisk = trisk_e,
                         nrisk = nrisk_e,
                         maxy = 100)


excsel_ipd_exe <- getIPD(prep = excsel_exe,
                            armID = 1)


plot(excsel_ipd_exe)

cv_mort_exe <- excsel_ipd_exe$IPD

# control arm of EXCSEL 

control <- 
  read.table("F:/GLP1_agonists/data/excsel/cvmort_control.txt")

summary(control)

# convert V2 to survival - 

control$V2 <- 100 - control$V2

summary(control)

# now to extract information 

nrisk_c <- c(7396, 7278, 6470, 4091, 2666, 907)

trisk_c <- c(0,1,2,3,4,5)

excsel_c <- preprocess(dat = control,
                         trisk = trisk_c,
                         nrisk = nrisk_c,
                         maxy = 100)


excsel_ipd_c <- getIPD(prep = excsel_c,
                         armID = 0)


plot(excsel_ipd_c)

cv_mort_control <- excsel_ipd_c$IPD

# combine dataset and then plot 

cvmort_excsel <- rbind(cv_mort_exe, cv_mort_control)


cvmort_excsel_s <- survfit(Surv(time, status) ~ treat,
                           data = cvmort_excsel)

ggsurvplot(cvmort_excsel_s,
           data = cvmort_excsel,
           xlim = c(0,5),
            ylim = c(0,0.18),
            break.y.by = 0.03,
           censor.size = 0,
          fun = "event",
          palette = c("blue","red"),
          linetype = c(2,1))


# convert time to months
# save dataset

cvmort_excsel$time <- cvmort_excsel$time*12

write_csv(cvmort_excsel,
  'F:\\GLP1_agonists\\data\\pooled_data\\cvmort_excsel.csv'
)


# Harmony Outcomes 


albi <- 
  read.table("F:/GLP1_agonists/data/Harmony/cvmort_albiglutide.txt")

summary(albi)

# convert V2 to survival - 

albi$V2 <- 100 - albi$V2

summary(albi)

# now to extract information 

nrisk_a <- c(4731, 4681, 4611, 4379, 3274, 2234, 1121)

trisk_a <- c(0, 4, 8, 12, 16, 20, 24 )

harmony_a <- preprocess(dat = albi,
                         trisk = trisk_a,
                         nrisk = nrisk_a,
                         maxy = 100)


harmony_ipd_a <- getIPD(prep = harmony_a,
                         armID = 1)


plot(harmony_ipd_a)

cv_mort_albi <- harmony_ipd_a$IPD


# control arm in Harmony Outcomes 


control <- 
read.table("F:/GLP1_agonists/data/Harmony/cvmort_control.txt")

summary(control)

# convert V2 to survival - 

control$V2 <- 100 - control$V2

summary(control)

# now to extract information 

nrisk_c <- c(4732, 4662, 4580, 4373, 3245, 2261, 1121)

trisk_c <- c(0, 4, 8, 12, 16, 20, 24 )

harmony_c <- preprocess(dat = control,
                        trisk = trisk_c,
                        nrisk = nrisk_c,
                        maxy = 100)


harmony_ipd_c <- getIPD(prep = harmony_c,
                        armID = 0)


plot(harmony_ipd_c)

cv_mort_harmony_c <- harmony_ipd_c$IPD

# combine dataasets

cvmort_harmony <- rbind(cv_mort_albi, cv_mort_harmony_c)

cvmort_hs <- survfit(Surv(time, status) ~ treat, 
                     data = cvmort_harmony)

ggsurvplot(
  cvmort_hs,
  data = cvmort_harmony,
  fun = "event",
  xlim = c(0,28),
  ylim = c(0, 0.18),
  break.y.by = 0.02,
  break.x.by = 4,
  censor.size = 0,
  palette = c("blue2","red")
)


# convert time to months from weeks

cvmort_harmony$time <- cvmort_harmony$time*4

write_csv(
  cvmort_harmony,
  'F:\\GLP1_agonists\\data\\pooled_data\\cvmort_harmony.csv'
)

# REWIND


dula <- 
read.table("F:/GLP1_agonists/data/rewind/cvmort_dulaglutide.txt")

summary(dula)

dula$V1 <- with(dula, ifelse(V1 < 0, 0 , V1))

dula$V2 <- with(dula, ifelse(V2 < 0, 0, V2))

# convert V2 to survival - 

dula$V2 <- 100 - dula$V2

summary(dula)

# now to extract information 

nrisk_d <- c(4949, 4866, 4773, 4663, 4556, 3887, 807)

trisk_d <- c(0, 1, 2, 3, 4, 5, 6)

rewind_d <- preprocess(dat = dula,
                        trisk = trisk_d,
                        nrisk = nrisk_d,
                        maxy = 100)


rewind_ipd_d <- getIPD(prep = rewind_d,
                        armID = 1)


plot(rewind_ipd_d)

cvmort_rewind_d <- rewind_ipd_d$IPD

# control arm of rewind 

control <- 
  
read.table(
  "F:/GLP1_agonists/data/rewind/cvmort_control.txt"
)

summary(control)

control$V1 <- with(control, ifelse(V1 < 0, 0 , V1))

# convert V2 to survival - 

control$V2 <- 100 - control$V2

summary(control)

# now to extract information 

nrisk_c <- c(4952, 4854, 4748, 4617, 4499, 3813, 802)

trisk_c <- c(0, 1, 2, 3, 4, 5, 6)

rewind_c <- preprocess(dat = control,
                       trisk = trisk_c,
                       nrisk = nrisk_c,
                       maxy = 100)


rewind_ipd_c <- getIPD(prep = rewind_c,
                       armID = 0)


plot(rewind_ipd_c)

cvmort_rewind_c <- rewind_ipd_c$IPD

# combine dataset

cvmort_rewind <- rbind(cvmort_rewind_d, cvmort_rewind_c)

# create survobject and plot 

cvmort_rewind_s <- survfit(Surv(time, status) ~ treat, 
                           data = cvmort_rewind)

plot(cvmort_rewind_s)

ggsurvplot(cvmort_rewind_s,
           data = cvmort_rewind,
           fun = "event",
           ylim = c(0, 0.18),
           censor.size = 0,
           break.y.by = 0.03)

# convert years to months 

cvmort_rewind$time <- cvmort_rewind*12

write_csv(
  cvmort_rewind,
  'F:\\GLP1_agonists\\data\\pooled_data\\cvmort_rewind.csv'
)

# PIONEER 6 

psema <- 
  read.table(
    "F:/GLP1_agonists/data/pioneer6/cvmort_semaglutide.txt")


summary(psema)


psema$V1 <- with(psema, ifelse(V1 < 0, 0, V1))

psema$V2 <- with(psema, ifelse(V2 < 0, 0, V2))

# convert V2 to survival - 

psema$V2 <- 100 - psema$V2

summary(psema)



# now to extract information 

pioneer6_sema <- preprocess(dat = psema,
                            totalpts = 1591,
                         maxy = 100)


pioneer6_ipd_sema <- getIPD(prep = pioneer6_sema,
                         armID = 1)


plot(pioneer6_ipd_sema)

cvmort_pioneer6_sema <- pioneer6_ipd_sema$IPD


# control arm for pioneer6 


p_control <- 
  read.table(
    "F:/GLP1_agonists/data/pioneer6/cvmort_control.txt")


summary(p_control)


p_control$V1 <- with(p_control, ifelse(V1 < 0, 0, V1))

# convert V2 to survival - 

p_control$V2 <- 100 - p_control$V2

summary(p_control)


# now to extract information 

pioneer6_control <- preprocess(dat = p_control,
                            totalpts = 1592,
                            maxy = 100)


pioneer6_ipd_control <- getIPD(prep = pioneer6_control,
                            armID = 0)


cvmort_pioneer6_control <- pioneer6_ipd_control$IPD

cvmort_pioneer6 <- rbind(cvmort_pioneer6_sema, cvmort_pioneer6_control)

# plot graph 

cvmort_pioneer6_s <- survfit(Surv(time, status) ~ treat, 
                             data = cvmort_pioneer6)


ggsurvplot(cvmort_pioneer6_s,
           data = cvmort_pioneer6,
           xlim = c(0, 83),
           ylim = c(0, 0.1),
           break.x.by = 9,
           break.y.by = 0.01,
           censor.size = 0,
           palette = c("gray","blue2"),
          fun = "event")

# convert weeks to months 
# save combined dataset

cvmort_pioneer6$time <- cvmort_pioneer6$time/4

write_csv(cvmort_pioneer6,
          'F:\\GLP1_agonists\\data\\pooled_data\\cvmort_pioneer6.csv')


