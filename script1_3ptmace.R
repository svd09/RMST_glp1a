# Study: Pooled analysis of d-RMST for GLP1 agonists
# 
# This script contains code to extract IPD data from the curves for 3pt MACE
# Create a combined dataset and then save it to this folder.

library(easypackages)

libraries(c("tidyverse","IPDfromKM",'survival',
    "survRMST","flexsurv","broom","rstpm2","survminer"))
    

# ELIXA - 
# for lixisenatide 

nrisk_l <- c(3034, 2785, 1558, 484)

trisk_l <- c(0,1,2,3.4)

lix <- read.table("F:\\GLP1_agonists\\data\\elixa\\3pmace_lixisenatide.txt")

lix$V2 <- abs(lix$V2)

lix <- lix %>% arrange(V1, V2)

# y axis data needs to be as KM plot not cumulative event plot.
# reformat this if needed.


lix$V2 <- 100 - lix$V2

summary(lix)

res_lix <- preprocess(dat = lix,
                      trisk = trisk_l,
                      nrisk = nrisk_l,
                      maxy = 100)




elixa_ipd_lix <- getIPD(prep = res_lix,
                        armID = 1)


summary(elixa_ipd_lix)

plot(elixa_ipd_lix)


df_elixa_lix <- elixa_ipd_lix$IPD


# for control 

nrisk_c <- c(3034,2759,1566,476)

trisk_c <- c(0,1,2,3.4)

control <- read.table("F:\\GLP1_agonists\\data\\elixa\\3pmace_control.txt")


control$V2 <- 100 - control$V2

elixa_control <- preprocess(dat = control,
                            trisk = trisk_c,
                            nrisk = nrisk_c,
                            maxy = 100)


elixa_ipd_control <- getIPD(prep = elixa_control,
                            armID = 0)

df_elixa_control <- elixa_ipd_control$IPD

# combine both df 

elixa_df <- rbind(df_elixa_control, df_elixa_lix)


# change the time from years to months 

elixa_df$time <- elixa_df$time*12

# create survival obj and then plot.

elixa_plot <- survfit(Surv(time, status) ~ factor(treat), data = elixa_df)


# create plot 

ggsurvplot(elixa_plot, fun = "event",
           data = elixa_df, censor.size = 0,
           break.x.by = 12, ylim = c(0, 0.2))


# plot looks good compared to the original. 
# for the analysis , am going to use months.
# so no change needed here.



write_csv(elixa_df, 
          "F:\\GLP1_agonists\\data\\pooled_data\\3ptmace_elixa.csv")

# EXCSEL trial 



# for exenatide

nrisk_e <- c(7356, 7101, 6893, 6580, 5912, 4475, 3595, 3053, 2281, 1417, 727)

trisk_e <- c(0,0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5)

exe <- read.table("F:\\GLP1_agonists\\data\\excsel\\3ptmace_exenatide.txt")


# y axis data needs to be as KM plot not cumulative event plot.
# reformat this if needed.

exe$V2 <- 100 - exe$V2

summary(exe)

res_exe <- preprocess(dat = exe,
                      trisk = trisk_e,
                      nrisk = nrisk_e,
                      maxy = 100)




excsel_ipd_exe <- getIPD(prep = res_exe,
                         armID = 1)


summary(excsel_ipd_exe)

plot(excsel_ipd_exe)


df_excsel_exe <- excsel_ipd_exe$IPD


# for control 

nrisk_c <- c(7396, 7120, 6897, 6565, 5908, 4468, 3365, 2961, 2209, 1366, 687)

trisk_c <- c(0,0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5)

control <- read.table("F:\\GLP1_agonists\\data\\excsel\\3ptmace_control.txt")


control$V2 <- 100 - control$V2

excsel_control <- preprocess(dat = control,
                             trisk = trisk_c,
                             nrisk = nrisk_c,
                             maxy = 100)


excsel_ipd_control <- getIPD(prep = excsel_control,
                             armID = 0)

df_excsel_control <- excsel_ipd_control$IPD

# combine both df 

excsel_df <- rbind(df_excsel_control, df_excsel_exe)


# create survival obj and then plot.

excsel_plot <- survfit(Surv(time, status) ~ factor(treat), data = excsel_df)


# create plot 

ggsurvplot(excsel_plot, fun = "event",
           data = excsel_df, censor.size = 0, break.y.by = 0.03)

# plot looks good compared to the original.
# for pooling purposes, am going to convert all time to months.
# so changed to mnonths

excsel_df$time <- excsel_df$time*12

# save this data.


write_csv(excsel_df, 
          "F:\\GLP1_agonists\\data\\pooled_data\\3ptmace_excsel.csv")

# Harmony outcomes 

# Albiglutide arm


nrisk_a <- c(4731, 4239, 1064)

trisk_a <- c(0,1,2)

alb <- read.table("F:\\GLP1_agonists\\data\\Harmony\\3ptmace_albiglutide.txt")

glimpse(alb)

summary(alb)

# convert survival from 100

alb$V2 <- 100 - alb$V2

summary(alb)


harmony_a <- preprocess(dat = alb,
                             trisk = trisk_a,
                             nrisk = nrisk_a,
                             maxy = 100)




harmony_ipd_alb <- getIPD(prep = harmony_a,
                             armID = 1)


plot(harmony_ipd_alb)


df_harmony_alb <- harmony_ipd_alb$IPD

# control arm 


nrisk_c <- c(4731,4208, 1030)

trisk_c <- c(0,1,2)

h_con <- read.table("F:\\GLP1_agonists\\data\\Harmony\\3ptmace_control.txt")

glimpse(h_con)

summary(h_con)

# convert survival from 100

h_con$V2 <- 100 - h_con$V2

summary(h_con)


harmony_c <- preprocess(dat = h_con,
                        trisk = trisk_c,
                        nrisk = nrisk_c,
                        maxy = 100)


harmony_ipd_con <- getIPD(prep = harmony_c,
                          armID = 0)


plot(harmony_ipd_con)


df_harmony_con <- harmony_ipd_con$IPD

# combine the df 

df_harmony <- rbind(df_harmony_alb, df_harmony_con)

# survobject and then plot 
# convert to months 

df_harmony$time <- df_harmony$time*12


harmony_s <- survfit(Surv(time, status) ~ treat, data = df_harmony)

ggsurvplot(harmony_s, fun = "event",
           break.x.by = 4, break.y.by = 0.02,
           censor.size = 0, ylim = c(0,0.16),
           palette = c("blue","red"))

# time already in months, so no change needed.
# now to save the dataset.


write_csv(df_harmony, 
          "F:\\GLP1_agonists\\data\\pooled_data\\3ptmace_harmony.csv")


# LEADER 

# Liraglutide 


nrisk_lir <- c(4668, 4593,4496, 4400, 4280, 4172, 4072, 3982, 1562)

trisk_lir <- c(0,0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4)

lir <- read.table("F:\\GLP1_agonists\\data\\leader\\3ptmace_liraglutide.txt")

glimpse(lir)

summary(lir)

# convert survival from 100

lir$V2 <- 100 - lir$V2

summary(lir)


leader_l <- preprocess(dat = lir,
                        trisk = trisk_lir,
                        nrisk = nrisk_lir,
                        maxy = 100)




leader_ipd_l <- getIPD(prep = leader_l,
                          armID = 1)


plot(leader_ipd_l)


df_leader_lir <- leader_ipd_l$IPD

# control arm of leader 



nrisk_c <- c(4672, 4588,4473, 4352, 4237, 4123, 4010, 3914, 1543)

trisk_c <- c(0,0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4)

leader_c <- read.table("F:\\GLP1_agonists\\data\\leader\\3ptmace_control.txt")

glimpse(leader_c)

summary(leader_c)

# convert survival from 100

leader_c$V2 <- 100 - leader_c$V2

summary(leader_c)


leader_control <- preprocess(dat = leader_c,
                       trisk = trisk_c,
                       nrisk = nrisk_c,
                       maxy = 100)




leader_ipd_c <- getIPD(prep = leader_control,
                       armID = 0)


plot(leader_ipd_c)


df_leader_control <- leader_ipd_c$IPD


# combine data and then convert time to months 
# time converted to months, plan to pool data as months


df_leader <- rbind(df_leader_lir, df_leader_control)

df_leader$time <- df_leader$time*12

# survobj and plot 

leader_s <- survfit(Surv(time, status) ~ treat, data = df_leader)


ggsurvplot(leader_s, fun = "event",
           break.x.by = 6, break.y.by = 0.05,
           censor.size = 0, ylim = c(0,0.2),
           palette = c("blue","red"),
           xlim = c(0,54))

# looks good 
# save data

write_csv(df_leader, 
          "F:\\GLP1_agonists\\data\\pooled_data\\3ptmace_leader.csv")

# pioneer6

# semaglutide


sema <- 
read.table("F:\\GLP1_agonists\\data\\pioneer6\\3ptmace_semaglutide.txt")

glimpse(sema)

summary(sema)

# convert survival from 100

sema$V2 <- with(sema, ifelse(V2 < 0, -1*V2, V2))

summary(sema)

sema$V2 <- 100 - sema$V2

summary(sema)


pioneer_sema <- preprocess(dat = sema,
                             totalpts = 1591,
                             maxy = 100)


pioneer_ipd_sema <- getIPD(prep = pioneer_sema,
                       armID = 1)


plot(pioneer_ipd_sema)

df_pioneer_sema <- pioneer_ipd_sema$IPD


# pioneer6 control 



p_control <- 
  read.table("F:\\GLP1_agonists\\data\\pioneer6\\3ptmace_control.txt")

glimpse(p_control)

summary(p_control)

# convert survival from 100

p_control$V2 <- with(p_control, ifelse(V2 < 0, -1*V2, V2))

summary(p_control)

p_control$V2 <- 100 - p_control$V2

summary(p_control)


pioneer_control <- preprocess(dat = p_control,
                           totalpts = 1592,
                           maxy = 100)


pioneer_ipd_control <- getIPD(prep = pioneer_control,
                           armID = 0)


plot(pioneer_ipd_control)

df_pioneer_control <- pioneer_ipd_control$IPD


df_pioneer <- rbind(df_pioneer_sema, df_pioneer_control)

summary(df_pioneer)

# create survobj and then look at plot

df_pioneer$time <- (df_pioneer$time*365)/7

# convert time to weeks 

pioneer_s <- survfit(Surv(time, status) ~ treat, data = df_pioneer)


ggsurvplot(pioneer_s, fun = "event",
           break.x.by = 9, break.y.by = 0.01,
           censor.size = 0, ylim = c(0,0.1),
           palette = c("gray","blue"),
           xlim = c(0,83))

# figure looks fine, save dataset in weeks.
# am going to convert all time to months now.

df_pioneer$time <- df_pioneer$time/4

# converted to months, save dataaset


write_csv(df_pioneer, 
  "F:\\GLP1_agonists\\data\\pooled_data\\3ptmace_pioneer.csv")

# rewind 

# dulaglutide 


dula <- 
  read.table("F:\\GLP1_agonists\\data\\rewind\\3ptmace_dulaglutide.txt")

glimpse(dula)

summary(dula)

# convert survival from 100


dula$V2 <- with(dula, ifelse(V2 < 0, -1*V2, V2))

dula$V1 <- with(dula, ifelse(V1 < 0, -1*V1, V1))

summary(dula)

dula$V2 <- 100 - dula$V2

summary(dula)

nrisk_d <- c(4949, 4815, 4670, 4521, 4369, 3686, 741)
trisk_d <- c(0,1,2,3,4,5,6)

rewind_dula <- preprocess(dat = dula,
                          nrisk = nrisk_d,
                          trisk = trisk_d,
                           maxy = 100)


rewind_ipd_dula <- getIPD(prep = rewind_dula,
                           armID = 1)


plot(rewind_ipd_dula)

df_rewind_dula <- rewind_ipd_dula$IPD


# rewind control 


nrisk_c <- c(4952, 4791, 4625, 4437, 4275, 3575, 742)

trisk_c <- c(0,1,2,3,4, 5, 6)

control <- 
read.table("F:\\GLP1_agonists\\data\\rewind\\3ptmace_control.txt")


summary(control)

control$V2 <- 100 - control$V2

rewind_control <- preprocess(dat = control,
                             trisk = trisk_c,
                             nrisk = nrisk_c,
                            # totalpts = 4952,
                            maxy = 100)

rewind_ipd_control <- getIPD(prep = rewind_control,
                            armID = 0)

df_rewind_control <- rewind_ipd_control$IPD

# combine datasets 

df_rewind <- rbind(df_rewind_dula, df_rewind_control)

# survobj and then plot 

rewind_s <- survfit(Surv(time, status) ~ treat, data = df_rewind)


ggsurvplot(rewind_s, fun = "event",
           break.x.by = 1, break.y.by = 0.03,
           censor.size = 0, ylim = c(0,0.18),
           palette = c("blue","red"),
           xlim = c(0,6))

# convert years to months

df_rewind$time <- df_rewind$time*4

# save dataset.

write_csv(df_rewind, 
"F:\\GLP1_agonists\\data\\pooled_data\\3ptmace_rewind.csv")

# sustain 

# semaglutide 

semag <- 
  read.table("F:\\GLP1_agonists\\data\\sustain6\\3ptmace_semaglutide.txt")


summary(semag)

semag$V2 <- 100 - semag$V2

sustain6_semag <- preprocess(dat = semag,totalpts = 1648,
                             maxy = 100)

sustain6_ipd_semag <- getIPD(prep = sustain6_semag,
                             armID = 1)

df_sustain6_semag <- sustain6_ipd_semag$IPD

# control for sustain6

control <- 
  read.table("F:\\GLP1_agonists\\data\\sustain6\\3ptmace_control.txt")


summary(control)

control$V2 <- 100 - control$V2

sustain6_control <- preprocess(dat = control,
                               totalpts = 1648,
                             maxy = 100)

sustain6_ipd_control <- getIPD(prep = sustain6_control,
                             armID = 0)

df_sustain6_control <- sustain6_ipd_control$IPD


df_sustain6 <- rbind(df_sustain6_semag, df_sustain6_control)

# create surv obj and then plot 
# convert to weeks for the plot.

df_sustain6$time <- (df_sustain6$time*365)/7



sustain_s <- survfit(Surv(time, status) ~ treat, 
                     data = df_sustain6)


ggsurvplot(sustain_s, fun = "event",
           break.x.by = 8, break.y.by = 0.01,
           censor.size = 0, ylim = c(0,0.10),
           palette = c("gray","skyblue"),
           xlim = c(0,104))

# now to convert to months for pooling data.

df_sustain6$time <- df_sustain6$time/4

write_csv(df_sustain6,
          "F:\\GLP1_agonists\\data\\pooled_data\\3ptmace_sustain6.csv" )


# amplitude_o


efpe <- 
  read.table("F:\\GLP1_agonists\\data\\amplitude_o\\3ptmace_efpeglenatide.txt")

summary(efpe)

efpe$V2 <- 100 - efpe$V2

trisk_e <- c(0,1,2,3,4)

nrisk_e <- c(2717, 2644, 2587, 2503, 594)



amp_efpe <- preprocess(dat = efpe,
                       trisk = trisk_e,
                       nrisk = nrisk_e,
                             maxy = 100)

amp_ipd_efpe <- getIPD(prep = amp_efpe,
                             armID = 1)

df_amp_efpe <- amp_ipd_efpe$IPD

# amplitude control arm 



control <- 
  read.table("F:\\GLP1_agonists\\data\\amplitude_o\\3ptmace_control.txt")

summary(control)

control$V2 <- 100 - control$V2

trisk_c <- c(0,1,2,3,4)

nrisk_c <- c(1359, 1311, 1258, 1213, 278)



amp_control <- preprocess(dat = control,
                       trisk = trisk_c,
                       nrisk = nrisk_c,
                       maxy = 100)

amp_ipd_control <- getIPD(prep = amp_control,
                       armID = 0)

df_amp_control <- amp_ipd_control$IPD


# create dataaset 

df_amplitude <- rbind(df_amp_efpe, df_amp_control)

# convert time to months 

df_amplitude$time <- df_amplitude$time*12

# survobj and then plot 


amp_s <- survfit(Surv(time, status) ~ treat, data = df_amplitude)



ggsurvplot(amp_s, fun = "event",
           break.x.by = 6, break.y.by = 0.03,
           censor.size = 0, ylim = c(0,0.15),
           palette = c("blue3","red"),
           linetype  = c(2,1),
           xlim = c(0,24))



# save dataset.

write_csv(
  df_amplitude,
  "F:\\GLP1_agonists\\data\\pooled_data\\3ptmace_amplitude.csv" )


