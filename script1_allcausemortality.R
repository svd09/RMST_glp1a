# script for all-cause mortality 
# extract the data from the information 
# save the datasets after combining the data.

library("tidyverse")
library("survival")
library("IPDfromKM")
library("metaRMST")
library("broom")
library("survminer")

# only 2 trials have information regarding all-cause mortality.

# LEADER 

# liraglutide arm 


li <- read.table("F:/GLP1_agonists/data/leader/allmort_liraglutide.txt")

summary(li)

# convert V2 to survival - 


li$V2 <- 100 - li$V2

summary(li)

# now to extract information 

nrisk_l <- c(4668, 4641, 4599, 4558, 4505, 4445, 4382, 4322, 1723)

trisk_l <- c(0,6, 12, 18, 24, 32, 36, 42, 48)

leader_l <- preprocess(dat = li,
                       trisk = trisk_l,
                       nrisk = nrisk_l,
                       maxy = 100)




leader_ipd_li <- getIPD(prep = leader_l,
                        armID = 1)


plot(leader_ipd_li)

all_mort_li <- leader_ipd_li$IPD

# LEADER control arm 



control <- 
  read.table("F:/GLP1_agonists/data/leader/allmort_control.txt")

summary(control)

# convert V2 to survival - 

control$V2 <- 100 - control$V2

summary(control)

# now to extract information 

nrisk_c <- c(4672, 4648, 4601, 4546, 4479, 4407, 4338, 4268, 1709)

trisk_c <- c(0,6, 12, 18, 24, 32, 36, 42, 48)

leader_c <- preprocess(dat = control,
                       trisk = trisk_c,
                       nrisk = nrisk_c,
                       maxy = 100)


leader_ipd_c <- getIPD(prep = leader_c,
                       armID = 0)


plot(leader_ipd_c)

allmort_control <- leader_ipd_c$IPD


# combine data and then plot


df_leader_allmort <- rbind(all_mort_li, allmort_control)



allmort_leader_s <- 
  survfit(Surv(time, status) ~ treat, data = df_leader_allmort)


ggsurvplot(allmort_leader_s, 
           fun = "event",
           ylim = c(0,0.20),
           xlim = c(0,54),
           break.x.y = 6,
           break.y.by = 0.05,
           censor.size = 0,
           palette = c("gray","blue2"),
           break.x.by = 6)


# save dataset.


  write_csv(df_leader_allmort,
            'F:\\GLP1_agonists\\data\\pooled_data\\allmort_leader.csv')

  
# EXCSEL - 
  
exe <- 
    read.table("F:/GLP1_agonists/data/excsel/allmort_exenatide.txt")
  
  summary(exe)
  
  # convert V2 to survival - 
  
exe$V2 <- 100 - exe$V2
  
  summary(exe)
  
  # now to extract information 
  
nrisk_e <- c(7356, 7234, 6433, 4095, 2698, 907)
  
  trisk_e <- c(0,1,2,3,4,5)
  
  excsel_exe <- preprocess(dat = exe,
                           trisk = trisk_e,
                           nrisk = nrisk_e,
                           maxy = 100)
  
  
  excsel_ipd_exe <- getIPD(prep = excsel_exe,
                           armID = 1)
  
  
  plot(excsel_ipd_exe)
  
all_mort_exe <- excsel_ipd_exe$IPD
  
# control arm of EXCSEL 
  
  control <- 
    read.table("F:/GLP1_agonists/data/excsel/allmort_control.txt")
  
  summary(control)
  
  # convert V2 to survival - 
  
  control$V2 <- 100 - control$V2
  
  summary(control)
  
  control$V2 <- with(control, ifelse(V2 > 100, 100, V2))
  
  # now to extract information 
  
  nrisk_c <- c(7396, 7278, 6470, 4091, 2666, 892)
  
  trisk_c <- c(0,1,2,3,4,5)
  
  excsel_c <- preprocess(dat = control,
                         trisk = trisk_c,
                         nrisk = nrisk_c,
                         maxy = 100)
  
  
  excsel_ipd_c <- getIPD(prep = excsel_c,
                         armID = 0)
  
  
  plot(excsel_ipd_c)
  
all_mort_control <- excsel_ipd_c$IPD
  
  # combine dataset and then plot 
  
  allmort_excsel <- rbind(all_mort_exe, all_mort_control)
  
  
  allmort_excsel_s <- survfit(Surv(time, status) ~ treat,
                             data = allmort_excsel)
  
  ggsurvplot(allmort_excsel_s,
             data = allmort_excsel,
             xlim = c(0,5),
             ylim = c(0,0.18),
             break.y.by = 0.03,
             censor.size = 0,
             fun = "event",
             palette = c("blue","red"),
             linetype = c(2,1))
  
  
  # convert time to months
  # save dataset
  
  allmort_excsel$time <- allmort_excsel$time*12
  
  write_csv(allmort_excsel,
            'F:\\GLP1_agonists\\data\\pooled_data\\allmort_excsel.csv'
  )
  
  
  