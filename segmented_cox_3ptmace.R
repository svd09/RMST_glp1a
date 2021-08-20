##############################################
##  MULTISEGMENTED COX MODEL FOR THE PAPER  ##
##############################################
# 
# AUTHOR: SALIL DEO
# 
# EMAIL: svd14@case.edu
# 
# DATE: 2021-08-20
# 
# SCRIPT PURPOSE: CREATE A MULTISEGMENTED COX MODEL AT TIME POINTS FOR THE PAPER 
#   
# SCRIPT OTHER DETAILS: Create a multisegmented cox model for the paper. Pool using a random effects model and 
# calculate heterogeneity of pooled estimate. Plan to pool @ 12, 24, 36 and 48 months.
# Plan to calculate the CPH HR using the trials according to data available in the trials.
# This will also support the time varying effect observed by the delta-rmst method.
# It will also provide a relative effect estimate as opposed to the absolute effect estimate.


# SET WORKING DIRECTORY ---

# IMPORT PACKAGE LIBRARIES NEEDED ---

cat("IMPORTING THESE PACKAGES... \n\n", sep = "")
packages <- c("tidyverse","survival","metafor")
n_packages <- length(packages)

# install missing packages 

new.packages <- packages[!(packages %in% installed.packages())]

if(length(new.packages)){
  install.packages(new.packages)
}

# load all libraries needed

for(n in 1:n_packages){
  cat("Loading Library #", n, "of", n_packages, "...currently loading: ", packages[n], "\n", sep = "")
  lib_load <- paste("library(\"", packages[n], "\")", sep = "")
  eval(parse(text = lib_load))
}

# SETTING OPTIONS ---

cat("SETTING OPTIONS... \n\n", sep = "")
options(scipen = 999)
options(encoding = "UTF-8")

#---------------------------------------------------
# START CODING HERE
#---------------------------------------------------

# get the data from all the trials here



# get all the data.

amplitude <- 
  read_csv('F:/GLP1_agonists/data/pooled_data/3ptmace_amplitude.csv')



elixa <- 
  read_csv('F:/GLP1_agonists/data/pooled_data/3ptmace_elixa.csv')



sustain <- 
  read_csv('F:/GLP1_agonists/data/pooled_data/3ptmace_sustain6.csv')



harmony <- 
  read_csv('F:/GLP1_agonists/data/pooled_data/3ptmace_harmony.csv')



excsel <- 
  read_csv('F:/GLP1_agonists/data/pooled_data/3ptmace_excsel.csv')



leader <- 
  read_csv('F:/GLP1_agonists/data/pooled_data/3ptmace_leader.csv')



rewind <- 
  read_csv('F:/GLP1_agonists/data/pooled_data/3ptmace_rewind.csv')



pioneer <- 
  read_csv('F:/GLP1_agonists/data/pooled_data/3ptmace_pioneer.csv')


# now the data is loaded 

glimpse(leader)

# work with each trial separately and obtain the HR for each trial at the time points depending upon the data in each trial

# AMPLITUDE / HR at end of 12 months

amplitude12 = survival::survSplit(Surv(time, status) ~ ., 
                                  data = amplitude,
                                  cut = c(12),
                                  episode = "timegroup")

glimpse(amplitude12)


amplitude_s12 = coxph(Surv(tstart, time, status) ~ treat,
                      data = amplitude12[amplitude12$timegroup == 1, ])


summary(amplitude_s12)

# AMPLITUDE / HR at end of 24 months


amplitude24 = survival::survSplit(Surv(time, status) ~ ., 
                                  data = amplitude,
                                  cut = c(24),
                                  episode = "timegroup")

glimpse(amplitude24)


amplitude_s24 = coxph(Surv(tstart, time, status) ~ treat,
                      data = amplitude24[amplitude24$timegroup == 1, ])

summary(amplitude_s24)


# ELIXA


# ELIXA / HR at end of 12 months


elixa12 = survival::survSplit(Surv(time, status) ~ ., 
                                  data = elixa,
                                  cut = c(12),
                                  episode = "timegroup")

glimpse(elixa12)


elixa_s12 = coxph(Surv(tstart, time, status) ~ treat,
                      data = elixa12[elixa12$timegroup == 1, ])


summary(elixa_s12)


# ELIXA / HR at end of 24 months


elixa24 = survival::survSplit(Surv(time, status) ~ ., 
                              data = elixa,
                              cut = c(24),
                              episode = "timegroup")

glimpse(elixa24)


elixa_s24 = coxph(Surv(tstart, time, status) ~ treat,
                  data = elixa24[elixa24$timegroup == 1, ])


summary(elixa_s24)

# ELIXA at 36 months 



elixa36 = survival::survSplit(Surv(time, status) ~ ., 
                              data = elixa,
                              cut = c(36),
                              episode = "timegroup")

glimpse(elixa36)


elixa_s36 = coxph(Surv(tstart, time, status) ~ treat,
                  data = elixa36[elixa36$timegroup == 1, ])


summary(elixa_s36)

# SUSTAIN6

# SUSTAIN6 at 12 months 


sustain12 = survival::survSplit(Surv(time, status) ~ ., 
                              data = sustain,
                              cut = c(12),
                              episode = "timegroup")

glimpse(sustain12)


sustain_s12 = coxph(Surv(tstart, time, status) ~ treat,
                  data = sustain12[sustain12$timegroup == 1, ])


summary(sustain_s12)

# SUSTAIN6 at 24 months 

sustain24 = survival::survSplit(Surv(time, status) ~ ., 
                                data = sustain,
                                cut = c(24),
                                episode = "timegroup")

glimpse(sustain24)


sustain_s24 = coxph(Surv(tstart, time, status) ~ treat,
                    data = sustain24[sustain24$timegroup == 1, ])


summary(sustain_s24)

# Harmony Outcomes

# Harmony Outcomes at 12 months

harmony12 = survival::survSplit(Surv(time, status) ~ ., 
                                data = harmony,
                                cut = c(12),
                                episode = "timegroup")

glimpse(harmony12)


harmony_s12 = coxph(Surv(tstart, time, status) ~ treat,
                    data = harmony12[harmony12$timegroup == 1, ])


summary(harmony_s12)

# Harmony Outcomes at 24 months

harmony24 = survival::survSplit(Surv(time, status) ~ ., 
                                data = harmony,
                                cut = c(24),
                                episode = "timegroup")

glimpse(harmony24)


harmony_s24 = coxph(Surv(tstart, time, status) ~ treat,
                    data = harmony24[harmony24$timegroup == 1, ])


summary(harmony_s24)

# EXCSEL 

# EXCSEL at 12 months

excsel12 = survival::survSplit(Surv(time, status) ~ ., 
                                data = excsel,
                                cut = c(12),
                                episode = "timegroup")

glimpse(excsel12)


excsel_s12 = coxph(Surv(tstart, time, status) ~ treat,
                    data = excsel12[excsel12$timegroup == 1, ])


summary(excsel_s12)


# EXCSEL at 24 months


excsel24 = survival::survSplit(Surv(time, status) ~ ., 
                               data = excsel,
                               cut = c(24),
                               episode = "timegroup")

glimpse(excsel24)


excsel_s24 = coxph(Surv(tstart, time, status) ~ treat,
                   data = excsel24[excsel24$timegroup == 1, ])

summary(excsel_s24)


# EXCSEL at 36 months


excsel36 = survival::survSplit(Surv(time, status) ~ ., 
                               data = excsel,
                               cut = c(36),
                               episode = "timegroup")

glimpse(excsel36)


excsel_s36 = coxph(Surv(tstart, time, status) ~ treat,
                   data = excsel36[excsel36$timegroup == 1, ])

summary(excsel_s36)

# EXCSEL at 48 months 

excsel48 = survival::survSplit(Surv(time, status) ~ ., 
                               data = excsel,
                               cut = c(48),
                               episode = "timegroup")

glimpse(excsel48)


excsel_s48 = coxph(Surv(tstart, time, status) ~ treat,
                   data = excsel48[excsel48$timegroup == 1, ])

summary(excsel_s48)

# LEADER 

# LEADER at 12 months

leader12 = survival::survSplit(Surv(time, status) ~ ., 
                               data = leader,
                               cut = c(12),
                               episode = "timegroup")

glimpse(leader12)


leader_s12 = coxph(Surv(tstart, time, status) ~ treat,
                   data = leader12[leader12$timegroup == 1, ])


summary(leader_s12)

# LEADER at 24 months

leader24 = survival::survSplit(Surv(time, status) ~ ., 
                               data = leader,
                               cut = c(24),
                               episode = "timegroup")

glimpse(leader24)


leader_s24 = coxph(Surv(tstart, time, status) ~ treat,
                   data = leader24[leader24$timegroup == 1, ])


summary(leader_s24)

# LEADER at 36 months

leader36 = survival::survSplit(Surv(time, status) ~ ., 
                               data = leader,
                               cut = c(36),
                               episode = "timegroup")

glimpse(leader36)


leader_s36 = coxph(Surv(tstart, time, status) ~ treat,
                   data = leader36[leader36$timegroup == 1, ])


summary(leader_s36)

# LEADER at 48 months

leader48 = survival::survSplit(Surv(time, status) ~ ., 
                               data = leader,
                               cut = c(48),
                               episode = "timegroup")

glimpse(leader48)


leader_s48 = coxph(Surv(tstart, time, status) ~ treat,
                   data = leader48[leader48$timegroup == 1, ])


summary(leader_s48)

# REWIND 

# REWIND at 12 months


rewind12 = survival::survSplit(Surv(time, status) ~ ., 
                               data = rewind,
                               cut = c(12),
                               episode = "timegroup")

glimpse(rewind12)


rewind_s12 = coxph(Surv(tstart, time, status) ~ treat,
                   data = rewind12[rewind12$timegroup == 1, ])


summary(rewind_s12)

# REWIND at 24 months


rewind24 = survival::survSplit(Surv(time, status) ~ ., 
                               data = rewind,
                               cut = c(24),
                               episode = "timegroup")

glimpse(rewind24)


rewind_s24 = coxph(Surv(tstart, time, status) ~ treat,
                   data = rewind24[rewind24$timegroup == 1, ])


summary(rewind_s24)

# PIONEER6 

# PIONEER6 at 12 months


summary(pioneer)

pioneer12 = survival::survSplit(Surv(time, status) ~ ., 
                               data = pioneer,
                               cut = c(12),
                               episode = "timegroup")

glimpse(pioneer12)


pioneer_s12 = coxph(Surv(tstart, time, status) ~ treat,
                   data = pioneer12[pioneer12$timegroup == 1, ])


summary(pioneer_s12)

