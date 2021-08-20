# add into the supplemental section 
# do plotting of KM with RP model to graphically demonstrate model fit.


library(easypackages)

libraries(c("tidyverse","IPDfromKM",'survival',
            "flexsurv","broom","rstpm2","survminer",
            "metaRMST", "ggthemes","ckbplotr","broom",
            "jskm"))


# get all the data.

amplitude <- 
  read_csv('F:/GLP1_agonists/data/pooled_data/3ptmace_amplitude.csv')

amplitude$trialID <- 1

elixa <- 
  read_csv('F:/GLP1_agonists/data/pooled_data/3ptmace_elixa.csv')

elixa$trialID <- 2

sustain <- 
  read_csv('F:/GLP1_agonists/data/pooled_data/3ptmace_sustain6.csv')

sustain$trialID <- 3

harmony <- 
  read_csv('F:/GLP1_agonists/data/pooled_data/3ptmace_harmony.csv')

harmony$trialID <- 4

excsel <- 
  read_csv('F:/GLP1_agonists/data/pooled_data/3ptmace_excsel.csv')

excsel$trialID <- 5

leader <- 
  read_csv('F:/GLP1_agonists/data/pooled_data/3ptmace_leader.csv')

leader$trialID <- 6

rewind <- 
  read_csv('F:/GLP1_agonists/data/pooled_data/3ptmace_rewind.csv')

rewind$trialID <- 7

pioneer <- 
  read_csv('F:/GLP1_agonists/data/pooled_data/3ptmace_pioneer.csv')

pioneer$trialID <- 8


# fit RP model with 3 spline terms and TVC to fit the model for each study.

# Amplitude 


kmfit <- survfit(Surv(time, status) ~ treat, data = amplitude)

amp <- rstpm2::stpm2(Surv(time, status) ~ treat, 
                     df = 3,tvc = list(treat = 1),  data = amplitude)

amp_fit <- predict(amp, newdata = data.frame(treat = c(0:1)),
                   grid = T, full = T, se.fit = T,
                   type = "surv")

amp_fit$cif <- 1 - amp_fit$Estimate

amp.p <- jskm::jskm(sfit = kmfit,
           cumhaz = T,legend = F,
           linecols = "black")

amp.p2 = p + ylim(0,0.25)



# convert amp to CIF 

amp.p3 = amp.p2 + geom_line(data = amp_fit, aes(x = time, y = cif, 
                                        color = factor(treat)),
                    size = 1) +
  labs(x = "Months since Randomisation", 
  y = "Cumulative Incidence ") + ggtitle('Amplitude-O') + 
  theme(plot.title=element_text(hjust = 0.5, vjust = - 20))
  
amp.p3                  


# ELIXA 

kmfit_elixa <- survfit(Surv(time, status) ~ treat, data = elixa)

e_stpm2 <- rstpm2::stpm2(Surv(time, status) ~ treat, 
                     df = 3,tvc = list(treat = 1),  data = elixa)

elixa_fit <- predict(e_stpm2, newdata = data.frame(treat = c(0:1)),
                   grid = T, full = T, se.fit = T,
                   type = "surv")

elixa_fit$cif <- 1 - elixa_fit$Estimate

e.p <- jskm::jskm(sfit = kmfit_elixa,
                    cumhaz = T,legend = F,
                    linecols = "black",marks = F)

e.p

e.p2 = e.p + ylim(0,0.25)

e.p2

# convert amp to CIF 

e.p3 = e.p2 + geom_line(data = elixa_fit, aes(x = time, y = cif, 
                                                color = factor(treat)),
                            size = 1) +
  labs(x = "Months since Randomisation", 
       y = "Cumulative Incidence ") + ggtitle('ELIXA') + 
  theme(plot.title=element_text(hjust = 0.5, vjust = - 20))

e.p3                  

# SUSTAIN

kmfit_sustain <- survfit(Surv(time, status) ~ treat, 
                         data = sustain)

sustain_stpm2 <- rstpm2::stpm2(Surv(time, status) ~ treat, 
                         df = 3,tvc = list(treat = 1),  data = sustain)

sustain_fit <- predict(sustain_stpm2, newdata = data.frame(treat = c(0:1)),
                     grid = T, full = T, se.fit = T,
                     type = "surv")

sustain_fit$cif <- 1 - sustain_fit$Estimate

sustain.p <- jskm::jskm(sfit = kmfit_sustain,
                  cumhaz = T,legend = F,
                  linecols = "black",marks = F)

sustain.p

sustain.p2 = sustain.p + ylim(0,0.25)

sustain.p2

# convert amp to CIF 

sustain.p3 = sustain.p2 + geom_line(data = sustain_fit, aes(x = time, y = cif, 
                                              color = factor(treat)),
                        size = 1) +
  labs(x = "Months since Randomisation", 
       y = "Cumulative Incidence ") + ggtitle('Sustain 6') + 
  theme(plot.title=element_text(hjust = 0.5, vjust = - 20))

sustain.p3                 


# Harmony Outcomes 

kmfit_harmony <- survfit(Surv(time, status) ~ treat, 
                         data = harmony)

harmony_stpm2 <- rstpm2::stpm2(Surv(time, status) ~ treat, 
              df = 3,tvc = list(treat = 1),  data = harmony)

harmony_fit <- predict(harmony_stpm2, newdata = data.frame(treat = c(0:1)),
                       grid = T, full = T, se.fit = T,
                       type = "surv")

harmony_fit$cif <- 1 - harmony_fit$Estimate

harmony.p <- jskm::jskm(sfit = kmfit_harmony,
                        cumhaz = T,legend = F,
                        linecols = "black",marks = F)

harmony.p

harmony.p2 = harmony.p + ylim(0,0.25)

harmony.p2

# convert amp to CIF 

harmony.p3 = harmony.p2 + 
  geom_line(data = harmony_fit, aes(x = time, y = cif, 
                                                            color = factor(treat)),
                                    size = 1) +
  labs(x = "Months since Randomisation", 
       y = "Cumulative Incidence ") + ggtitle('Harmony \n Outcomes') + 
  theme(plot.title=element_text(hjust = 0.5, vjust = - 20))

harmony.p3   

library(patchwork)

plot <- amp.p3 + e.p3

plot

plot[[2]] <- plot[[2]] + theme(axis.title.y = element_blank())

plot

plot2 <- sustain.p3 + harmony.p3

plot2[[2]] <- plot2[[2]] + theme(axis.title.y = element_blank())

plot2

plot3 <- plot/plot2

plot3

ggsave(
  plot = plot3,
  filename = "F:\\GLP1_agonists\\analysis\\results\\RPfit_4trials.pdf",
  device = "pdf",
  height = 5,
  width = 8,
  units = "in"
)

# create plots for next figure panel 

# EXCSEL

kmfit_excsel <- survfit(Surv(time, status) ~ treat, 
                         data = excsel)

excsel_stpm2 <- rstpm2::stpm2(Surv(time, status) ~ treat, 
                df = 3,tvc = list(treat = 1),  data = excsel)

excsel_fit <- predict(excsel_stpm2, newdata = data.frame(treat = c(0:1)),
                       grid = T, full = T, se.fit = T,
                       type = "surv")

excsel_fit$cif <- 1 - excsel_fit$Estimate

excsel.p <- jskm::jskm(sfit = kmfit_excsel,
                        cumhaz = T,legend = F,
                        linecols = "black",marks = F)

excsel.p

excsel.p2 = excsel.p + ylim(0,0.25) 

excsel.p2

# convert to CIF 

excsel.p3 = excsel.p2 + 
  geom_line(data = excsel_fit, aes(x = time, y = cif, 
                                    color = factor(treat)),
            size = 1) +
  labs(x = "Months since Randomisation", 
       y = "Cumulative Incidence ") + ggtitle('EXCSEL') + 
  theme(plot.title=element_text(hjust = 0.5, vjust = - 20))

excsel.p3   

# LEADER

kmfit_leader <- survfit(Surv(time, status) ~ treat, 
                        data = leader)

leader_stpm2 <- rstpm2::stpm2(Surv(time, status) ~ treat, 
                df = 3,tvc = list(treat = 1),  data = leader)

leader_fit <- predict(leader_stpm2, newdata = data.frame(treat = c(0:1)),
                      grid = T, full = T, se.fit = T,
                      type = "surv")

leader_fit$cif <- 1 - leader_fit$Estimate

leader.p <- jskm::jskm(sfit = kmfit_leader,
                       cumhaz = T,legend = F,
                       linecols = "black",marks = F)

leader.p

leader.p2 = leader.p + ylim(0,0.25) 

leader.p2

# convert to CIF 

leader.p3 = leader.p2 + 
  geom_line(data = leader_fit, aes(x = time, y = cif, 
                                   color = factor(treat)),
            size = 1) +
  labs(x = "Months since Randomisation", 
       y = "Cumulative Incidence ") + ggtitle('LEADER') + 
  theme(plot.title=element_text(hjust = 0.5, vjust = - 20))

leader.p3 

# REWIND

kmfit_rewind <- survfit(Surv(time, status) ~ treat, 
                        data = rewind)

rewind_stpm2 <- rstpm2::stpm2(Surv(time, status) ~ treat, 
                              df = 3,tvc = list(treat = 1),  
                              data = rewind)

rewind_fit <- predict(rewind_stpm2, newdata = data.frame(treat = c(0:1)),
                      grid = T, full = T, se.fit = T,
                      type = "surv")

rewind_fit$cif <- 1 - rewind_fit$Estimate

rewind.p <- jskm::jskm(sfit = kmfit_rewind,
                       cumhaz = T,legend = F,
                       linecols = "black",marks = F)

rewind.p

rewind.p2 = rewind.p + ylim(0,0.25) 

rewind.p2

# convert to CIF 

rewind.p3 = rewind.p2 + 
  geom_line(data = rewind_fit, aes(x = time, y = cif, 
                                   color = factor(treat)),
            size = 1) +
  labs(x = "Months since Randomisation", 
       y = "Cumulative Incidence ") + ggtitle('REWIND') + 
  theme(plot.title=element_text(hjust = 0.5, vjust = - 20))

rewind.p3 

# Pioneer 6 

kmfit_pioneer <- survfit(Surv(time, status) ~ treat, 
                        data = pioneer)

pioneer_stpm2 <- rstpm2::stpm2(Surv(time, status) ~ treat, 
                              df = 3,tvc = list(treat = 1),  
                              data = pioneer)

pioneer_fit <- predict(pioneer_stpm2, newdata = data.frame(treat = c(0:1)),
                      grid = T, full = T, se.fit = T,
                      type = "surv")

pioneer_fit$cif <- 1 - pioneer_fit$Estimate

pioneer.p <- jskm::jskm(sfit = kmfit_pioneer,
                       cumhaz = T,legend = F,
                       linecols = "black",marks = F)

pioneer.p

pioneer.p2 = pioneer.p + ylim(0,0.10) 

pioneer.p2

# convert to CIF 

pioneer.p3 = pioneer.p2 + 
  geom_line(data = pioneer_fit, aes(x = time, y = cif, 
                                   color = factor(treat)),
            size = 1) +
  labs(x = "Months since Randomisation", 
       y = "Cumulative Incidence ") + ggtitle('Pioneer 6') + 
  theme(plot.title=element_text(hjust = 0.5, vjust = - 20))

pioneer.p3 

# combine the 4 plots to present panel plot

plot <- excsel.p3 + leader.p3

plot

plot[[2]] <- plot[[2]] + theme(axis.title.y = element_blank())

plot

plot2 <- rewind.p3 + pioneer.p3

plot2[[2]] <- plot2[[2]] + theme(axis.title.y = element_blank())

plot2

plot3 <- plot/plot2

plot3

ggsave(
  plot = plot3,
  filename = "F:\\GLP1_agonists\\analysis\\results\\RPfit_4trials_panel2.pdf",
  device = "pdf",
  height = 5,
  width = 8,
  units = "in"
)