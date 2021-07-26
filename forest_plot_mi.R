# forest plot for myocardial infarction
# run the parametric model and then combine to plot.
# INCOMPLETE...


sustain <- 
  read_csv('F:/GLP1_agonists/data/pooled_data/mi_sustain6.csv')

sustain$trialID <- 1

harmony <- 
  read_csv('F:/GLP1_agonists/data/pooled_data/mi_harmony.csv')

harmony$trialID <- 2

leader <- 
  read_csv('F:/GLP1_agonists/data/pooled_data/mi_leader.csv')

leader$trialID <- 3

rewind <- 
  read_csv('F:/GLP1_agonists/data/pooled_data/mi_rewind.csv')

rewind$trialID <- 4



df <- rbind(rewind, leader, sustain, 
            harmony)


df$triaID <- df$trialID
df$Time <- df$time
df$Event <- df$status
df$Arm <- df$treat

df2 <- df %>% select(trialID, Time, Event, Arm)



