# now to pool all cause mortality for the study
# 2 trials - LEADER, EXCSEL 

library(survival)
library(tidyverse)
library(ggthemes)
library(metaRMST)


leader <- read_csv(
  'F:\\GLP1_agonists\\data\\pooled_data\\allmort_leader.csv'
)

excsel <- read_csv(
  'F:\\GLP1_agonists\\data\\pooled_data\\allmort_excsel.csv'
)

# combine trials and rename variables for model fit.

leader$trialID <- 1

excsel$trialID <- 2


df <- rbind( leader, excsel )


df$triaID <- df$trialID
df$Time <- df$time
df$Event <- df$status
df$Arm <- df$treat

df2 <- df %>% select(trialID, Time, Event, Arm)


# univariate model for 48 months 



res_uni <- metaRMSTD(df2, 
                     time_horizons = c(24,48),
                     MA_method ="uni")

res_uni

result_univ_km <- tbl_df(res_uni$result)

result_univ_km

result_km <- result_univ_km %>% select(time_horizon,
                   

                            
                                       Estimate,
                                       upper,
                                       lower,
                                       pval)

result_km

# univariate flexible model 

res_flex <- metaRMSTD(df2,
                      time_horizons = c(24,48),
                      MA_method ="uni_flex")

df_res <- tbl_df(res_flex$result)

result_flex <- df_res %>% select(time_horizon,Estimate,lower, upper, pval)

result_flex


# plot both 


result_flex$method <- 'parametric_model'
result_km$method <- 'KM integration'


result_comb <- rbind(result_flex, result_km)


t <- ggplot(data = result_comb) + 
  geom_point(aes(x = factor(time_horizon),
                 y = Estimate,
                 group = method,
                 color = method),
             size = 2, position = position_dodge2(width = 0.5),
             shape = "square") 

t2 = t +
  geom_linerange(aes( x = factor(time_horizon),
                      ymin = lower,
                      ymax = upper,
                      group = method,
                      color = method),
                 position = position_dodge2(width = 0.5))

t2


t3 <- t2 + ylab("RMST Difference (Months)") + 
  xlab("Months since Randomisation") + 
  theme_minimal() + 
  scale_y_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5)) + 
  geom_hline(yintercept = 0, linetype = 3)

t3

