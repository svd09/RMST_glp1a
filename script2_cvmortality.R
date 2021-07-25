# script to pool cv mortality data from all trials 

# get libraries

library(tidyverse)
library(metaRMST)
library(ggthemes)
library(survival)

# get the data


# get all the data.

sustain6 <- 
  read_csv('F:/GLP1_agonists/data/pooled_data/cvmort_sustain6.csv')

sustain6$trialID <- 1

harmony <- 
  read_csv('F:/GLP1_agonists/data/pooled_data/cvmort_harmony.csv')

harmony$trialID <- 2

excsel <- 
  read_csv('F:/GLP1_agonists/data/pooled_data/cvmort_excsel.csv')

excsel$trialID <- 3

leader <- 
  read_csv('F:/GLP1_agonists/data/pooled_data/cvmort_leader.csv')

leader$trialID <- 4

rewind <- 
  read_csv('F:/GLP1_agonists/data/pooled_data/cvmort_rewind.csv')

rewind$trialID <- 5

pioneer6 <- 
  read_csv('F:/GLP1_agonists/data/pooled_data/cvmort_pioneer6.csv')

pioneer6$trialID <- 6


df <- rbind(pioneer6, rewind, 
            leader, excsel, 
            sustain6, harmony)


df$triaID <- df$trialID
df$Time <- df$time
df$Event <- df$status
df$Arm <- df$treat

df2 <- df %>% select(trialID, Time, Event, Arm)

# do a univariate KM model to see trial data.



res_uni <- metaRMSTD(df2, 
                     time_horizons = c(12,24,36,48),
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

# flexible parametric model 

res_flex <- metaRMSTD(df2,
                 time_horizons = c(12,24,36,48),
                 MA_method ="uni_flex")

df_res <- tbl_df(res_flex$result)

result_flex <- df_res %>% select(time_horizon,Estimate,lower, upper, pval)

result_flex

# combine both results and then plot.


result_flex$method <- 'parametric_model'
result_km$method <- 'KM integration'


result_comb <- rbind(result_flex, result_km)


t <- ggplot(data = result_comb) + 
  geom_point(aes(x = time_horizon,
                 y = Estimate,
                 group = method,
                 color = method),
             size = 2, position = position_dodge2(width = 3),
             shape = "square") 

t2 = t +
  geom_linerange(aes( x = time_horizon,
                      ymin = lower,
                      ymax = upper,
                      group = method,
                      color = method),
                 position = position_dodge2(width = 3))

t2


t3 <- t2 + ylab("RMST Difference (Months)") + 
  xlab("Months since Randomisation") + 
  theme_minimal() + scale_x_continuous(breaks = c(12,24,36,48)) + 
  scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.2)) + 
  geom_hline(yintercept = 0, linetype = 3)

t3

# analysis removing Pioneer6 - 

df_exc_pio <- df2 %>% filter(trialID != c(1,6))



res_uni_exc_pio <- metaRMSTD(df_exc_pio, 
                     time_horizons = c(12,24,36,48),
                     MA_method ="uni")

res_uni_exc_pio

result_univ_exc_pio_km <- tbl_df(res_uni_exc_pio$result)

result_univ_exc_pio_km

result_exc_pio_km <- result_univ_exc_pio_km %>% select(time_horizon,
                                       
                                       Estimate,
                                       upper,
                                       lower,
                                       pval)

result_exc_pio_km

# Wald test for sub-group comparison
# compare parametric method and KM method using the Wald test 

source('wald_test_MA.R')


km_w <- result_univ_km %>% select(time_horizon,
                                  Estimate,
                                  SE)


flex_w <- df_res %>% select(time_horizon,
                            Estimate,
                            SE)


flex_w2 <- flex_w %>% rename(
  Estimate2 = Estimate,
  SE2 = SE
) 

table <- left_join(km_w, flex_w2, by = "time_horizon")

table

table2 <- table %>% filter(time_horizon %in% c(24,48))

table2$wald_test_p <- with(table2, 
                           wald_test_MA(Estimate = 
                                          Estimate, 
                                        SE = SE, 
                                        Estimate2 = Estimate2, 
                                        SE2 = SE2))


table2