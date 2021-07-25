# mi pooled results 

library(tidyverse)
library(metaRMST)
library(ggthemes)

# combine data and then model 



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

# univariate model at 24, 48 months

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

# univariate model with flexible parametric modeling.


res <- metaRMSTD(df2,time_horizons = c(24,48),
                 MA_method ="uni_flex")

df_res <- tbl_df(res$result)

result_flex <- df_res %>% select(time_horizon,Estimate,lower, upper, pval)

result_flex

# combine data and then plot.
# combine plot and then graph according to earlier method.


result_flex$Method <- 'Parametric Model'
result_km$Method <- 'Sensitvity Analysis (Kaplan Meier Method)'

result_comb <- rbind(result_flex, result_km)


t <- ggplot(data = result_comb) + 
  geom_point(aes(x = factor(time_horizon),
                 y = Estimate,
                 group = Method,
                 color = Method),
             size = 2, position = position_dodge2(width = 0.5),
             shape = "square") + 
        scale_fill_manual(values = c("blue","red")) 

t2 = t +
  geom_linerange(aes( x = factor(time_horizon),
                      ymin = lower,
                      ymax = upper,
                      group = Method,
                      color = Method),
                 position = position_dodge2(width = 0.5)) +
  scale_color_manual(values = c("blue","red")) 
t2


t3 <- t2 + ylab("RMST Difference (Months)") + 
  xlab("Months since Randomisation") + 
  theme_minimal() + 
  scale_y_continuous(breaks = c(0, 0.1, 0.2, 0.3, 
                                0.4, 0.5,0.6,0.7,
                                0.8,0.9,1)) + 
  geom_hline(yintercept = 0, linetype = 3) +
  
  theme(legend.position = "bottom",
        legend.title = element_text(face = "bold")) 


t3

ggsave(plot = t3, 
filename = "F:/GLP1_agonists/analysis/results/pooled_rmstd_mi.pdf",
width = 8,
              height = 5,
              units = "in")


# Wald test for subgroup analyses 

source("wald_test_MA.R")

d1 <- result_univ_km %>% select(Estimate, SE, time_horizon)

d2 <- df_res %>% select(time_horizon, Estimate, SE) %>%
  rename(
    Estimate2 = Estimate,
    SE2 = SE) 

df <- left_join(d1, d2, by = "time_horizon")

df$wald_test_p = with(df, wald_test_MA(Estimate = 
                                Estimate,
                                SE = SE,
                                Estimate2 = Estimate2,
                                SE2 = SE2
                              ))

df