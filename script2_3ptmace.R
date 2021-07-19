# script to obtain rmst and then pool for 3 point MACE outcome
# primary outcome for all the studies


library(easypackages)

libraries(c("tidyverse","IPDfromKM",'survival',
            "flexsurv","broom","rstpm2","survminer",
            "metaRMST", "ggthemes"))


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


df <- rbind(pioneer, rewind, leader, excsel, sustain, amplitude, 
            harmony, elixa)


df$triaID <- df$trialID
df$Time <- df$time
df$Event <- df$status
df$Arm <- df$treat

df2 <- df %>% select(trialID, Time, Event, Arm)

# run the model metaRMST using df and time @ 12,24,36,48 fitting
# the flexible RP model and with extrapolation of curves.

res <- metaRMSTD(df2, 
                 time_horizons = c(12,24,36,48),
                 MA_method ="uni_flex")

df_res <- tbl_df(res$result)

result <- df_res %>% select(time_horizon,Estimate,lower, upper, pval)

result

plot <- RMSTcurves(df, time_horizons = c(12, 24, 36, 48),
           MA_mvma = F,MA_mvma_boot = F,MA_uni = T,MA_uni_flex = T)

RMSTplot(plot,
         ylim = c(-0.1, 1),
         trial_legend = F,
         estimates = T,
         MA_legend = F,
         xlim = c(0,48),
         ylab = "RMST Difference (Months)",
         xlab = "Months since Randomisation")
legend(x = "topleft",
       legend = c(
         "AMPLITUDE-O",
         "ELIXA",
         "SUSTAIN-6",
         "Harmony Outcomes",
         "EXCSEL",
         "LEADER",
         "REWIND",
         "PIONEER-6"),
       fill = c("red","blue","green","orange",
               "purple","yellow","brown","gray"),
       
       bty = "n",
       cex = 0.8
       )

# create plot using ggplot

result$time_horizon <- factor(result$time_horizon)

a <- ggplot() + geom_pointrange(data = result, aes(x = time_horizon,
                               y = Estimate,
                               ymin = lower, 
                               ymax = upper))

a

a2 <- a + theme_minimal()

a3 <- a2 + xlab("Months Since Randomisation") + 
  ylab("RMST Difference(Months)") + 
  annotate(geom = "label",
           x = 1, y = 0.25, 
           label = "0.034(0.014 - 0.056)",
           fill = "lightblue",
           color = "black") + 
  
  annotate(geom = "label",
           x = 2, y = 0.50, 
           label = "0.158(0.083 - 0.233)",
           fill = "orange",
           color = "black") +

  annotate(geom = "label",
x = 3, y = 0.75, 
label = "0.368(0.178 - 0.558)",
fill = "lightgray",
color = "black") + 
  
  annotate(
    geom = "label",
    x = 4, y = 1.25, 
    label = "0.627(0.270 - 0.984)",
    fill = "#CC79A7",
    color = "black"
  )
  

ggsave(plot = a3, 
       filename = "F:/GLP1_agonists/analysis/results/pooled_rmstd.pdf",
       width = 8,
       height = 5,
       units = "in")

# code to create a continuous plot.

time <- seq(from = 0, to = 48, by = 0.1)

cont_res <- metaRMSTD(df2, 
                 time_horizons = time,
                 MA_method = "uni")

df_cont_res <- tbl_df(cont_res$result)

cont_result <- df_cont_res %>% select(time_horizon,Estimate,lower, upper, pval)

cont_result2 <- cont_result %>% drop_na()

glimpse(cont_result2)


### need to work on the graph
### add ribbon, color, axes, other details

ggplot(data = cont_result2, aes(x = time_horizon, y = Estimate)) + 
  geom_line() + 
  geom_ribbon(data = cont_result2, aes(ymin = lower, ymax = upper, fill = "gray")) +
  scale_alpha(0.2)

