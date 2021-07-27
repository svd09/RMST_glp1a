# script to obtain rmst and then pool for 3 point MACE outcome
# primary outcome for all the studies


library(easypackages)

libraries(c("tidyverse","IPDfromKM",'survival',
            "flexsurv","broom","rstpm2","survminer",
            "metaRMST", "ggthemes","ckbplotr"))


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

# run the model metaRMST using df and time @ 12, 24, 36, 48 and 
# fit the univariate model using integration of the KM curves.

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



# run the model metaRMST using df and time @ 12,24,36,48 fitting
# the flexible RP model and with extrapolation of curves.


res <- metaRMSTD(df2,time_horizons = c(12,24,36,48),
                 MA_method ="uni_flex")

df_res <- tbl_df(res$result)

result_flex <- df_res %>% select(time_horizon,Estimate,lower, upper, pval)

result_flex

result_flex$method <- 'parametric_model'
result_km$method <- 'KM integration'


result_comb <- rbind(result_flex, result_km)


# create a plot to compare both estimates.
# save this plot and going to include in the supplemental section.

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
  scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.2))
  
t3
  
ggsave(plot = t3,
       filename = "F:/GLP1_agonists/analysis/results/compare_models.pdf",
       height = 5,
       width = 8,
       units = "in")

  
  

### 


plot <- RMSTcurves(df, time_horizons = c(12, 24, 36, 48),
           MA_mvma = F,MA_mvma_boot = F,MA_uni = T,MA_uni_flex = T)


pdf("F:\\GLP1_agonists\\analysis\\results\\rmstplot.pdf",
    height = 5,
    width = 8)


RMSTplot(plot,
         ylim = c(-0.1, 2),
         trial_legend = F,
         estimates = T,
         MA_legend = F,
         xlim = c(0,48),
         ylab = "RMST Difference (Months)",
         xlab = "Months since Randomisation",
         col = c("red","blue","green","orange",
                 "purple","yellow","brown","gray"))
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
       cex = 0.8)


dev.off()

# only km method...

plot_km <- RMSTcurves(df, time_horizons = c(12, 24, 36, 48),
        MA_mvma = F,MA_mvma_boot = F,MA_uni = T,MA_uni_flex = F)



RMSTplot(plot_km,
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
       cex = 0.8)



# parametric method.


plot_rp <- RMSTcurves(df, time_horizons = c(12, 24, 36, 48),
                      MA_mvma = F,MA_mvma_boot = F,MA_uni = F,MA_uni_flex = T)



RMSTplot(plot_rp,
         # ylim = c(-0.1, 1),
         trial_legend = F,
         estimates = T,
         MA_legend = F,
         xlim = c(0,48),
         ylab = "RMST Difference (Months)",
         xlab = "Months since Randomisation",)
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
       cex = 0.8)




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

time <- seq(from = 0, to = 48, by = 3)

cont_res <- metaRMSTD(df2, 
                 time_horizons = time,
                 MA_method = "uni")

df_cont_res <- tbl_df(cont_res$result)

cont_result <- df_cont_res %>% select(time_horizon,Estimate,lower, upper, pval)

cont_result2 <- cont_result %>% drop_na()

glimpse(cont_result2)


# graph for whole time period.
# univariate model using KM curves.

p <- ggplot(data = cont_result2, aes(x = time_horizon, y = Estimate)) + 
  geom_line() + 
  geom_line(data = cont_result2, aes(x = time_horizon, y = lower), color = "red",
            linetype  = 2) +
  geom_line(data = cont_result2, aes(x = time_horizon, y = upper), color = "red",
            linetype = 2)
            
p2 <- p + theme_minimal()

p3 <- p2 + scale_x_continuous(breaks = c(12,24,36,48)) + 
  ylab("RMST Difference (Months)") + xlab("Months since Randomisation")


# Wald test to compare results for sensitivity analysis.

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

table$zscore =  ((table$Estimate*table$Estimate) - 
        (table$Estimate2*table$Estimate2))/
  sqrt((table$SE*table$SE) + (table$SE2*table$SE2))


table$Wald_test <- 2*pnorm(-abs(table$zscore))

# create a function so that can run that for the other values.
# this function can be used for the other tests.

wald_test_MA <- function(Estimate, SE, Estimate2, SE2){
  
  a = (Estimate*Estimate) - (Estimate2*Estimate2)
  
  b = sqrt((SE*SE) + (SE2*SE2))
  
  z = a/b
  
  p = 2*pnorm(-abs(z))
  
  return(p)
}

# am going to plot forest plot for the 48 months estimate for this data.

dfm <- res$rmstd_est %>% tbl_df()

dfse <- res$se_rmstd_est %>% tbl_df()


glimpse(dfm)

glimpse(dfse)

table <- tibble(
  Study = c("Amplitude-O","ELIXA","Sustain6", 
            "Harmony Outcomes",'EXCEL',"LEADER","REWIND","PIONEER6",
            "Pooled Estimate"),
  estimate = c(dfm$RMSTD_est_at_48,0.627),
  stderr= c(dfse$se_RMSTD_est_at_48, 0.182)
)

table

table$shape = rep(15, 9)

table$color = rep("black",9)

table$fill = rep("black",9)

table$variable = c('a','b','c','d','e','f','g','h','i')

table <- data.frame(table)

table$row.level.labels = rep("heading",9)



# now to create the forest plot 

rowlabels <- data.frame(heading = rep('Study',9),
                     variable = c('a','b','c','d','e','f','g','h','i'),
                     row.labels.levels = rep("heading",9))


forestplot1 <- 
  make_forest_plot(panels = list(table),
  col.key          = "Study",
  row.labels       = ckbplotr_row_labels,
#  row.labels.levels = c("heading", "subheading", "label"),
  rows             = c("Diff RMST"),
  exponentiate     = FALSE,
  panel.names      = c("RMST"),
  blankrows        = c(0, 1, 0, 1),
  scalepoints      = TRUE,
  pointsize        = 3,
  shape            = "shape",
  colour           = "colour",
#  col.bold         = "bold",
#  col.diamond      = "diamond",
#  ciunder          = "ciunder"
)


fp1 <- make_forest_plot(
  panels = list(table),
  col.key = "variable",
  row.labels = rowlabels,
#  col.estimate = estimate,
# col.stderr = stderr,
  row.labels.levels = c("heading"),
  rows = "Study",
  exponentiate = F
  
)




