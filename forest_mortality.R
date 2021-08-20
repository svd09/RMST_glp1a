# forest plot for mortality. 
# combine plot for both CV and all cause mortality.
# 24, 48 months time frame.


library(easypackages)
libraries(c("metaRMST","patchwork","tidyverse",
            "ckbplotr"))

# first get the data set ready for CV mortality.


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

# run model for CV mortality.


res_flex <- metaRMSTD(df2,
                      time_horizons = c(24,48),
                      MA_method ="uni_flex")

df_res <- tbl_df(res_flex$result)

# now get the data to create the datasets.


dfm <- res_flex$rmstd_est %>% tbl_df()

dfse <- res_flex$se_rmstd_est %>% tbl_df()



table_24 <- tibble(
  study = c("Sustain6","Harmony Outcomes","EXCSEL",
            "LEADER","REWIND","Pioneer6", NA, "Pooled (RP)","Pooled (KM)"),
  estimate = c(dfm$RMSTD_est_at_24, NA,0.0409, 0.02),
  stderr = c(dfse$se_RMSTD_est_at_24, NA,0.0401, 0.02),
  variable = c("1","2","3","4","5","6",NA ,"Pooled (RP)","Pooled (KM)"),
  color = c("black","black","gray","black","gray","gray","NA","red","blue")
)

table_48 <- 
tibble(
  study = c("Sustain6","Harmony Outcomes","EXCSEL",
            "LEADER","REWIND","Pioneer6", NA, "Pooled (RP)", "Pooled (KM)"),
  estimate = c(dfm$RMSTD_est_at_48, NA, 0.163, 0.142),
  stderr = c(dfse$se_RMSTD_est_at_48, NA, 0.140, 0.153),
  variable = c("1","2","3","4","5","6",NA ,"Pooled (RP)", "Pooled (KM)"),
  color = c("gray","black","gray","black","gray","gray","NA","red","blue")
)



rowlabels = tibble(heading = c(rep("CV Mortality",9)),
                   subheading = c(rep("Trials Included",6),NA,NA,NA),
                   variable = c("1","2","3","4","5","6",NA, "Pooled (RP)", "Pooled (KM)"),
                   label = c("Sustain6","Harmony Outcomes","EXCSEL",
                                     "LEADER","REWIND","Pioneer6", NA, 
                             "Pooled Estimate (RP)", "Pooled (KM)"))




f <- make_forest_plot(panels = list(table_24, table_48),
                      
                      row.labels = rowlabels,
                      col.key = "variable",
                      rows = c("CV Mortality"),
                      row.labels.levels = c("heading","subheading",
                                            "label"),
                      panel.names = c("dRMST at 24 months","dRMST at 48 months"), 
                    
                      exponentiate = F,
                      col.right.heading = "RMST Diff.(95%CI)",
                      colour = "color",
                      nullval = 0,
                      xlim = c(-2,5),
                      xlab = "dRMST(95% CI)   favours GLP1 agonist")

f


# now am going to also run all cause mortality model and 
# then plot all together.



leader_m <- read_csv(
  'F:\\GLP1_agonists\\data\\pooled_data\\allmort_leader.csv'
)

excsel_m <- read_csv(
  'F:\\GLP1_agonists\\data\\pooled_data\\allmort_excsel.csv'
)

# combine trials and rename variables for model fit.

leader_m$trialID <- 1

excsel_m$trialID <- 2


df <- rbind( leader_m, excsel_m )


df$triaID <- df$trialID
df$Time <- df$time
df$Event <- df$status
df$Arm <- df$treat

df2_m <- df %>% select(trialID, Time, Event, Arm)

# run the model here.

# fixed model 

fixed <- metaRMSTD(df2_m,
                   time_horizons = c(24,48),
                   MA_method ="uni")
  



res_flex_m <- metaRMSTD(df2_m,
                      time_horizons = c(24,48),
                      MA_method ="uni_flex")


df_res_m <- tbl_df(res_flex_m$result)

df_m_est <- res_flex_m$rmstd_est %>% tbl_df()

df_m_se <- res_flex_m$se_rmstd_est %>% tbl_df()

df_m_est

df_m_se

table_24_m <- tibble(
  study = c("LEADER","EXCSEL", NA, "Pooled (PM)", "Pooled (KM)"),
  estimate = c(df_m_est$RMSTD_est_at_24, NA, 0.0553, 0.054),
  stderr = c(df_m_se$se_RMSTD_est_at_24, NA, 0.0265, 0.027),
  variable = c("1","2",NA,"Pooled (PM)", "Pooled (KM)"),
  color = c(rep("black",2),NA,"red", "blue")
)

table_48_m <- tibble(
  study = c("LEADER","EXCSEL", NA, "Pooled (PM)", "Pooled (KM)"),
  estimate = c(df_m_est$RMSTD_est_at_48, NA, 0.261, 0.265),
  stderr = c(df_m_se$se_RMSTD_est_at_48, NA, 0.0896, 0.090),
  variable = c("1","2",NA,"Pooled (PM)", "Pooled (KM)"),
  color = c(rep("black",2),NA,"red", "blue")
  
)


rowlabels_m = tibble(heading = c(rep("All cause Mortality",5)),
                   subheading = c(rep("Trials Included",2),NA,NA, NA),
                   variable = c("1","2",NA, "Pooled (PM)", "Pooled (KM)"),
                   label = c("LEADER", "EXCSEL", NA, "Pooled Estimate (PM)", "Pooled Estimate (KM)"))




f_m <- make_forest_plot(panels = list(table_24_m, table_48_m),
                      
                      row.labels = rowlabels_m,
                      col.key = "variable",
                      rows = c("All cause Mortality"),
                      row.labels.levels = c("heading","subheading",
                                            "label"),
                      panel.names = c("dRMST at 24 months","dRMST at 48 months"), 
                      
                      exponentiate = F,
                      col.right.heading = "RMST Diff.(95%CI)",
                      colour = "color",
                      nullval = 0,
                      xlim = c(-0.5,0.5),
                      xlab = "dRMST(95% CI)   favours GLP1 agonist")

f_m

#  AM GOING TO USE PATCHWORK TO COMBINE PLOTS 
# THIS IS THE FIGURE GOING INTO THE PAPER.

library(patchwork)


mort <- f$plot/f_m$plot

mort2 <- mort + plot_annotation(
  
  tag_levels = "A"
)


mort2



ggsave(mort2,
       filename = "F:\\GLP1_agonists\\analysis\\results\\mortality_comb.pdf",
       device = "pdf",
       height = 8,
       width = 8,
       units = "in")
















# see how to combine and plot 

t_24 <- rbind(table_24, table_24_m)

t_48 <- rbind(table_48, table_48_m)

rowlabels_c <- rbind(rowlabels, rowlabels_m)

rowlabels_c$variable2 <- c("1","2","3","4","5","6",NA,"8",
                           "9","10",NA,"12","13","14")

t_24$variable2 <- rowlabels_c$variable2

t_48$variable2 <- rowlabels_c$variable2

f_c <- make_forest_plot(panels = list(t_24, t_48),
                      
                      row.labels = rowlabels_c,
                      col.key = "variable2",
                      rows = c("CV Mortality", "All cause Mortality"),
                      row.labels.levels = c("heading","subheading",
                                            "label"),
                      panel.names = c("dRMST at 24 months","dRMST at 48 months"), 
                      
                      exponentiate = F,
                      col.right.heading = "RMST Diff.(95%CI)",
                      colour = "color",
                      nullval = 0,
                      xlim = c(-1.5,1.5),
                      xlab = "dRMST(95% CI)   favours GLP1 agonist")


ggsave(f_c$plot,
  filename = "F:\\GLP1_agonists\\analysis\\results\\mortality.pdf",
  device = "pdf",
  height = 5,
  width = 8,
  units = "in")