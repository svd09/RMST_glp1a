# create forest plot for 3pt MACE.
# created 4 panel forest plot for 3 pt MACE and then saved to folder.



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


# run the parametric model
# then combine data and then plot.



res <- metaRMSTD(df2,time_horizons = c(12,24,36,48),
                 MA_method ="uni_flex")

df_res <- tbl_df(res$result)

# get the data and then make a tibble


dfm <- res$rmstd_est %>% tbl_df()

dfse <- res$se_rmstd_est %>% tbl_df()



table_12 <- tibble(
  study = c("Amplitude-O","ELIXA","Sustain6", 
            "Harmony Outcomes",'EXCEL',"LEADER","REWIND","PIONEER6",NA,
            "Pooled Estimate (PM)", 'Pooled (KM)'),
  estimate = c(dfm$RMSTD_est_at_12,NA,0.0354, 0.0346),
  stderr = c(dfse$se_RMSTD_est_at_12,NA, 0.0109, 0.0114),
  variable = c("1","2","3","4","5","6","7","8",NA,"Pooled (PM)", "Pooled (KM)"),
  color = c(rep("black",8),NA,"red", "blue")
)


table_24 <- tibble(
  study = c("Amplitude-O","ELIXA","Sustain6", 
            "Harmony Outcomes",'EXCEL',"LEADER","REWIND","PIONEER6",NA,
            "Pooled Estimate (PM)", "Pooled (KM)"),
  estimate = c(dfm$RMSTD_est_at_24,NA,0.158, 0.154),
  stderr = c(dfse$se_RMSTD_est_at_24,NA, 0.0383, 0.037),
  variable = c("1","2","3","4","5","6","7","8",NA,"Pooled (PM)", "Pooled (KM)"),
  color = c("gray",rep("black",6),"gray",NA, "red", "blue")
)


table_36 <- tibble(
  study = c("Amplitude-O","ELIXA","Sustain6", 
            "Harmony Outcomes",'EXCEL',"LEADER","REWIND","PIONEER6",NA,
            "Pooled Estimate (PM)", "Pooled Estimate (KM)"),
  estimate = c(dfm$RMSTD_est_at_36,NA,0.368, 0.188),
  stderr = c(dfse$se_RMSTD_est_at_36,NA, 0.0971, 0.0919),
  variable = c("1","2","3","4","5","6","7","8",NA,"Pooled (PM)", "Pooled (KM)"),
  color = c("gray","black","gray","gray","black","black","gray","gray",NA,"red", "blue")
)


table_48 <- tibble(
  study = c("Amplitude-O","ELIXA","Sustain6", 
            "Harmony Outcomes",'EXCEL',"LEADER","REWIND","PIONEER6",NA,
            "Pooled Estimate (PM)", "Pooled Estimate (KM)"),
  estimate = c(dfm$RMSTD_est_at_48,NA,0.627, 0.381),
  stderr = c(dfse$se_RMSTD_est_at_48,NA, 0.182, 0.149),
  variable = c("1","2","3","4","5","6","7","8",NA,"Pooled (PM)", "Pooled (KM)"),
  color = c("gray","gray","gray","gray","black","black","gray","gray",NA,"red", "blue")
)

rowlabels = tibble(heading = c(rep("3 point MACE",11)),
                    subheading = c(rep("Trials Included",8),
                                   NA, NA, NA),
                    variable = c("1","2","3","4","5","6","7","8",NA, "Pooled (PM)", "Pooled (KM)"),
                    label = c("Amplitude-O","ELIXA","Sustain6", 
                              "Harmony Outcomes",'EXCEL',"LEADER","REWIND","PIONEER6",
                              NA, "Pooled Estimate (PM)", "Pooled Estimate (KM)"))





# now to plot the data...




f <- make_forest_plot(panels = list(table_12, table_24),
                      
                      row.labels = rowlabels,
                      col.key = "variable",
                      rows = c("3 point MACE"),
                      row.labels.levels = c("heading","subheading",
                                            "label"),
                      panel.names = c("dRMST at 12 months", "dRMST at 24 months"
                      ),
                      exponentiate = F,
                      col.right.heading = "RMST Diff.(95%CI)",
                      xlim = c(-1,3),
                      colour = "color",
                      nullval = 0,
                      xlab = " dRMST (95% CI)    favours GLP1agonists")


f



f2 <- make_forest_plot(panels = list(table_36, table_48),
                       
                       row.labels = rowlabels,
                       col.key = "variable",
                       rows = c("3 point MACE"),
                       row.labels.levels = c("heading","subheading",
                                             "label"),
                       panel.names = c("dRMST at 36 months", "dRMST at 48 months"
                       ),
                       exponentiate = F,
                       col.right.heading = "RMST Diff.(95%CI)",
                       colour = "color",
                       nullval = 0,
                       xlim = c(-1,3),
                       xlab = "dRMST (95% CI)   favours GLP1 agonists")

f2


# combine both plots using patchwork...
# both plots combined and saved as pdf.
# will need to change to tiff for the paper.

library(patchwork)

mace <- f$plot/f2$plot

mace2 <- mace + plot_annotation(
  title = "dRMST plot for 3 point MACE",
  subtitle = "Parametric Method",
  tag_levels = "A"
)

ggsave(mace2,
  filename = "F:\\GLP1_agonists\\analysis\\results\\3ptmace_same_axes.pdf",
  device = "pdf",
  height = 10,
  width = 10,
  units = "in")

# figure saved.