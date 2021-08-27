#' KMSubtractionEvaluateMatch
#'
#' This function evaluates of the quality of matching by inspecting Empirical cumulative distribution functions and the Kolmogorov−Smirnov Test of follow-up times between matched pairs as well as with Bland-Altman plots to explore discrepancies between matched pairs.
#' @param df_match A dataframe (data) from KMSubtractionMatch() containing time-to-event data of the matched cohorts "time" and "status"; the curve it originates from "curve"; as well as the "matchid".
#' @keywords
#' @export
#' @examples
#' data(cancer)
#' df_overall=colon
#' df_subgroup=colon[1:200,]
#' match=KMSubtractionMatch(df_overall, df_subgroup, matching="bipartite")
#'
#' KMSubtractionEvaluateMatch(df_match=match$data)

KMSubtractionEvaluateMatch=function(df_match){
  # events
  df_match_overall=subset(df_match, df_match$status==1 & !is.na( df_match$matchid) & df_match$curve=="overall")
  df_match_overall=df_match_overall[order(df_match_overall$matchid),]
  df_match_subgroup=subset(df_match, df_match$status==1 & !is.na(df_match$matchid) & df_match$curve=="subgroup")
  df_match_subgroup=df_match_subgroup[order(df_match_subgroup$matchid),]

  # bland altman plot, events
  df_matched_timediff=cbind(matchid=df_match_overall$matchid,
                            overall=df_match_overall$time,
                            subgroup=df_match_subgroup$time,
                            time_diff_months=(df_match_overall$time-df_match_subgroup$time)) %>% data.frame

  blandr_events=blandr.statistics(df_matched_timediff$overall, df_matched_timediff$subgroup , sig.level=0.95)

  plot_comparison_blandr_events=blandr.draw(df_matched_timediff$overall , df_matched_timediff$subgroup, ciDisplay = F , ciShading = T,   point_size = 0.1,
                                            plotTitle = paste("Bland-Altman plot\nDifferences in follow-up time\nbetween matched pairs with events"), y="Time difference") +
    # ylim(-5,5)+
    annotate("text", y=3, x=0,  hjust = 0,label = paste("Total follow-up time (Overall) = " , format(round(sum(df_matched_timediff$overall),3), nsmall=3),
                                                        "\nTotal follow-up time (Subgroup) = ", format(round(sum(df_matched_timediff$subgroup),3), nsmall=3),
                                                        "\nMean of absolute differences = ", format(round(mean(abs(blandr_events$differences)),5), nsmall=5),
                                                        sep=''), size = 3)

  # Empirical cumulative distribution function, events
  df_cum=cbind(fu=c(df_matched_timediff$overall,df_matched_timediff$subgroup), Cohort=c(rep("Overall (matched)", length(df_matched_timediff$overall)), rep(paste("Subgroup"), length(df_matched_timediff$subgroup)))) %>% data.frame()
  df_cum$fu=as.numeric(df_cum$fu)

  ks.test=ks.test(df_matched_timediff$overall, df_matched_timediff$subgroup)

  plot_cum_events=ggplot(df_cum, aes(fu, colour = Cohort)) +
    stat_ecdf()+
    theme(plot.title = element_text(hjust = 0.5))+
    scale_colour_manual(values=c("#000000", "#00b3aa"))+
    labs(title="Empirical cumulative distribution function\nbetween matched pairs with events" , x="Follow-up")+
    annotate("text", y=0.7, x=1,  hjust = 0,label = paste("Kolmogorov-Smirnov Test, p=", format(round(ks.test$p.value,3), nsmall=3),
                                                          sep=''), size = 3)+

    theme(legend.position="top")


  evaluate.events=data.frame(rbind(
                  'Kolmogorov-Smirnov Test, p'= format(round(ks.test$p.value,3), nsmall=3),
                  "Total follow-up time (Overall)" =  format(round(sum(df_matched_timediff$overall),3), nsmall=3),
                  "Total follow-up time (Subgroup)"= format(round(sum(df_matched_timediff$subgroup),3), nsmall=3),
                  "Mean of absolute differences" =format(round(mean(abs(blandr_events$differences)),5), nsmall=5)
                  ))


  # censor
  df_match_overall=subset(df_match, df_match$status==0 & !is.na( df_match$matchid) & df_match$curve=="overall")
  df_match_overall=df_match_overall[order(df_match_overall$matchid),]
  df_match_subgroup=subset(df_match, df_match$status==0 & !is.na(df_match$matchid) & df_match$curve=="subgroup")
  df_match_subgroup=df_match_subgroup[order(df_match_subgroup$matchid),]

  # bland altman plot, censor
  df_matched_timediff=cbind(matchid=df_match_overall$matchid,
                            overall=df_match_overall$time,
                            subgroup=df_match_subgroup$time,
                            time_diff_months=(df_match_overall$time- df_match_subgroup$time)) %>% data.frame

  blandr_events=blandr.statistics(df_matched_timediff$overall, df_matched_timediff$subgroup , sig.level=0.95)

  plot_comparison_blandr_censorships=blandr.draw(df_matched_timediff$overall , df_matched_timediff$subgroup, ciDisplay = F , ciShading = T,   point_size = 0.1,
                                                 plotTitle = paste("Bland-Altman plot\nDifferences in follow-up time\nbetween matched pairs with censorships"), y="Time difference") +
    # ylim(-5,5)+
    annotate("text", y=3, x=0,  hjust = 0,label = paste("Total follow-up time (Overall) = " , format(round(sum(df_matched_timediff$overall),3), nsmall=3),
                                                        "\nTotal follow-up time (Subgroup) = ", format(round(sum(df_matched_timediff$subgroup),3), nsmall=3),
                                                        "\nMean of absolute differences = ", format(round(mean(abs(blandr_events$differences)),5), nsmall=5),
                                                        sep=''), size = 3)

  # Empirical cumulative distribution function, events
  df_cum=cbind(fu=c(df_matched_timediff$overall,df_matched_timediff$subgroup), Cohort=c(rep("Overall (matched)", length(df_matched_timediff$overall)), rep(paste("Subgroup"), length(df_matched_timediff$subgroup)))) %>% data.frame()
  df_cum$fu=as.numeric(df_cum$fu)

  ks.test=ks.test(df_matched_timediff$overall, df_matched_timediff$subgroup)


  plot_cum_censorships=ggplot(df_cum, aes(fu, colour = Cohort)) +
    stat_ecdf()+
    theme(plot.title = element_text(hjust = 0.5))+
    scale_colour_manual(values=c("#000000", "#00b3aa"))+
    labs(title="Empirical cumulative distribution function\nbetween matched pairs with censorships" , x="Follow-up")+
    annotate("text", y=0.7, x=1,  hjust = 0,label = paste("Kolmogorov-Smirnov Test, p=", format(round(ks.test$p.value,3), nsmall=3),
                                                          sep=''), size = 3)+

    theme(legend.position="top")


  evaluate.censor=data.frame(rbind(
                  'Kolmogorov-Smirnov Test, p'= format(round(ks.test$p.value,3), nsmall=3),
                  "Total follow-up time (Overall)" =  format(round(sum(df_matched_timediff$overall),3), nsmall=3),
                  "Total follow-up time (Subgroup)"= format(round(sum(df_matched_timediff$subgroup),3), nsmall=3),
                  "Mean of absolute differences" =format(round(mean(abs(blandr_events$differences)),5), nsmall=5)
                  ))


  # KM plots
  km=survfit(Surv(time, status)~ curve, data=subset(df_match, !is.na( df_match$matchid)))
  plot_km=ggsurvplot(km,
             data = subset(df_match, !is.na( df_match$matchid)),
             size=0.7,
             risk.table = TRUE,
             censor.shape="|",
             censor.size = 1.2,
             conf.int = F,
             xlab = "Time",
             palette=c("#000000", "steelblue4"),
             ggtheme = theme(),
             risk.table.y.text = F,
             risk.table.fontsize=3,
             tables.theme = theme_cleantable(),
             title="Comparison of matched patients")


  plot_combined = ggarrange(ggarrange(plot_cum_events, plot_comparison_blandr_events, plot_cum_censorships, plot_comparison_blandr_censorships,
              nrow =2, ncol = 2),

              ggarrange(plot_km$plot, plot_km$table, nrow=2, ncol=1, heights=c(0.8,0.2)),

              nrow=1,ncol=2)

  # cleaning and formatting tables
  evaluate.events[,1]=as.numeric(evaluate.events[,1])
  evaluate.censor[,1]=as.numeric(evaluate.censor[,1])
  colnames(evaluate.events)=""
  colnames(evaluate.censor)=""

  # export
  out=list()
  out$'Evaluate match among patients with events'=evaluate.events
  out$'Evaluate match among patients with censorships'=evaluate.censor
  out$plot=plot_combined
  class(out)="KMSubtractionEvaluateMatch"
  out

}









