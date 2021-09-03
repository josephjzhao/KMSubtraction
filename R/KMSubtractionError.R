#' Conduct Monte Carlo simulations to evaluate the limits of error of `KMSubtractionMatch`
#'
#' This function conducts Monte Carlo simulations to evaluate the limits of error of KMSubtractionMatch() given parameters surrounding the reconstruction task required.
#' Follow-up time was modeled by a random weibull distribution of common shape parameter of 1.000 and scale parameter of 5.000.
#'
#' @param n Size of overall cohort.
#' @param subgroup.p Proportion of reported subgroup.
#' @param censor_overall.p Proportion of patients with censorship status in the overall cohort, which may be found from the object of KMSubtractionMatch().
#' @param censor_subgroup.p Proportion of patients with censorship status in the subgroup cohort, which may be found from the object of KMSubtractionMatch().
#' @param interval Number of number-at-risk table intervals.
#' @param missing.p Proportion of missing data from the opposing subgroup. Default is set at 0.00.
#' @param ncores Number of cores to be utilized for parallel processing. The number of cores in your device may be found using the function 'parallel::detectCores()' from the parallel package.
#' @param mc Number of Monte Carlo iterations. Default is set at 1000.
#' @keywords
#' 
#' @return Reconstructed and original survival data were compared by means of marginal Cox-proportional hazard models and restricted mean survival time difference (RMSTD). This function returns density plots per matching algorithm for both ln(HR) and RMST-Difference; as well as summary statistics in table format. 
#'
#' @export
#' @examples
#' Size of dataset, censorship proportion and subgroup proportion may be retrieved from the KMSubtractionMatch object, under Parameters.
#' 
#' data(cancer)
#' df_overall=colon
#' df_subgroup=colon[1:200,]
#' match=KMSubtractionMatch(df_overall, df_subgroup, matching="bipartite")
#' 
#' match$Parameters
#' 
#' n=match$Parameters[1,1]
#' subgroup.p=match$Parameters[3,1]
#' censor_overall.p=match$Parameters[4,1]
#' censor_subgroup.p=match$Parameters[5,1]
#' 
#' KMSubtractionError(n=n,
#' mc=1000,
#' censor_overall.p=censor_overall.p,
#' censor_subgroup.p=censor_subgroup.p,
#' subgroup.p=subgroup.p,
#' interval=8,
#' missing.p=0.02,
#' ncores=5)

KMSubtractionError=function(n,
                           censor_overall.p,
                           censor_subgroup.p,
                           subgroup.p,
                           interval=8,
                           missing.p=0,
                           mc=1000,
                           ncores=1,
                           fig.label=""){

  #### Prepare parallel processing settings ####
  # count the total number of iterations
  iterations=length(n)*length(censor_overall.p)*length(censor_subgroup.p)*length(missing.p)*length(interval)*mc
  
  set.seed(123)
  
  # set number of cores
  registerDoParallel(ncores)
  registerDoSNOW(makeSOCKcluster(ncores))
  
  # progress bar
  pb <- txtProgressBar(min=1, max=iterations, style=3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=progress)
  
  library(parallel)
  library(doParallel)
  library(foreach)
  library(doSNOW)
  
  #### Simulation ####
  df_analysis=NULL
  system.time({
    df_analysis=
      foreach (i.n = 1:length(n), .options.snow=opts, .combine=rbind) %:%
        foreach (i.censor_overall = 1:length(censor_overall.p), .options.snow=opts, .combine=rbind) %:%
          foreach (i.censor_subgroup = 1:length(censor_subgroup.p), .options.snow=opts, .combine=rbind) %:% 
            foreach (i.subgroup = 1:length(subgroup.p), .options.snow=opts, .combine=rbind) %:%
              foreach (i.missing = 1:length(missing.p), .options.snow=opts, .combine=rbind) %:%
                foreach (i.interval = 1:length(interval), .options.snow=opts, .combine=rbind) %:%
                  foreach (i.mc = 1:mc, .options.snow=opts, .combine=rbind) %dopar% {
        # tryCatch({
          library(survival)
          library(survminer)
          library(wakefield)
          library(svglite)
          library(magick)
          library(stringi)
          library(stringr)
          library(dplyr)
          library(tidyr)
          library(tidyverse)
          library(scales)
          library(svglite)
          library(rsvg)
          library(IPDfromKM)
          library(readr)
          library(RcppHungarian)
          library(dplyr)
          library(MatchIt)
          library(survRM2)
          library(KMSubtraction)
  
  
          if ((subgroup.p[i.subgroup]*censor_overall.p[i.censor_overall])/censor_subgroup.p[i.censor_subgroup]>1){
            stop("Combination of parameters not possible!")
          }
  
          #### 1_Generate simulation data ####
          
          # number of patients per subgroup
          n.subgroup=round(subgroup.p[i.subgroup]*n[i.n],0)
          n.opposingsubgroup=round((n[i.n]-n.subgroup),0)
          
          # number of censored patients in the opposing subgroup
          n.censoroverall=round(censor_overall.p[i.censor_overall]*n[i.n], 0)
          n.censorsubgroup=round(censor_subgroup.p[i.censor_subgroup]*n.subgroup,0)
          n.censoroppsingsubgroup=n.censoroverall-n.censorsubgroup
          censor_opposingsubgroup.p=n.censoroppsingsubgroup/n.opposingsubgroup
          
          df_overall=NULL
          df_overall$time=rweibull(n.subgroup+n.opposingsubgroup, shape=1,scale=5)
          
          df_overall$status=c(rep(0, n.censorsubgroup),
                      rep(1, n.subgroup-n.censorsubgroup),
                      rep(0, n.censoroppsingsubgroup),
                      rep(1, n.opposingsubgroup-n.censoroppsingsubgroup))
          
          df_overall$subgroup=c(rep(1, n.subgroup),
                        rep(0, n.opposingsubgroup))
          
          df_overall$subgroup[df_overall$subgroup==0]=r_sample_binary(length(df_overall$subgroup[df_overall$subgroup==0]), x = c(NA,0), prob = c(missing.p[i.missing], 1-missing.p[i.missing]), name = "Binary")
          
          df_overall=data.frame(df_overall)
  
          group=paste0("_n",n[i.n],
                       "_subgroup", subgroup.p[i.subgroup],
                       "_censor_overall", censor_overall.p[i.censor_overall],
                       "_censor_subgroup", censor_subgroup.p[i.censor_subgroup],
                       "_missing", missing.p[i.missing],
                       "_interval", interval[i.interval],
                       "_mc", i.mc)
  
          df_subgroup=subset(df_overall, df_overall$subgroup==1)
          
          km_overall=survfit(Surv(time, status) ~ 1, data=df_overall)
          km_subgroup=survfit(Surv(time, status) ~ 1, data=df_subgroup)
          
          df_overall_clicks=getxycoordinates(df_overall, label=paste0("overall",group))
          df_subgroup_clicks=getxycoordinates(df_subgroup, label=paste0("subgroup",group))
  
          #### 2_Reconstruction ####
          # at risk tables
          df_overall_risktable=ggsurvtable(km_overall, data=df_overall, break.time.by = (max(df_overall$time)/interval[i.interval]))$risk.table$data[,2:3]
          df_subgroup_risktable=ggsurvtable(km_subgroup, data=df_subgroup, break.time.by = (max(df_subgroup$time)/interval[i.interval]))$risk.table$data[,2:3]
          
          # overall
          df_overall_recon=getIPD(prep=preprocess(dat=df_overall_clicks, 
                                                  trisk=df_overall_risktable$time, 
                                                  nrisk=df_overall_risktable$n.risk, 
                                                  maxy=1), 
                                                  armID=1, tot.events=NULL)$IPD
          
          # subgroup
          df_subgroup_recon=getIPD(prep=preprocess(dat=df_subgroup_clicks, 
                                                   trisk=df_subgroup_risktable$time, 
                                                   nrisk=df_subgroup_risktable$n.risk, 
                                                   maxy=1), 
                                                   armID=1, tot.events=NULL)$IPD
          
  
          #### 3_Matching ####
          match.algo=c("bipartite", "maha", "logit")
          
          tbl=NULL
          for (i.match in match.algo){
            
            df_match=KMSubtractionMatch(df_overall_recon, df_subgroup_recon, matching=i.match)$data
            
            df_match$strata="matched"
            df_overall$strata="original"
            
            df_combined=bind_rows(df_match, df_overall)
            
            # survival analysis
            # cox
            sum.cox=coxph(formula = Surv(time, status) ~ strata, data=subset(df_combined, df_combined$subgroup==0)) %>% summary
            GT=coxph(formula = Surv(time, status) ~ strata, data=subset(df_combined, df_combined$subgroup==0)) %>% cox.zph
            
            df_combined_0=subset(df_combined, df_combined$subgroup==0)
            # RMST
            rmst=rmst2(df_combined_0$time, df_combined_0$status, ifelse(df_combined_0$strata=="matched", 0,1), tau=NULL)
            
            # collate
            tbl=rbind(tbl,c(# identifiers
                            n=n[i.n],
                            subgroup.p=subgroup.p[i.subgroup],
                            censorship_overall.p=censor_overall.p[i.censor_overall],
                            censorship_subgroup.p=censor_subgroup.p[i.censor_subgroup],
                            missing.p=missing.p[i.missing],
                            interval=interval[i.interval],
                            mc=i.mc,
                            
                            # survival outcomes
                            logrank=sum.cox$sctest[3],
                            HR=sum.cox$coefficients[2],
                            TE=sum.cox$coefficients[1],
                            se=sum.cox$coefficients[3],
                            GT.p=GT$table[1,3],
                            RMSTD=rmst$unadjusted.result[1,1],
                            
                            # matching algorithm
                            matching=i.match))
            
          }
  
          #### 4_Export summary statistics ####
          tbl
  
        # },
        # error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
      }
  })
  
  # clear parallel processing
  stopImplicitCluster()
  close(pb)
  stopCluster(makeSOCKcluster(ncores))
  
  #### Final df ####
  df_analysis=data.frame(df_analysis)
  df_analysis[,-grep("matching|outcome", colnames(df_analysis))]=apply(df_analysis[,-grep("matching|outcome", colnames(df_analysis))],2, as.numeric)
  
  #### Extract mean sd ####
  table_cox=NULL
  for (i in 1:length(unique(df_analysis$matching))){
    mean=subset(df_analysis, df_analysis$matching==unique(df_analysis$matching)[i])$TE %>% mean
    sd=subset(df_analysis, df_analysis$matching==unique(df_analysis$matching)[i])$TE %>% sd
    
    abs.mean=subset(df_analysis, df_analysis$matching==unique(df_analysis$matching)[i])$TE %>% abs %>% mean
    abs.sd=subset(df_analysis, df_analysis$matching==unique(df_analysis$matching)[i])$TE %>% abs %>% sd
  
    table_cox=rbind(table_cox,c(
      Matching=unique(df_analysis$matching)[i],
      mean,
      sd=sd,
      absTE=abs.mean,
      abs.sd=abs.sd
      ))
  
  }
  
  table_cox=data.frame(table_cox)
  table_cox[,-grep("Matching", colnames(table_cox))]=apply(table_cox[,-grep("Matching", colnames(table_cox))],2, as.numeric)
  
  colnames(table_cox)=c("Matching", "ln(HR), mean", "ln(HR), sd", "|ln(HR)|, mean", "|ln(HR)|, sd")
  
  table_rmst=NULL
  for (i in 1:length(unique(df_analysis$matching))){
    mean=subset(df_analysis, df_analysis$matching==unique(df_analysis$matching)[i])$RMSTD %>% mean
    sd=subset(df_analysis, df_analysis$matching==unique(df_analysis$matching)[i])$RMSTD %>% sd
    
    abs.mean=subset(df_analysis, df_analysis$matching==unique(df_analysis$matching)[i])$RMSTD %>% abs %>% mean
    abs.sd=subset(df_analysis, df_analysis$matching==unique(df_analysis$matching)[i])$RMSTD %>% abs %>% sd
  
    table_rmst=rbind(table_rmst,c(
      Matching=unique(df_analysis$matching)[i],
      mean,
      sd=sd,
      absRMSTD=abs.mean,
      abs.sd=abs.sd
      ))
  
  }
  
  table_rmst=data.frame(table_rmst)
  table_rmst[,-grep("Matching", colnames(table_rmst))]=apply(table_rmst[,-grep("Matching", colnames(table_rmst))],2, as.numeric)
  colnames(table_rmst)=c("Matching", "RMST-D, mean", "RMST-D, sd", "|RMST-D|, mean", "|RMST-D|, sd")
  
  #### Histograms ####
  # Mean error
  plot_cox=ggplot(df_analysis, aes(x=TE, fill=matching))+
    geom_density(alpha=0.7)+
    # geom_histogram(aes(y=..density..), position="identity", alpha=0.3)+
    scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
    geom_vline(xintercept = 0, linetype="dashed", col="red")+
    labs(title=paste("KMSubtraction vs original (marginal cox model)",fig.label),x="ln(HR)", fill="Matching algorithm")+
    geom_vline(data=table_cox, aes(xintercept=`ln(HR), mean`), linetype="dashed")
  
  plot_rmst=ggplot(df_analysis, aes(x=RMSTD, fill=matching))+
    geom_density(alpha=0.7)+
    # geom_histogram(aes(y=..density..), position="identity", alpha=0.3)+
    scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
    geom_vline(xintercept = 0, linetype="dashed", col="red")+
    labs(title=paste("KMSubtraction vs original (RMST)",fig.label),x="RMST-Differences", fill="Matching algorithm")+
    geom_vline(data=table_rmst, aes(xintercept=`RMST-D, mean`), linetype="dashed")
  
  # Mean absolute error
  plot_abscox=ggplot(df_analysis, aes(x=abs(TE), fill=matching))+
    geom_density(alpha=0.7)+
    # geom_histogram(aes(y=..density..), position="identity", alpha=0.3)+
    scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
    geom_vline(xintercept = 0, linetype="dashed", col="red")+
    labs(title=paste("KMSubtraction vs original (marginal Cox model)",fig.label),x="|ln(HR)|", fill="Matching algorithm")+
    geom_vline(data=table_cox, aes(xintercept=`|ln(HR)|, mean`), linetype="dashed")
  
  plot_absrmst=ggplot(df_analysis, aes(x=abs(RMSTD), fill=matching))+
    geom_density(alpha=0.7)+
    # geom_histogram(aes(y=..density..), position="identity", alpha=0.3)+
    scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
    geom_vline(xintercept = 0, linetype="dashed", col="red")+
    labs(title=paste("KMSubtraction vs original (RMST)",fig.label),x="|RMST-Differences|", fill="Matching algorithm")+
    geom_vline(data=table_rmst, aes(xintercept=`|RMST-D|, mean`), linetype="dashed")
  
  #### Convergence diagnostics ####
  
  CI_z <- function (x, ci = 0.95){
    `%>%` <- magrittr::`%>%`
    standard_deviation <- sd(x)
    sample_size <- length(x)
    Margin_Error <- abs(qnorm((1-ci)/2))* standard_deviation/sqrt(sample_size)
    df_out <- c(mean=mean(x), 
                # sample_size=length(x), 
                # sd=sd(x),
                # Margin_Error=Margin_Error,
                u.ci=(mean(x) - Margin_Error),
                l.ci=(mean(x) + Margin_Error))
    return(df_out)
  }
  
  out=NULL
  for (i in unique(df_analysis$mc)){
    for (i.matching in unique(df_analysis$matching)){
      df_temp=subset(df_analysis, df_analysis$mc<=i & df_analysis$matching==i.matching)
      m=CI_z(df_temp$TE %>% abs)
      out=rbind(out, c(mc=i, abslnhr=m, matching=i.matching))
    }
  }
  
  out=data.frame(na.omit(out))
  out[,-grep("matching|outcome", colnames(out))]=apply(out[,-grep("matching|outcome", colnames(out))],2, as.numeric)
  
  plot_convergence=ggplot(out, color=matching,fill=matching)+
    geom_line(aes(x=mc, y=abslnhr.mean,color=matching), size=0.5)+
    geom_ribbon(aes(x=mc,y=abslnhr.mean, ymin=abslnhr.l.ci, ymax=abslnhr.u.ci,fill=matching), alpha=0.1)+
    scale_color_manual(values=c("red3", "steelblue3", "black"))+
    scale_fill_manual(values=c("black", "black", "black"))+
    labs(title=paste("Convergence",fig.label), y="|ln(HR)|, mean", x="Iteration", color="Matching algorithm")+
    guides(fill = "none")

  #### Export ####
  out=list()
  
  # plots
  out$density.cox=plot_cox
  out$density.rmst=plot_rmst
  out$density.abscox=plot_abscox
  out$density.absrmst=plot_absrmst
  out$convergence=plot_convergence
  
  # tables
  out$table.cox=table_cox
  out$table.rmst=table_rmst
  
  # df
  out$data=df_analysis
  
  class(out)="KMSubtractionError"
  
  out
  
}


