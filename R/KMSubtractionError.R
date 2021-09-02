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
  
  df_analysis=NULL
  library(parallel)
  library(doParallel)
  library(foreach)
  library(doSNOW)
  
  #### Simulation ####
  system.time({
    df_analysis=
      foreach (i.n = 1:length(n), .options.snow=opts, .combine=rbind) %:%
      foreach (i.censor_overall = 1:length(censor_overall.p), .options.snow=opts, .combine=rbind) %:%
      foreach (i.censor_subgroup = 1:length(censor_subgroup.p), .options.snow=opts, .combine=rbind) %:% 
      foreach (i.subgroup = 1:length(subgroup.p), .options.snow=opts, .combine=rbind) %:%
      foreach (i.missing = 1:length(missing.p), .options.snow=opts, .combine=rbind) %:%
      foreach (i.interval = 1:length(interval), .options.snow=opts, .combine=rbind) %:%
      foreach (i.mc = 1:mc, .options.snow=opts, .combine=rbind) %dopar% {
        tryCatch({
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
  
          # number of patients per subgroup
          n.subgroup=round(subgroup.p[i.subgroup]*n[i.n],0)
          n.opposingsubgroup=round((n[i.n]-n.subgroup),0)
          
          # number of censored patients in the opposing subgroup
          n.censoroverall=round(censor_overall.p[i.censor_overall]*n[i.n], 0)
          n.censorsubgroup=round(censor_subgroup.p[i.censor_subgroup]*n.subgroup,0)
          n.censoroppsingsubgroup=n.censoroverall-n.censorsubgroup
          censor_opposingsubgroup.p=n.censoroppsingsubgroup/n.opposingsubgroup
          
          # generating simulated data
          df=NULL
          df$time=rweibull(n.subgroup+n.opposingsubgroup, shape=1,scale=5)
          
          df$status=c(rep(0, n.censorsubgroup),
                      rep(1, n.subgroup-n.censorsubgroup),
                      rep(0, n.censoroppsingsubgroup),
                      rep(1, n.opposingsubgroup-n.censoroppsingsubgroup))
          
          df$subgroup=c(rep(1, n.subgroup),
                        rep(0, n.opposingsubgroup))
          
          df$subgroup[df$subgroup==0]=r_sample_binary(length(df$subgroup[df$subgroup==0]), x = c(NA,0), prob = c(missing.p[i.missing], 1-missing.p[i.missing]), name = "Binary")
          
          df=data.frame(df)
  
          group=paste0("_n",n[i.n],
                       "_subgroup", subgroup.p[i.subgroup],
                       "_censor_overall", censor_overall.p[i.censor_overall],
                       "_censor_subgroup", censor_subgroup.p[i.censor_subgroup],
                       "_missing", missing.p[i.missing],
                       "_interval", interval[i.interval],
                       "_mc", i.mc)
  
          # data
          df_overall=df
          df_subgroup=subset(df, df$subgroup==1)
  
          # plots
          scale_points=matrix(c(0,max(df$time),0,max(df$time),0,0,1,1), ncol=2, nrow=4)
  
          km_overall=survfit(Surv(time, status) ~ 1, data=df)
          km_subgroup=survfit(Surv(time, status) ~ 1, data=subset(df, df$subgroup==1))
  
          # as svg
          svglite(paste0("df_overall_curve",group,".svg"), width=5, height=5)
          plot(km_overall, col=c("#f042ed"), conf.int = F, lwd=1, xlim=c(0,max(df$time)))
          points(scale_points, col="#f04242", pch=16, cex=0.3)
          dev.off()
  
          svglite(paste0("df_subgroup_curve",group,".svg"), width=5, height=5)
          plot(km_subgroup, col=c("#f042ed"), conf.int = F, lwd=1, xlim=c(0,max(df$time)))
          points(scale_points, col="#f04242", pch=16, cex=0.3)
          dev.off()
  
  
          #### 2_Pixelization ####
          # overall
          im=image_read_svg(paste0("df_overall_curve",group,".svg"), width=1000)
          unlink(paste0("df_overall_curve",group,".svg"))
          df_im=image_raster(im,tidy=T) %>% filter(col=="#f04242ff" | col=="#f042edff")
          df_im$y=-df_im$y # flipping the y axis
  
          # rescale
          df_im$y=rescale(df_im$y, to=c(0,1))
          df_im$x=rescale(df_im$x, to=c(0,max(df$time)))
  
          # removing the scaling points
          df_scaled=subset(df_im, df_im$col=="#f042edff")
  
          # cleaning
          df_scaled$col=NULL
          rownames(df_scaled)=NULL
  
          df_overall_clicks=df_scaled
  
          # subgroup
          im=image_read_svg(paste0("df_subgroup_curve",group,".svg"), width=1000)
          unlink(paste0("df_subgroup_curve",group,".svg"))
          df_im=image_raster(im,tidy=T) %>% filter(col=="#f04242ff" | col=="#f042edff")
          df_im$y=-df_im$y # flipping the y axis
  
          # rescale
          df_im$y=rescale(df_im$y, to=c(0,1))
          df_im$x=rescale(df_im$x, to=c(0,max(df$time)))
  
          # removing the scaling points
          df_scaled=subset(df_im, df_im$col=="#f042edff")
  
          # cleaning
          df_scaled$col=NULL
          rownames(df_scaled)=NULL
  
          df_subgroup_clicks=df_scaled
  
          #### 3_Reconstruction ####
          # at risk tables
          df_overall_risktable=survminer::ggsurvtable(km_overall, data=df, break.time.by = (max(df$time)/interval[i.interval]))$risk.table$data[,2:3]
          df_subgroup_risktable=survminer::ggsurvtable(km_subgroup, data=df_subgroup, break.time.by = (max(df_subgroup$time)/interval[i.interval]))$risk.table$data[,2:3]
          
          # overall
          preprocess=preprocess(dat=df_overall_clicks, trisk=df_overall_risktable$time, nrisk=df_overall_risktable$n.risk, maxy=1)
          est_ipd=getIPD(prep=preprocess, armID=1, tot.events=NULL)
          df_overall_recon=est_ipd$IPD
  
          # subgroup
          preprocess=preprocess(dat=df_subgroup_clicks, trisk=df_subgroup_risktable$time, nrisk=df_subgroup_risktable$n.risk, maxy=1)
          est_ipd=getIPD(prep=preprocess, armID=1, tot.events=NULL)
          df_subgroup_recon=est_ipd$IPD
  
          #### 4_Matching ####
          df_match_BP=KMSubtractionMatch(df_overall_recon, df_subgroup_recon, matching="bipartite")$data
          df_match_Maha=KMSubtractionMatch(df_overall_recon, df_subgroup_recon, matching="maha")$data
          df_match_Logit=KMSubtractionMatch(df_overall_recon, df_subgroup_recon, matching="logit")$data
  
          # compare
          df_original=df_overall
  
          # identify which data set is which
          df_match_BP$strata="matched"
          df_match_Maha$strata="matched"
          df_match_Logit$strata="matched"
          df_original$strata="original"
  
          df_combined_BP=bind_rows(df_match_BP, df_original)
          df_combined_Maha=bind_rows(df_match_Maha, df_original)
          df_combined_Logit=bind_rows(df_match_Logit, df_original)
  
          # compare within cox model
          sum.cox_BP=coxph(formula = Surv(time, status) ~ strata, data=subset(df_combined_BP, df_combined_BP$subgroup==0)) %>% summary
          sum.cox_Maha=coxph(formula = Surv(time, status) ~ strata, data=subset(df_combined_Maha, df_combined_Maha$subgroup==0)) %>% summary
          sum.cox_Logit=coxph(formula = Surv(time, status) ~ strata, data=subset(df_combined_Logit, df_combined_Logit$subgroup==0)) %>% summary
  
  
          # GT test
          GT_BP=coxph(formula = Surv(time, status) ~ strata, data=subset(df_combined_BP, df_combined_BP$subgroup==0)) %>% cox.zph
          GT_Maha=coxph(formula = Surv(time, status) ~ strata, data=subset(df_combined_Maha, df_combined_Maha$subgroup==0)) %>% cox.zph
          GT_Logit=coxph(formula = Surv(time, status) ~ strata, data=subset(df_combined_Logit, df_combined_Logit$subgroup==0)) %>% cox.zph
  
          # compare within rmst model
          df_BP_subgroup0=subset(df_combined_BP, df_combined_BP$subgroup==0)
          df_Maha_subgroup0=subset(df_combined_Maha, df_combined_Maha$subgroup==0)
          df_Logit_subgroup0=subset(df_combined_Logit, df_combined_Logit$subgroup==0)
  
          rmst_BP=rmst2(df_BP_subgroup0$time, df_BP_subgroup0$status, ifelse(df_BP_subgroup0$strata=="matched", 0,1), tau=NULL)
          rmst_Maha=rmst2(df_Maha_subgroup0$time, df_Maha_subgroup0$status, ifelse(df_Maha_subgroup0$strata=="matched", 0,1), tau=NULL)
          rmst_Logit=rmst2(df_Logit_subgroup0$time, df_Logit_subgroup0$status, ifelse(df_Logit_subgroup0$strata=="matched", 0,1), tau=NULL)
  
          # rmst per arm
          df_original_subgroup0_1=subset(df_original, df_original$subgroup==0)
          df_original_subgroup0_2=subset(df_original, df_original$subgroup==0)
          df_original_subgroup0_1$strata=0
          df_original_subgroup0_2$strata=1
          df_original_subgroup0_12=bind_rows(df_original_subgroup0_1, df_original_subgroup0_2)
  
          df_BP_subgroup0_1=subset(df_BP_subgroup0, df_BP_subgroup0$strata=="matched")
          df_BP_subgroup0_2=subset(df_BP_subgroup0, df_BP_subgroup0$strata=="matched")
          df_BP_subgroup0_1$strata=0
          df_BP_subgroup0_2$strata=1
          df_BP_subgroup0_12=bind_rows(df_BP_subgroup0_1, df_BP_subgroup0_2)
  
          df_Maha_subgroup0_1=subset(df_Maha_subgroup0, df_Maha_subgroup0$strata=="matched")
          df_Maha_subgroup0_2=subset(df_Maha_subgroup0, df_Maha_subgroup0$strata=="matched")
          df_Maha_subgroup0_1$strata=0
          df_Maha_subgroup0_2$strata=1
          df_Maha_subgroup0_12=bind_rows(df_Maha_subgroup0_1, df_Maha_subgroup0_2)
  
          df_Logit_subgroup0_1=subset(df_Logit_subgroup0, df_Logit_subgroup0$strata=="matched")
          df_Logit_subgroup0_2=subset(df_Logit_subgroup0, df_Logit_subgroup0$strata=="matched")
          df_Logit_subgroup0_1$strata=0
          df_Logit_subgroup0_2$strata=1
          df_Logit_subgroup0_12=bind_rows(df_Logit_subgroup0_1, df_Logit_subgroup0_2)
  
          rmst1_original=rmst2(df_original_subgroup0_12$time, df_original_subgroup0_12$status, df_original_subgroup0_12$strata, tau=NULL)
          rmst1_BP=rmst2(df_BP_subgroup0_12$time, df_BP_subgroup0_12$status, df_BP_subgroup0_12$strata, tau=NULL)
          rmst1_Maha=rmst2(df_Maha_subgroup0_12$time, df_Maha_subgroup0_12$status, df_Maha_subgroup0_12$strata, tau=NULL)
          rmst1_Logit=rmst2(df_Logit_subgroup0_12$time, df_Logit_subgroup0_12$status, df_Logit_subgroup0_12$strata, tau=NULL)
  
  
          #### 5_Export summary statistics ####
          rbind(
            c(
              n=n[i.n],
              subgroup.p=subgroup.p[i.subgroup],
              censorship_overall.p=censor_overall.p[i.censor_overall],
              censorship_subgroup.p=censor_subgroup.p[i.censor_subgroup],
              missing.p=missing.p[i.missing],
              interval=interval[i.interval],
              mc=i.mc,
              logrank=sum.cox_BP$sctest[3],
              HR=sum.cox_BP$coefficients[2],
              TE=sum.cox_BP$coefficients[1],
              se=sum.cox_BP$coefficients[3],
              GT.p=GT_BP$table[1,3],
              RMSTD=rmst_BP$unadjusted.result[1,1],
              matching="Bipartite"),
  
            c(
              n=n[i.n],
              subgroup.p=subgroup.p[i.subgroup],
              censorship_overall.p=censor_overall.p[i.censor_overall],
              censorship_subgroup.p=censor_subgroup.p[i.censor_subgroup],
              missing.p=missing.p[i.missing],
              interval=interval[i.interval],
              mc=i.mc,
              logrank=sum.cox_Maha$sctest[3],
              HR=sum.cox_Maha$coefficients[2],
              TE=sum.cox_Maha$coefficients[1],
              se=sum.cox_Maha$coefficients[3],
              GT.p=GT_Maha$table[1,3],
              RMSTD=rmst_Maha$unadjusted.result[1,1],
              matching="Mahalanobis"),
  
            c(
              n=n[i.n],
              subgroup.p=subgroup.p[i.subgroup],
              censorship_overall.p=censor_overall.p[i.censor_overall],
              censorship_subgroup.p=censor_subgroup.p[i.censor_subgroup],
              missing.p=missing.p[i.missing],
              interval=interval[i.interval],
              mc=i.mc,
              logrank=sum.cox_Logit$sctest[3],
              HR=sum.cox_Logit$coefficients[2],
              TE=sum.cox_Logit$coefficients[1],
              se=sum.cox_Logit$coefficients[3],
              GT.p=GT_Logit$table[1,3],
              RMSTD=rmst_Logit$unadjusted.result[1,1],
              matching="Logistic")
          )
  
  
  
        },
        error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
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
    summary=unlist(summary(subset(df_analysis, df_analysis$matching==unique(df_analysis$matching)[i])$TE))
    
    abs.mean=subset(df_analysis, df_analysis$matching==unique(df_analysis$matching)[i])$TE %>% abs %>% mean
    abs.sd=subset(df_analysis, df_analysis$matching==unique(df_analysis$matching)[i])$TE %>% abs %>% sd
    abs.summary=unlist(summary(abs(subset(df_analysis, df_analysis$matching==unique(df_analysis$matching)[i])$TE)))
  
    table_cox=rbind(table_cox,c(
      Matching=unique(df_analysis$matching)[i],
      summary,
      sd=sd,
      absTE=abs.summary,
      abs.sd=abs.sd
      ))
  
  }
  
  table_cox=data.frame(table_cox)
  table_cox[,-grep("Matching", colnames(table_cox))]=apply(table_cox[,-grep("Matching", colnames(table_cox))],2, as.numeric)
  
  table_rmst=NULL
  for (i in 1:length(unique(df_analysis$matching))){
    mean=subset(df_analysis, df_analysis$matching==unique(df_analysis$matching)[i])$RMSTD %>% mean
    sd=subset(df_analysis, df_analysis$matching==unique(df_analysis$matching)[i])$RMSTD %>% sd
    summary=unlist(summary(subset(df_analysis, df_analysis$matching==unique(df_analysis$matching)[i])$RMSTD))
    
    abs.mean=subset(df_analysis, df_analysis$matching==unique(df_analysis$matching)[i])$RMSTD %>% abs %>% mean
    abs.sd=subset(df_analysis, df_analysis$matching==unique(df_analysis$matching)[i])$RMSTD %>% abs %>% sd
    abs.summary=unlist(summary(abs(subset(df_analysis, df_analysis$matching==unique(df_analysis$matching)[i])$RMSTD)))
  
    table_rmst=rbind(table_rmst,c(
      Matching=unique(df_analysis$matching)[i],
      summary,
      sd=sd,
      absRMSTD=abs.summary,
      abs.sd=abs.sd
      ))
  
  }
  
  table_rmst=data.frame(table_rmst)
  table_rmst[,-grep("Matching", colnames(table_rmst))]=apply(table_rmst[,-grep("Matching", colnames(table_rmst))],2, as.numeric)
  
  
  #### Histograms ####
  
  plot_cox=ggplot(df_analysis, aes(x=TE, fill=matching))+
    geom_density(alpha=0.7)+
    # geom_histogram(aes(y=..density..), position="identity", alpha=0.3)+
    scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
    geom_vline(xintercept = 0, linetype="dashed", col="red")+
    labs(title=paste("Imputed vs original (marginal cox model)",fig.label),x="ln(HR)", fill="Matching algorithm")+
    geom_vline(data=table_cox, aes(xintercept=Mean), linetype="dashed")
  
  
  plot_rmst=ggplot(df_analysis, aes(x=RMSTD, fill=matching))+
    geom_density(alpha=0.7)+
    # geom_histogram(aes(y=..density..), position="identity", alpha=0.3)+
    scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
    geom_vline(xintercept = 0, linetype="dashed", col="red")+
    labs(title=paste("Imputed vs original (RMST)",fig.label),x="RMST-Differences", fill="Matching algorithm")+
    geom_vline(data=table_rmst, aes(xintercept=Mean), linetype="dashed")
  
  
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
  
  out$density.cox=plot_cox
  out$density.rmst=plot_rmst
  out$convergence=plot_convergence
  
  out$table.cox=table_cox
  out$table.rmst=table_rmst
  
  out$data=df_analysis
  
  class(out)="KMSubtractionError"
  
  out
  
}

