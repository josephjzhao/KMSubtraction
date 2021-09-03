#' Match patients between subgroup and an overall cohorts.
#'
#' This function conducts matching of patients between subgroup and an overall cohort.
#' @param df_overall Dataframe of the overall cohort containing time-to-event data. Prepared with columns "time" and "status".
#' @param df_subgroup Dataframe of the subgroup cohort containing time-to-event data. Prepared with columns "time" and "status".
#' @param matching Choice of matching algorithm. The default matching algorithm is the minimum cost bipartite matching with the Hungarian algorithm ("bipartite"), users may choose to use Mahalanobis distance matching ("maha") or nearest neighbor matching with distances determined by logistic regression ("logit").
#' @keywords
#' @export
#' @examples
#' data(cancer)
#' df_overall=colon
#' df_subgroup=colon[1:200,]
#' KMSubtractionMatch(df_overall, df_subgroup, matching="bipartite")

KMSubtractionMatch=function(df_overall,
                               df_subgroup,
                               matching="bipartite"
                               ){
  
  if (nrow(df_overall)<nrow(df_subgroup)){
    stop("Number of patients in overall cohort is less than the number of patients in the subgroup cohort.")
  }
  
  # Conduct matching
  if (matching=="bipartite"){
    #events
    overall=subset(df_overall, df_overall$status==1 )$time
    subgroup=subset(df_subgroup, df_subgroup$status==1 )$time
    matrix.events=sapply(subgroup, function(x){x-overall}, simplify=T) %>% abs()
    match.events=HungarianSolver(matrix.events)
    match.events$pairs[,2][match.events$pairs[,2]==0]=NA

    #censorships
    overall=subset(df_overall, df_overall$status==0 )$time
    subgroup=subset(df_subgroup, df_subgroup$status==0 )$time
    matrix.censorships=sapply(subgroup, function(x){x-overall}, simplify=T) %>% abs()
    match.censorships=HungarianSolver(matrix.censorships)
    match.censorships$pairs[,2][match.censorships$pairs[,2]==0]=NA

    df_match=rbind(
      rbind(cbind(subset(df_overall, df_overall$status==1), matchid=match.events$pairs[,2], curve="overall"), cbind(subset(df_subgroup, df_subgroup$status==1), matchid=1, curve="subgroup")),# events
      rbind(cbind(subset(df_overall, df_overall$status==0), matchid=match.censorships$pairs[,2], curve="overall"), cbind(subset(df_subgroup, df_subgroup$status==0), matchid=1, curve="subgroup")) # censors
    )

    df_match$subgroup=ifelse(is.na(df_match$matchid), 0, 1)
    matching="Minimum cost bipartite matching with the Hungarian algorithm"

  }

  if (matching=="logit"){

    df_overall$curve=0
    df_subgroup$curve=1

    df=bind_rows(df_overall, df_subgroup)

    match.events=matchit(curve~time, method="nearest", data=subset(df, df$status==1), ratio=1)
    match.censorships=matchit(curve~time, method="nearest", data=subset(df, df$status==0), ratio=1)

    df_match=rbind(cbind(subset(df, df$status==1), matchid=match.events$subclass),
                   cbind(subset(df, df$status==0), matchid=match.censorships$subclass))

    df_match$subgroup=ifelse(is.na(df_match$matchid), 0, 1)

    df_match$curve=ifelse(df_match$curve==0, "overall", "subgroup")
    matching="Nearest neighbor matching with distances determined by logistic regression"
  }

  if (matching=="maha"){

    df_overall$curve="overall"
    df_subgroup$curve="subgroup"

    df=bind_rows(df_overall, df_subgroup)

    match.events=matchit(curve~time, distance="mahalanobis", data=subset(df, df$status==1), ratio=1)
    match.censorships=matchit(curve~time, distance="mahalanobis", data=subset(df, df$status==0), ratio=1)

    df_match=rbind(cbind(subset(df, df$status==1), matchid=match.events$subclass),
                   cbind(subset(df, df$status==0), matchid=match.censorships$subclass))

    df_match$subgroup=ifelse(is.na(df_match$matchid), 0, 1)
    matching="Mahalanobis distance matching"
  }

  # Parameters of match
  Parameters=data.frame(rbind('Size of overall cohort'=nrow(df_overall),
              'Size of subgroup cohort'=nrow(df_subgroup),
              'Proportion of subgroup'=round(nrow(df_subgroup)/nrow(df_overall), 3),
              'Proportion of censorship, overall cohort'=round(prop.table(table(df_overall$status))[[1]], 3),
              'Proportion of censorship, subgroup cohort'=round(prop.table(table(df_subgroup$status))[[1]], 3)
              ))

  # Cleaning and formatting tables
  Parameters[,1]=as.numeric(Parameters[,1])
  colnames(Parameters)=""

  # Export
  out=list()
  out$'Matching algorithm'=matching
  out$Parameters=Parameters
  class(out)="KMSubtractionMatch"
  out$data=df_match
  out$data_unreportedsubgroup=subset(df_match, df_match$subgroup==0)

  out
}



