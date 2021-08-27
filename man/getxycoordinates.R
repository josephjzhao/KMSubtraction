#' getxycoordinates
#'
#' This function converts a Kaplan Meier curve into xy coordinates. It returns a dataframe containing xy coordinates of the Kaplan Meier curve created with time-to-event data from the dataframe provided.
#' @param df The dataframe containing time and status data required to create the Kaplan Meier curve.
#' @keywords
#' @export
#' @examples

getxycoordinates = function(df){
  scale_points=matrix(c(0,max(df$time),0,max(df$time),0,0,1,1), ncol=2, nrow=4)

  km_overall=survfit(Surv(time, status) ~ 1, data=df)

  # as svg
  svglite(paste0("curve.svg"), width=5, height=5)
  plot(km_overall, col=c("#f042ed"), conf.int = F, lwd=1, xlim=c(0,max(df$time)))
  points(scale_points, col="#f04242", pch=16, cex=0.3)
  dev.off()

  # overall
  im=image_read_svg(paste0("curve.svg"), width=1000)
  unlink(paste0("curve.svg"))
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

  return(df_scaled)

}

