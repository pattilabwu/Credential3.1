#' Extract peak-table from an xcmsSet object
#'
#'@param xs a xcmsSet object
#'@param intchoice the intensity options provided by xcms. Options: "into"(default) -raw peak area, "intb" -baseline corrected peak area, "maxo" - peak height
#'@param sampling the way to sampling intensity in each feature group. If 1, sampling mean intensities of all peaks. If 2, sampling the maximum intensity among all peaks.
#'@import data.table
#'@export
#'@return a list includes two data.tables. featuretable is the detailed featuretable. credTable is the simplified format credentialing accepts in the next step.

xsGroupExtract = function(xs,intchoice="into",sampling=1){

  cat("Extracting peaks based on mzmed,rtmed and",intchoice,"...\n" )

  #vetting start
  #xs=xs_1
  #intchoice="into"
  #sampling=1
  #vetting end
  # Extract the intensity based on intchoice and sampling

  xS=data.table(xs@groups) ## vetting

  for(i in seq.int(xS[, .N])){

    idx_tmp = xs@groupidx[[i]]
    peak_tmp = xs@peaks[idx_tmp,] # extract all the intensities and take the average
    if(sampling==1){
      xS[i,"intmean":= mean(peak_tmp[,intchoice])] #intm is the mean of into's of all sample peaks
      intindex = "intmean"
    } else if(sampling==2){ # extract the maximum intensity
      xS[i,"intmax":= max(peak_tmp[,intchoice])]
      intindex = "intmax"
    } else{
      stop("Wrong sampling number. 1--mean intensity, or 2-- maximum intensity", call. = TRUE)
    }
    xS[i,"rtdiff" := max(peak_tmp[,"rt"]) - min(peak_tmp[,"rt"])] # check peak width
    xS[i,"mzppm":= (max(peak_tmp[,"mzmax"])-min(peak_tmp[,"mzmin"]))*10E6/mean(peak_tmp[,"mz"])] # check mz ppm error
  }

  # dataset clean-up
  xS[,"cc":=seq.int(xS[,.N])] #assign a referenece number to each feature
  if(sampling==1){
  dt <- xS[,c("cc","mzmed","rtmed","intmean")]
  } else if(sampling==2)
  dt <- xS[,c("cc","mzmed","rtmed","intmax")]

  colnames(dt) <- c("cc","mz", "rt", "i")

  extractFeatures = list(featuretable=xS,credTable=dt)

  cat("\nFeatures are successfully extracted.")
  cat("\nThe peak width range is between",min(xS[,rtdiff]),"and",max(xS[,rtdiff]),"with an average of",mean(xS[,rtdiff],na.rm = TRUE))
  cat("\nThe mass error in ppm is between",min(xS[,mzppm]),"and",max(xS[,mzppm]),"with an average of",mean(xS[,mzppm]),na.rm =TRUE)
  return(extractFeatures)

}
