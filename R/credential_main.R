#' #' Credentialing main function
#'
#' This is the main function for credentialing.
#'
#' @param peaktable1 the features from the first ratio dataset.
#' @param peaktable2 the features from the second ratio dataset.
#' @param ppm mass error tolerance for isotopologue searching and credentialed features searching.
#' @param rtwin retention time window for 1st round credentialing
#' @param rtcom retention time window for 2nd round credentialing
#' @param ratio1 ratio of designated 12C/13C ratio in the first peaktable, by default set to 1/1
#' @param ratio2 ratio of designated 12C/13C ratio in the first peaktable, by default set to 1/2
#' @param ratio_tol  a factor that controls the ratio range to pass the 1st round credentialing, The default value is 0.1.
#' @param ratios_tol a factor that controls the ratio range to pass the 2nd round credentialing, The default value is 0.1.
#' @param cd unit mass difference among isotopologues. by defalut is 13C-12C = 13.00335 - 12.
#' @param charges charge states of the ions to be considered when searching isotope pairs. The default value is 1:2.
#' @param mpc A range of mass per carbon to be considered when searching possible isotope pairs and identifying U12C and U13C isotope peaks. The default value is c(12,120)
#' @param maxnmer maximum number of ion aggregation to be considered, The default value is 4.
#' @param export Boolean value. Whether to export excel files for the credentialed features. The default is FALSE.
#' @keywords credentialing
#' @import data.table matrixStats
#' @return A list of the credentialed features before and after 1st and 2nd round filtering by the function.
#' @export

credentialing = function(peaktable1, peaktable2, ppm, rtwin, rtcom, ratio1= 1/1, ratio2 = 1/2, ratio_tol=0.1, ratios_tol = 0.1, cd=13.00335-12, charges= 1:2, mpc=c(12,120), maxnmer=4){

  #initiation
  peaktable1 = data.table(peaktable1)
  peaktable2 = data.table(peaktable2)
  # 1st round credentialing

  cat("\n1st round Credentialing on the first peaktable...\n")

  # find possibel isotope head and tails in different charge states
  knots_f1 = findknots(peaktable1, .zs = charges, ppmwid = ppm, rtwid = rtwin, cd = cd)
  # Resolving issues with merged knots
  knots_f1 = fixmergedquipu(knots_f1,peaktable1)
  # heuristic search for knots that satisfying credentialing filters
  credentials_f1 = credentialknots(knots_f1$knot, ppmwid = ppm, rtwid = rtwin, mpc = mpc, ratio = ratio1, ratio.lim = ratio_tol, maxnmer = maxnmer, cd = cd)
  # labelled credentialed features at original peaktable
  ft_1 = peaktable1[knots_f1$cc_knot[credentials_f1$knot_quipu[!is.na(quipu)],on="knot"],,on="cc"]

  # calculate the ratios within each credentialed peak group(quipu)
  credentials_f1$quipu[,"ratio" := 0]

  credentials_f1$quipu = Calratio(ft_1,credentials_f1$quipu)


  cat("\n1st round Credentialing on the second peaktable...\n")

  #credentialing with ratio2 combination
  knots_f2 = findknots(peaktable2, .zs = charges, ppmwid = ppm, rtwid = rtwin, cd = cd)
  knots_f2 = fixmergedquipu(knots_f2,peaktable2)
  credentials_f2 = credentialknots(knots_f2$knot, ppmwid = ppm, rtwid = rtwin, mpc = mpc, ratio = ratio2, ratio.lim = ratio_tol, maxnmer = maxnmer, cd = cd)

  ft_2 = peaktable2[knots_f2$cc_knot[credentials_f2$knot_quipu[!is.na(quipu)],on="knot"],,on="cc"]

  credentials_f2$quipu[,"ratio" := 0]
  credentials_f2$quipu = Calratio(ft_2,credentials_f2$quipu)

  # rearrangement of credentialed features
  dt1 = ft_1[credentials_f1$quipu[,c("quipu","ratio")][!is.na(quipu)],,on="quipu"]
  dt2= ft_2[credentials_f2$quipu[,c("quipu","ratio")][!is.na(quipu)],,on="quipu"]

  dtR1 = dt1[order(dt1[,"quipu"], dt1[,"mz"]),]
  dtR2 = dt2[order(dt2[,"quipu"], dt2[,"mz"]),]

  # 2nd round filtering

  match_cf = matchcredfeature(dtR1,dtR2,ppm=ppm,drt=rtcom,ratio=ratio1/ratio2,ratio_tol = ratios_tol)

  # result showup


  # data output

  credentialing = list(CredentialedFeature1R1=ft_1, CredentialedFeature2R1=ft_2, knots1=knots_f1, knots2=knots_f2, credentialedKnots1=credentials_f1,credentialedKnots2=credentials_f2,
                       CredentialedFeatureGroups = match_cf$Credentialed_FeatureGroups, CredentialedFeatureR2=match_cf$Credentialed_Features, CredentialedFeatureR2F = match_cf$Credentialed_Features_Filtered,
                       CredentialedFeature1N2 = match_cf$NomatchFeatures_Group1, CredentialedFeature2N2 = match_cf$NomatchFeatures_Group2)


  return(credentialing)
}
