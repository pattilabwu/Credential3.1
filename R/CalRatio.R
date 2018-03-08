
#Calratio

calrat = function(dt){

  li = which.min(dt$mz)
  hi = which.max(dt$mz)
  ro = dt$i[li]/dt$i[hi]  #ratio between smallest and largest mz (12C/13C)
  return(ro)
}

Calratio = function(features_quipu,Credentials_quipu){

  for( i in seq_len(nrow(Credentials_quipu[!is.na(quipu)]))){

    idx = which(features_quipu[,"quipu"] ==i)
    Credentials_quipu[quipu == i,"ratio"] = calrat(features_quipu[idx])
  }
  return(Credentials_quipu)

}

