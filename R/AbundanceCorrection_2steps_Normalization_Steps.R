
  Correction_2steps_Normalization=function(x,d2){

  ### ****** STEP C1 : Corrections for cell sorting error ******
  ### ****** Supplemental Methods 4.2 ******

    #(a) replace all zeros in the crossover table to obtain 0.003% in proportion

    matrix=x$table

    for (r in 1:nrow(matrix)){
      tab=matrix[r,]
      nb_zero=length(tab[tab == 0])
      if (nb_zero >0){
       percent=0.003/100*nb_zero
       p=(percent)/(1-percent)*sum(tab[tab>0])
       remp=p/nb_zero
       tab[tab==0]=remp
       }
    matrix[r,]=tab
    }

    #(b) multiplication with the VCN values

    confusion=as.matrix(matrix)
    VCN=matrix(nrow=5,ncol=5,rep(0,5*5))
    diag(VCN)=x$VCN
    confusion=confusion%*%VCN

    #(c) proportion calul
    Confusion_prof_row = prop.table(confusion,1)

    #(d) apply the correction

    data_to_correct = d2[,1:7]
    data_values <- data_to_correct[,3:7]

    data_correctedC = data.frame()
    data_correctedC = data_to_correct[,1:2]
    data_correctedC[,3:7] = as.matrix(data_values) %*% Confusion_prof_row
    colnames(data_correctedC) = colnames(data_to_correct)


    ### ****** STEP C2 : Normalization for unbalanced sampling (DNA and VCN input) ******
    ### ****** Supplemental Methods 4.3 ******


  DNA=(x$DNA/0.0066*x$VCN)
  DNA_ref= DNA/((1e5))

  weights_Ini_data = DNA_ref
  M=matrix(0,ncol=5,nrow=5)
  diag(M)=weights_Ini_data^(-1) 
  M=as.matrix(M)

  data_correctedC2=  as.matrix(data_correctedC[,3:7])%*% M

  #Put the values in a data frame to return
  data_corrected_final=data_correctedC[,1:2]
  data_corrected_final[,3:7]=data_correctedC2  
  colnames(data_corrected_final) = colnames(data_to_correct)
  
  return(list('data_correctedC'=data_correctedC,'data_corrected_final'=data_corrected_final))

  }
