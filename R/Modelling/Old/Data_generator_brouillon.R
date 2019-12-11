### Table printer
Table_printer=function(data,step,var){
  
if (var=='tube'){data=data[,c('tube','posid')]}
colnames(data)=c("celltype","posid")
  
is_val=data %>% group_by(posid) %>% dplyr::summarise(cellcount=n())
indic_values=data.frame(NBofIS=nrow(is_val),Mean=mean(is_val$cellcount),Var=var(is_val$cellcount))
#print(kable(indic_values,caption=paste(step,'- indicators'))%>%
       # kable_styling())

## Print the cells proportions
prop=data %>% group_by(celltype) %>% dplyr::summarise(n=n()) %>% mutate(freq=n/sum(n))  %>% slice(match(celltype, cells))
#print(kable(t(prop),caption=paste(step,'- Cell types proportions:'))%>%
       # kable_styling())

return(list(prop=prop,indic_values=indic_values))
}


## Data creation
Dataset_creation=function(prob_type,cst_var,cst_ab=1e3,crossOverReports=crossOverReports,mu=3,indic=indic,loop=FALSE){
  # prob_type: definir .... TODO
  # cst_var, mu are simulation constants : see Step 2 in paragraph 8.2 of Supplementary methods
  #********* Modelisation *****************
  
  #Fix the result for reproductibility 
  if(loop == FALSE) {set.seed(123)}
  
  #Type of IS parameters
  prob_cells=crossOverReports$cellCount
  type_cells=list()
  type_cells[['GMBKT']]=1:5
  type_cells[['GM']]=2:3
  type_cells[['T']]=5
  
  # Parameters to create the total pop. dataframe
  pop_size= indic$Chao1 #Nb of IS to generate
  mu_data= indic$mean_cells #Mean of the abundance
  v_data = indic$var #Variance of the abundance
  cst_mean= mu*10^3
  
  ## calcul seuil bas pour que (en espérance il y ait au moins une cellule de chaque type)
  seuil_bas = 1 / min(crossOverReports$cellCount)  
  
  
  mu_sim = mu_data * cst_mean 
  v_sim = v_data * cst_mean^2 *cst_var
  

  
  ## calcul of the second parameter of the rnegbin (theta)
  
  #1) overdispersion
  sur_dispersion = function(mu, v ){
    sur_dis = v/mu -1
    return(sur_dis)
  }
  surdis=sur_dispersion(mu_sim + seuil_bas,v_sim)
  
  #2) theta
  theta_from_mu_surdis = function(mu,surdis){
    return(mu/surdis)
  }
  theta=theta_from_mu_surdis(mu_sim,surdis)
  
  
  #***IS creation***
  pop_df = tibble(posid = c(1:pop_size) , gene = rep(NA,pop_size) , IS_abund = ceiling(seuil_bas) + rnegbin(n=pop_size,mu=mu_sim  ,theta=theta))
  types = sample(c('GMBKT','GM','T'),size= pop_size , prob = prob_type, replace = T) 
  pop_df = mutate(pop_df , type = types)
  pop_df = mutate(pop_df , Bcells = 'NA', Monocytes = 'NA', Granulocytes = 'NA', NKcells = 'NA', Tcells= 'NA')
  
  for (typ in unique(types)){
    pop_df_filtered = filter(pop_df, type == typ )
    
    all_cells =do.call(rbind,lapply(1:nrow(pop_df_filtered),function(i){
      val_bycell = t(rmultinom(n=1,size=pop_df_filtered$IS_abund[i],prob_cells))
      
      #erase the cell not present in the type
      num=type_cells[[typ]]
      val_bycell[-num]=0
      val_bycell
      
    }))
    
    pop_df[pop_df$type == typ,cells] = all_cells
  }
  
  rm(pop_df_filtered)
  pop_df = mutate(pop_df , Bcells = as.numeric(Bcells),
                  Monocytes =as.numeric(Monocytes),
                  Granulocytes =as.numeric(Granulocytes),
                  NKcells = as.numeric(NKcells),
                  Tcells = as.numeric(Tcells))
  rm(all_cells)
  #Show the 3 first lines
  
  #print(kable(head(pop_df,3),caption='Total population : all blood')%>%
  # kable_styling())
  
  # Show differents indicators
  indic_values=data.frame(NBofIS=nrow(pop_df),Mean=mean(rowSums(pop_df[cells])),Var=var(rowSums(pop_df[cells])))
  # print(kable(indic_values,caption='Total population indicators')%>%
  #   kable_styling())
  
  prop=colSums(pop_df[,cells])/sum(pop_df[,cells])
  #print(kable(t(prop),caption='Cell types proportions:')%>%
  #  kable_styling())
  
  #keep some values
  mean=c(indic$mean_cells,indic_values$Mean)
  var=c(indic$var,indic_values$Var)
  
  ## ***(B) Add the bias*** ##
  
  #cd ~
  #touch .Renviron
  #open .Renviron
  #R_MAX_VSIZE=100Gb 
  
  
  ## **1) Sparse sampling**
  
  ## Change the format of the data to work on each cell (long format) and
  ##### Sparse sampling on the abundance --> blood test effect#####
  data_B=pop_df  %>% dplyr::select(posid,Bcells) %>% uncount(Bcells)  %>% mutate(celltype='Bcells') 
  n_cells = nrow(data_B)
  data_B = data_B %>% sample_n(.,size=n_cells/cst_ab)
  data_M=pop_df  %>% dplyr::select(posid,Monocytes) %>% uncount(Monocytes)  %>% mutate(celltype='Monocytes') 
  n_cells = nrow(data_M)
  data_M = data_M %>% sample_n(size=n_cells/cst_ab)
  data_G=pop_df  %>% dplyr::select(posid,Granulocytes) %>% uncount(Granulocytes)  %>% mutate(celltype='Granulocytes') 
  n_cells = nrow(data_G)
  data_G = data_G %>% sample_n(size=n_cells/cst_ab)
  data_K=pop_df  %>% dplyr::select(posid,NKcells) %>% uncount(NKcells)  %>% mutate(celltype='NKcells') 
  n_cells = nrow(data_K)
  data_K = data_K %>% sample_n(size=n_cells/cst_ab)
  data_T=pop_df  %>% dplyr::select(posid,Tcells) %>% uncount(Tcells)  %>% mutate(celltype='Tcells') 
  n_cells = nrow(data_T)
  data_T = data_T %>% sample_n(size=n_cells/cst_ab)
  
  bias1 = bind_rows(data_B,data_M,data_G,data_K,data_T)
  
  rm(data_B,data_M,data_G,data_K,data_T,n_cells)
  
  
  
  #head(bias1 %>% group_by(posid) %>% count(celltype))
  
  
  ## Print the indicators
  tab=Table_printer(data=bias1,step='Bias 1',var="celltype")
  
  mean=c(mean,tab$indic_values$Mean)
  var=c(var,tab$indic_values$Var)
  ## **2) Contamination bias**  
  data_wide = bias1 %>% group_by(posid) %>% count(celltype) %>% spread(celltype,n) %>% replace_na(list(Bcells  = 0, Granulocytes = 0, Monocytes = 0, NKcells =0, Tcells = 0))
  obs_conta = matrix(0,nrow=nrow(data_wide),ncol=5)
  names(obs_conta) = cells
  for (x in 1:5){
    #select the cell type and put a filter on it
    cell=cells[x]
    conta_prob=conta_prop[,x]
    data_cell = data_wide %>% dplyr::select(posid,!!cell)
    obs_conta = obs_conta + t(sapply(pull(data_cell,cell),rmultinom,n=1,prob=conta_prob))
  }
  
  data_wide$Bcells_conta = obs_conta[,1] 
  data_wide$Monocytes_cont = obs_conta[,2] 
  data_wide$Granulocytes_cont = obs_conta[,3]
  data_wide$NKcells_cont = obs_conta[,4] 
  data_wide$Tcells_cont = obs_conta[,5] 
  data_wide$cellcount_cont = rowSums(obs_conta )
  
  
  ### HERE#####
  data_B=data_wide  %>% dplyr::select(posid,Bcells_conta) %>% uncount(Bcells_conta)  %>% mutate(tube='Bcells') 
  data_M=data_wide  %>% dplyr::select(posid,Monocytes_cont) %>% uncount(Monocytes_cont)  %>% mutate(tube='Monocytes') 
  data_G=data_wide  %>% dplyr::select(posid,Granulocytes_cont) %>% uncount(Granulocytes_cont)  %>% mutate(tube='Granulocytes') 
  data_K=data_wide  %>% dplyr::select(posid,NKcells_cont) %>% uncount(NKcells_cont)  %>% mutate(tube='NKcells') 
  data_T=data_wide  %>% dplyr::select(posid,Tcells_cont) %>% uncount(Tcells_cont)  %>% mutate(tube='Tcells') 
  
  bias2 = bind_rows(data_B,data_M,data_G,data_K,data_T)
  
  rm(data_B,data_M,data_G,data_K,data_T,obs_conta)
  
  
  ## Print the indicators
  tab=Table_printer(data=bias2,step='Bias 2',var="tube")
  
  #keep some values
  mean=c(mean,tab$indic_values$Mean)
  var=c(var,tab$indic_values$Var)
  prop2=tab$prop
  
  ## **3) DNA quantity**
  
  DNA_freq=(DNA/0.0066)*VCN
  DNA_freq=DNA_freq/sum(DNA_freq)
  
  #reorder prop2
  
  colnames(prop2)=c("tube",'n','freq')
  prop2= prop2 %>% 
    slice(match(tube, cells)) %>% mutate(expected=diag(DNA_freq))
  
  chosen_one=prop2[which.min(prop2$freq),]
  nbcells=round(chosen_one$n/chosen_one$expected)
  prop2$nb_to_select=round(prop2$expected*nbcells)
  
  
  bias3=do.call(rbind,lapply(1:5,function(x){
    
    cell=cells[x]
    data_filter=filter(bias2,tube==cell)
    
    #select a certain frequence
    freq=filter(prop2,tube==cell)
    freq=freq$nb_to_select
    data_filtered= data_filter %>% ungroup() %>% sample_n(round(freq)) 
  }))
  
  
  ## Print the indicators
  tab=Table_printer(data=bias3,step='Bias 3',var="tube")
  
  #keep some values
  mean=c(mean,tab$indic_values$Mean)
  var=c(var,tab$indic_values$Var)
  
  pdataset=colSums(intSites_new_concat_sup0[,cells])/sum(intSites_new_concat_sup0[,'cellCount'])
  #print(kable(pdataset,caption='Real Dataset - Cell types proportions (pdataset):')%>%
  #      kable_styling())
  
  
  ## Go back to a wide format
  liste_sampled = mutate(bias3, posid = as.numeric(posid))
  sampled = liste_sampled %>% group_by(posid,tube) %>% dplyr::summarise(nombre = n()) %>% as.data.frame()
  data_sim = sampled %>% spread(tube,nombre) %>% arrange(posid)
  data_sim[is.na(data_sim)] <- 0
  data_sim$cellCount = rowSums(data_sim[,2:6])
  data_sim=left_join(data_sim,pop_df[,c('posid','type','gene')])
  data_sim= dplyr::select(data_sim,posid,gene,Bcells,Monocytes,Granulocytes,NKcells,Tcells,cellCount,type)
  
  return(list(pop_df=pop_df,data_sim=data_sim,repart=cbind(mean,var)))
  

}


