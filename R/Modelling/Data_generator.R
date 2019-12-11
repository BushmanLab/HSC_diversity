## Function that create the simulated dataset

Dataset_creation=function(prop_type,cst_var,mu=3,cst_sampling=1e3,crossOverReports=crossOverReports,indic=indic,conta_prop=conta_prop,loop=FALSE){
  
  # * prop_type: proportion for each type of IS (GMBKT,GM,T) to define an homogeneous or heterogeneous HSPC population.
  # * cst_var, mu are simulation constants : see Step 2 in paragraph 9.2 of Supplementary methods.
  # * cst_sampling is the % of total cells to select to simulate the random sparse sampling.
  # * crossOverReports,indic and conta_prop are real values giving the cells proportion, the contamination matrix etc.
  # * loop is a boolean value. By default the value is FALSE to use a seed for reproductibility.
  #   In the case of the MONTE-CARLO loop the value will be TRUE.
  
  ####################################
  ####  Whole blood IS population #### 
  ####################################
  
  ########  STEP 1 to STEP3 clonal tracking modelling - Supplemental methods 9.2 #########
  
  #Fix the result for reproductibility 
  if(loop == FALSE) {
    set.seed(123)
    print(paste(c("cst_variance :",cst_var,'||',
                  "proportion of GMBKT IS:",prop_type[1],'||',
                  "proportion of GM IS:",prop_type[2],'||',
                  "proportion of T IS:",prop_type[3],'||',
                  "cst for sampling :",cst_sampling),collapse = " "))}
  
  #Nb of IS to generate (see STEP 1 clonal tracking modelling - Supplemental methods 9.2)
  pop_size= indic$Chao1
  
  #Abundance parameters (see STEP 2 clonal tracking modelling - Supplemental methods 9.2)
  mu_data= indic$mean_cells
  v_data = indic$var 
  cst_mean= mu*10^3
  mu_sim = mu_data * cst_mean 
  v_sim = v_data * cst_mean^2 *cst_var
  
  ## Calcul of a low_threshold to avoid missing cell type
  low_threshold = 1 / min(crossOverReports$cellCount)  
  
  ## Calcul of the second parameter of the rnegbin (theta)
  
  #1) overdispersion
  sur_dispersion = function(mu, v ){
    sur_dis = v/mu -1
    return(sur_dis)
  }
  surdis=sur_dispersion(mu_sim + low_threshold,v_sim)
  
  #2) theta
  theta_from_mu_surdis = function(mu,surdis){
    return(mu/surdis)
  }
  theta=theta_from_mu_surdis(mu_sim,surdis)
  
  #Get each cell proportion according the type (see STEP 3 clonal tracking modelling - Supplemental methods 9.2)
  prob_cells=crossOverReports$cellCount
  type_cells=list()
  type_cells[['GMBKT']]=1:5
  type_cells[['GM']]=2:3
  type_cells[['T']]=5

  # ~ STEP 1-2-3 Clonal tracking modelling : whole blood simulation ~
  
  #Create a dataframe with a posid and a column 'gene' to match the real values format for the pipeline.
  #IS_abund correspond to the simulated abundance obtained with the rnegbin function.
  pop_df = tibble(posid = c(1:pop_size) , gene = rep(NA,pop_size) , IS_abund = ceiling(low_threshold) + rnegbin(n=pop_size,mu=mu_sim,theta=theta))
  
  #Determine the type of IS (GMBKT-GM-T) randomly determined with the prop_type frequencies 
  types = sample(c('GMBKT','GM','T'),size= pop_size , prob = prop_type, replace = T) 
  
  #Regroup the data and create empty column to receive the cell type repartition
  pop_df = mutate(pop_df , type = types)
  pop_df = mutate(pop_df , Bcells = 'NA', Monocytes = 'NA', Granulocytes = 'NA', NKcells = 'NA', Tcells= 'NA')
  
  #Cell type repartition for each IS
  
  ##loop & filter on the IS type
  for (typ in unique(types)){
    pop_df_filtered = filter(pop_df, type == typ )
    all_cells =do.call(rbind,lapply(1:nrow(pop_df_filtered),function(i){
    
    ##The total abundance is distributed with a multinomial function, according the prob_cells  
    val_bycell = t(rmultinom(n=1,size=pop_df_filtered$IS_abund[i],prob_cells))
      
    #erase the cell not present in the type
    num=type_cells[[typ]]
    val_bycell[-num]=0
    val_bycell
    }))
    pop_df[pop_df$type == typ,cells] = all_cells
  }
  
  rm(pop_df_filtered)
  rm(all_cells)
  
  ## The pop_df dataset represents the whole blood data.
  ## It is composed of a posid, a gene column, a total abundance (IS_abund), an IS type (GMBKT,GM or T) and an abundance for each cell type.
  pop_df = mutate(pop_df , Bcells = as.numeric(Bcells),
                  Monocytes =as.numeric(Monocytes),
                  Granulocytes =as.numeric(Granulocytes),
                  NKcells = as.numeric(NKcells),
                  Tcells = as.numeric(Tcells))

  
  #Calcul of differents indicators on this dataset.
  indic_values=data.frame(NBofIS=nrow(pop_df),Mean=mean(rowSums(pop_df[cells])),Var=var(rowSums(pop_df[cells])))
  prop=colSums(pop_df[,cells])/sum(pop_df[,cells])
  mean=c(indic$mean_cells,indic_values$Mean)
  var=c(indic$var,indic_values$Var)
  
  ####################################
  ########  Errors modelling ######### 
  ####################################
  
  ########  STEP 4 clonal tracking modelling - Supplemental methods 9.2 #########
  
  #We simulate the different errors introduced during the procedure#.
  
  ## **Error 1) Random sparse sampling **
  
  ## Change the format of the data to work on each cell (long format) and
  #  select only a fraction of the total cells to imitate the blood test effect (n_cells/cst_sampling).
  data_B=pop_df  %>% dplyr::select(posid,Bcells) %>% uncount(Bcells)  %>% mutate(celltype='Bcells') 
  n_cells = nrow(data_B)
  data_B = data_B %>% sample_n(.,size=n_cells/cst_sampling)
  data_M=pop_df  %>% dplyr::select(posid,Monocytes) %>% uncount(Monocytes)  %>% mutate(celltype='Monocytes') 
  n_cells = nrow(data_M)
  data_M = data_M %>% sample_n(size=n_cells/cst_sampling)
  data_G=pop_df  %>% dplyr::select(posid,Granulocytes) %>% uncount(Granulocytes)  %>% mutate(celltype='Granulocytes') 
  n_cells = nrow(data_G)
  data_G = data_G %>% sample_n(size=n_cells/cst_sampling)
  data_K=pop_df  %>% dplyr::select(posid,NKcells) %>% uncount(NKcells)  %>% mutate(celltype='NKcells') 
  n_cells = nrow(data_K)
  data_K = data_K %>% sample_n(size=n_cells/cst_sampling)
  data_T=pop_df  %>% dplyr::select(posid,Tcells) %>% uncount(Tcells)  %>% mutate(celltype='Tcells') 
  n_cells = nrow(data_T)
  data_T = data_T %>% sample_n(size=n_cells/cst_sampling)
  
  bias1 = bind_rows(data_B,data_M,data_G,data_K,data_T)
  
  rm(data_B,data_M,data_G,data_K,data_T,n_cells)
  
  #Calcul of differents indicators on this dataset.
  is_val=bias1 %>% dplyr::select(celltype,posid)  %>% group_by(posid) %>% dplyr::summarise(cellcount=n())
  indic_values=data.frame(NBofIS=nrow(is_val),Mean=mean(is_val$cellcount),Var=var(is_val$cellcount))
  mean=c(mean,indic_values$Mean)
  var=c(var,indic_values$Var)
  
  ## **Error 2) Contamination during cell sorting**  
  
  # Change into a count table format (wide format) with a row per IS.
  data_wide = bias1 %>% group_by(posid) %>% count(celltype) %>% spread(celltype,n) %>% replace_na(list(Bcells  = 0, Granulocytes = 0, Monocytes = 0, NKcells =0, Tcells = 0))
  
  #Get the contamination proportion
  obs_conta = matrix(0,nrow=nrow(data_wide),ncol=5)
  names(obs_conta) = cells
  for (x in 1:5){
    #select the cell type and put a filter on it
    cell=cells[x]
    #select the contamination % of this cell type
    conta_prob=conta_prop[,x]
    data_cell = data_wide %>% dplyr::select(posid,!!cell)
    #apply the contamination
    obs_conta = obs_conta + t(sapply(pull(data_cell,cell),rmultinom,n=1,prob=conta_prob))
  }
  
  #Add the contaminated values ("the cells in the tube") next to the original values
  data_wide$Bcells_cont = obs_conta[,1] 
  data_wide$Monocytes_cont = obs_conta[,2] 
  data_wide$Granulocytes_cont = obs_conta[,3]
  data_wide$NKcells_cont = obs_conta[,4] 
  data_wide$Tcells_cont = obs_conta[,5] 
  data_wide$cellcount_cont = rowSums(obs_conta)
  
  #Change into a "long-format" : 1 row per cell.
  #We now consider the "tube" cells instead of the original cells.
  data_B=data_wide  %>% dplyr::select(posid,Bcells_cont) %>% uncount(Bcells_cont)  %>% mutate(tube='Bcells') 
  data_M=data_wide  %>% dplyr::select(posid,Monocytes_cont) %>% uncount(Monocytes_cont)  %>% mutate(tube='Monocytes') 
  data_G=data_wide  %>% dplyr::select(posid,Granulocytes_cont) %>% uncount(Granulocytes_cont)  %>% mutate(tube='Granulocytes') 
  data_K=data_wide  %>% dplyr::select(posid,NKcells_cont) %>% uncount(NKcells_cont)  %>% mutate(tube='NKcells') 
  data_T=data_wide  %>% dplyr::select(posid,Tcells_cont) %>% uncount(Tcells_cont)  %>% mutate(tube='Tcells') 
  
  bias2 = bind_rows(data_B,data_M,data_G,data_K,data_T)
  rm(data_B,data_M,data_G,data_K,data_T,obs_conta)
  
  #Calcul of differents indicators on this dataset.
  is_val=bias2 %>% dplyr::select(tube,posid)  %>% group_by(posid) %>% dplyr::summarise(cellcount=n())
  indic_values=data.frame(NBofIS=nrow(is_val),Mean=mean(is_val$cellcount),Var=var(is_val$cellcount))
  mean=c(mean,indic_values$Mean)
  var=c(var,indic_values$Var)
 
  ## **Error 3) Unbalanced sampling between cell types (VCN and DNA input)**
  DNA=crossOverReports$DNA
  VCN=matrix(nrow=5,ncol=5,rep(0,5*5))
  diag(VCN)=crossOverReports$VCN
  DNA_freq=(DNA/0.0066)*VCN
  DNA_freq=DNA_freq/sum(DNA_freq)
  
  #Get the frequency of each cell type.
  prop2=bias2 %>% dplyr::select(tube,posid) %>% group_by(tube) %>%
    dplyr::summarise(n=n()) %>% mutate(freq=n/sum(n))  %>% slice(match(tube, cells))
  colnames(prop2)=c("tube",'n','freq')
  prop2= prop2 %>% 
    slice(match(tube, cells)) %>% mutate(expected=diag(DNA_freq))
  
  #Get the smaller cell type.
  chosen_one=prop2[which.min(prop2$freq),]
  #Estimate the nb of cells to take.
  nbcells=round(chosen_one$n/chosen_one$expected)
  prop2$nb_to_select=round(prop2$expected*nbcells)
  
  #Adjust the quantity of cells.
  bias3=do.call(rbind,lapply(1:5,function(x){
    cell=cells[x]
    data_filter=filter(bias2,tube==cell)
    #select a certain frequency
    freq=filter(prop2,tube==cell)
    freq=freq$nb_to_select
    data_filtered= data_filter %>% ungroup() %>% sample_n(round(freq)) 
  }))
  
  
  #Calcul of differents indicators on this dataset.
  is_val= bias3 %>% dplyr::select(tube,posid)  %>% group_by(posid) %>% dplyr::summarise(cellcount=n())
  indic_values=data.frame(NBofIS=nrow(is_val),Mean=mean(is_val$cellcount),Var=var(is_val$cellcount))
  mean=c(mean,indic_values$Mean)
  var=c(var,indic_values$Var)
  
  pdataset=colSums(intSites_new_concat_sup0[,cells])/sum(intSites_new_concat_sup0[,'cellCount'])

  ## Go back to a wide format (a row per IS).
  liste_sampled = mutate(bias3, posid = as.numeric(posid))
  sampled = liste_sampled %>% group_by(posid,tube) %>% dplyr::summarise(nombre = n()) %>% as.data.frame()
  data_sim = sampled %>% spread(tube,nombre) %>% arrange(posid)
  data_sim[is.na(data_sim)] <- 0
  data_sim$cellCount = rowSums(data_sim[,2:6])
  data_sim=left_join(data_sim,pop_df[,c('posid','type','gene')])
  data_sim= dplyr::select(data_sim,posid,gene,Bcells,Monocytes,Granulocytes,NKcells,Tcells,cellCount,type)
  ## data_sim is the dataset representing the real Initial data with the different errors introduced during the procedure.
  
  return(list(pop_df=pop_df,data_sim=data_sim,repart=cbind(mean,var)))
}


