
### ****** Clustering of HSPC lineage output ******###
### ****** Supplemental Methods 7 ******###

########################################################
### Functions used to apply the k-means clustering  ### 
#######################################################

Transformation <- function(data_prop_init,correction){
  
  data = data_prop_init[,grepl("_corrected_prop",colnames(data_prop_init))]
  
  Moyenne_Geo= geometric.mean(t(data)) 
  dataprop_mg = round(data/Moyenne_Geo , 4)
  
  if (correction == "logCLR")  {
    dataprop_logCLR = matrix(NA, nrow(data), ncol(data))
    for(i in 1:nrow(data)) {
      for(j in 1:ncol(data)) {
        if (dataprop_mg[i,j]<=1)
          dataprop_logCLR[i,j] = round(-log(1 - log(dataprop_mg[i,j]))^2 , 4)
        else dataprop_logCLR[i,j] = round((log(dataprop_mg[i,j]))^2 ,4 )
      }
    }
    VFdataprop_transf  = cbind(data_prop_init , dataprop_logCLR)
  }
  else{
    VFdataprop_transf  = cbind(data_prop_init , log(dataprop_mg))
  }
  nom = names(data)
  names(VFdataprop_transf) = c(names(data_prop_init),paste0(strsplit(nom,"_corrected_prop"),paste0("_transf_",correction)))
  return(VFdataprop_transf)
}

#--------------------------------------------------#

## from CLR or logCLR corrected data to proportions
return_to_prop = function(tableau,transf="logCLR"){
  n_clust = nrow(tableau)
  n_cells = ncol(tableau)
  tableau_prop = matrix(nrow=n_clust,ncol=(n_cells))
  if (transf=="logCLR"){
    for (i in 1:n_clust){
      for (j in 1:n_cells){
        if (tableau[i,j]>0){
          tableau_prop[i,j] = exp( sqrt(tableau[i,j]) )}
        if (tableau[i,j]<=0){
          tableau_prop[i,j] = exp(1 - exp( sqrt(-tableau[i,j])))}
      }
      tableau_prop[i,] = tableau_prop[i,] / sum(tableau_prop[i,])
    }
  }
  if (transf=="CLR"){
    for (i in 1:n_clust){
      for (j in 1:n_cells){
        tableau_prop[i,j] = exp( tableau[i,j]) }
      tableau_prop[i,] = tableau_prop[i,] / sum(tableau_prop[i,])            
    }
  }
  tableau_prop = data.frame(tableau_prop)
  names(tableau_prop) = names(tableau)
  return(tableau_prop)
}

#--------------------------------------------------#


  ############################
  ### clustering by Kmeans ###
  ############################


  KmeansClustering = function(data,transf="logCLR",scaling="FALSE",elbow = "TRUE",HD){
  
  ###Seed for reproducible results
  set.seed(321)
  
  # --- Data transformation --- # 
  data_transformed = Transformation(data,transf)
  
  data_clustering = data_transformed[,grepl("_transf",names(data_transformed))] 
  if (scaling == TRUE){
    data_clustering <- as.matrix(scale(data_clustering)) 
  }
  
  
  # --- Choice of the nb of clusters --- # 
  #Find the optimal number of clusters with 2 methods 
  if (elbow == TRUE){
    # Compute bss (var between the clusters), tss (var total), Rsquare values for k = 1 to k = 15 clusters.
    k.max <- 15
    resk <- sapply(1:k.max, 
                   function(k){kmeans(data_clustering, k, nstart=50,iter.max = 15)})
    
    bss=unlist(resk["betweenss",])
    tss=unlist(resk["totss",])
    Rsquare = bss/tss   #(1) Will be plot later for the elbow criterion
    
    n = nrow(data_clustering)
    d = ncol(data_clustering)
    pen = sqrt(n*d*(1:k.max))
    contrast =(1-bss/tss)
    nbclust_estim=which.min(contrast+(2/n)*pen) #(2) Estimation of the optimal nb of clusters 
    #print(nbclust_estim)
  }
  # --------------------- #
  
  
  # --- Clustering --- # # Clusterize data with the best k 
  nbclust=nbclust_estim
  
  # ************* change the number of clusters of some patients for more accurate results *************
  if (unique(data$key) == 'WAS4|m48|Blood') {nbclust=6}
  if (unique(data$key) == 'WAS4|m12|Blood') {nbclust=4}
  if (unique(data$key) == 'bS/bS|m24|Blood') {nbclust=5}
  if (unique(data$key) == 'pat|mxx|Blood') {nbclust=4}
  
  
  #use the kmeans algorithm
  clust = kmeans(data_clustering, nbclust, nstart=50, iter.max = 15)
  # ------------------ #
  
  ### save data in intSites_withKmeans
  intSites_withKmeans = data_transformed
  intSites_withKmeans$cluster_num = clust$cluster
  
  ### save clusters centers into a file order by abundance
  prop = return_to_prop(clust$centers)
  nom = names(data_clustering)
  names(prop) = paste0(strsplit(nom,paste0("_transf_",transf)),"_prop")
  tableau = data.frame(cbind(clust$centers,prop))
  
  tableau = dplyr::mutate(tableau,clust_size = clust$size)
  tableau = dplyr::mutate(tableau,clust_num = c(1:nbclust))
  tableau$key = unique(data$key)
  tableau = arrange(tableau,-clust_size)
  
  # --- calculate the distances between the centroids and the blood characteristics--- #
  ## Kullback-Leibler Distance
  kl = function(prop1, prop2){
    sum(prop1 * log(prop1/prop2))
  }
  
  ##Data preparation
  celltypes = c('G','M','B','K','T')
  compositions = list()     
  
  key=unique(data$key) 
  tableau = tableau 
  clus = tableau$clust_num
  tableau_dis = tableau[,grepl("_prop",colnames(tableau))]
  tableau_dis=tableau_dis[c(3,2,1,4,5)] #change the order of the cell types (GMBKT)
  
  initial_comp=HD
  nom = paste(celltypes,collapse = "")
  compositions[[nom]] = initial_comp
  
  ##Create the different blood references composition (GMBKT,GMBK,GMBT...,K,T)
  for (m in 1:4){
    aeffacer = combn(1:5,m)
    for (k in 1:ncol(aeffacer)){
      comp = initial_comp
      comp[aeffacer[,k]] = initial_comp[aeffacer[,k]]/10
      comp = comp / sum(comp)
      nom = paste(celltypes[-aeffacer[,k]],collapse = "")
      compositions[[nom]] = comp
    }
  }
  
  ##Calculation of the distances
  distances = matrix(ncol=nbclust,nrow=length(compositions))
  rownames(distances) = names(compositions)
  colnames(distances) = clus
  for (i in 1:nrow(tableau_dis)){
    distances[,i] = unlist(lapply(compositions,kl,tableau_dis[i,]))
  }
  
  distances = as.data.frame(distances)
  distances$key=rep(key,nrow(distances))
  distances$distance =  names(compositions)
  distances = melt(distances,id=c('key','distance'))
  
  ##Find the closest reference composition to the cluster centroid
  dist_min=distances %>% 
    group_by(.dots=c('key','variable')) %>% 
    slice(which.min(value))
  
  ##Define the name of each cluster according to its closest reference
  cluster = rep(NA,nbclust)
  for (i in 1:nrow(tableau)){
    cluster[i] = dist_min$distance[unique(dist_min$key)==unique(tableau$key) & dist_min$variable==tableau$clust_num][i] 
  }
  tableau$cluster = cluster
  # ----------------------------------------------------------------------------------------------------------- #
  
  ##save the clusters names in intSites_withKmeans
  intSites_withKmeans = dplyr::mutate(intSites_withKmeans , 
                                      cluster = plyr::mapvalues(cluster_num, 
                                                                from = tableau$clust_num, 
                                                                to = tableau$cluster))
  
  ##Return the results
  return(list(tableau=tableau ,
              intSites_withKmeans = intSites_withKmeans , 
              R2 = (clust$betweenss/clust$totss)*100,
              elbow = Rsquare,
              nbclust_estim=nbclust_estim
  ))
}
