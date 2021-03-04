##
## Funciones auxiliares para trabajar con singleCellExperiments
##
##


## Utility ####

# Funcion para calcular el disparity p-value de un ranking-list
# para reconocer los top-ranked significativoscon alpha < alpha_min
# si alpha_min==0 se devuelve un vector de pvalues
disparity<-function(x,alpha_min=0.05){
  x[is.na(x)]<-0
  xx <- x/sum(x)
  if(alpha_min==0){
    #cat('alpha_min=0, returning disparity-pvalues.\n')
    xx <- (1-xx)^(length(x)-1)
  }else{
    alphaOK <- (1-xx)^(length(x)-1) < alpha_min
    xx[!alphaOK]<-0
    xx[alphaOK] <- x[alphaOK]
  }
  return(xx)
}

# input: a sparse matrix with named rows and columns (dimnames)
# returns: a data frame representing triplets (r, c, x) suitable for writing to a CSV file
sparse2triples <- function(m) {
  SM = summary(m)
  D1 = m@Dimnames[[1]][SM[,1]]
  D2 = m@Dimnames[[2]][SM[,2]]
  data.frame(row=D1, col=D2, x=m@x)
}

sparse.cor4 <- function(x,y=NULL){
  n <- nrow(x)
  if(is.null(y)){
    cMeans <- colMeans(x)
    covmat <- (as.matrix(crossprod(x)) - n*tcrossprod(cMeans))/(n-1)
    sdvec <- sqrt(diag(covmat)) 
    cormat <- covmat/tcrossprod(sdvec)
  }else{
    xMeans <- colMeans(x)  
    yMeans <- colMeans(y)
    covmat<-((t(x)%*%y)- n*tcrossprod(xMeans,yMeans))/(n-1)
    
    xcovmat <- (as.matrix(crossprod(x)) - n*tcrossprod(xMeans))/(n-1)
    xsdvec <- sqrt(diag(xcovmat)) 
    ycovmat <- (as.matrix(crossprod(y)) - n*tcrossprod(yMeans))/(n-1)
    ysdvec <- sqrt(diag(ycovmat)) 
    
    cormat <- covmat/tcrossprod(xsdvec,ysdvec)
    cormat[is.nan(as.matrix(cormat))]<-0
  }
  list(cov=covmat,cor=cormat)
}

sparse.cor <- function(x,y=NULL){
  n <- nrow(x)
  cMeans <- colMeans(x)
  covmat <- (as.matrix(crossprod(x)) - n*tcrossprod(cMeans))/(n-1)
  sdvec <- sqrt(diag(covmat)) 
  cormat <- covmat/tcrossprod(sdvec)
  list(cov=covmat,cor=cormat)
}

#Calcula rowMeans en sparse marix sin contar los ceros
rowMeans_drop0 <- function (dgCMat) {
  RowInd <- dgCMat@i + 1
  sapply(split(dgCMat@x, RowInd), mean)
}

#idem con colMeans
colMeans_drop0 <- function (dgCMat) {
  nnz_per_col <- diff(dgCMat@p)
  ColInd <- rep.int(1:ncol(dgCMat), nnz_per_col)
  sapply(split(dgCMat@x, ColInd), mean)
}

#Calcula rowMSum en sparse marix sin contar los ceros
rowSums_drop0 <- function (dgCMat) {
  RowInd <- dgCMat@i + 1
  sapply(split(dgCMat@x, RowInd), sum)
}

#idem con colSums
colSums_drop0 <- function (dgCMat) {
  nnz_per_col <- diff(dgCMat@p)
  ColInd <- rep.int(1:ncol(dgCMat), nnz_per_col)
  sapply(split(dgCMat@x, ColInd), sum)
}

# Toma una lista de genes como simbolos (nombres) y devuelve sus egIDs, o viceversa
# mapIds va de key a columna
# NOTA: algunos elementos de la lista (mt-Co1, mt-Cytb, mt-Nd1) tienen el simbolo oficial pero 
# no aparecen asi en la base de datos, ni como alias. Otros elementos (3110021N24Rik) se desconocen
sym2eg <- function(genes,inverse=FALSE,keytype=c("ALIAS","ENSEMBL")[1]){
  require(org.Mm.eg.db)
  
  if (inverse){
    genes_t <- mapIds(org.Mm.eg.db, keys=genes, column=c("ENSEMBL"), keytype="ENTREZID")
  }else{
    genes_t <- mapIds(org.Mm.eg.db, keys=genes, column=c("ENTREZID"), keytype="ENSEMBL")
    
    # algunos nombres son secuencias:
    genes_na    <- names(genes_t[is.na(genes_t)])
    genes_na_id <- try(mapIds(org.Mm.eg.db, keys=genes_na, column=c("ENTREZID"), keytype="ACCNUM"),silent=TRUE)
    genes_t[names(genes_na_id)] <- genes_na_id
    
    #paste0(substring("FRMPD2",1,1),tolower(substring("FRMPD2",2)))
  }
  return(genes_t)
}


## PCA ####

# Devuelve una lista ordenada de nombres de genes segun su peso sobre los 
# componentes principales ncomp de X
# Inputs:
#  sce: singleCellExperiment con PCA (reducedDim(sca,"PCA")) definido
#  ncomp: numeric or vector. Si es un entero, indica la componente de interes. 
#         Si es un vector de 'm' componentes se consideran en conjunto las 'm' componentes
#         de manera agregada
#  retw: Si TRUE devuelve el peso ("importancia") de cada gen
genes_by_weight <- function(sce,ncomp=1,retw=FALSE){
  
  if(!"PCA"%in%reducedDimNames(sce)){
    stop("PCA reducedDim not found.\n")
  }
  
  p <- attr(reducedDim(sce,"PCA"),"rotation")[,ncomp]

    # lambdas <- (results$sdev)^2
  # w <- lambdas[ncomp]/sum(lambdas)
  w <- attr(reducedDim(sce,"PCA"),"percentVar")[ncomp]
  
  if (is.null(dim(p))){ # si solo hay un componente principal
    p <- p*w
    weights <- abs(p)
  }else{ # si hay mas de uno, es matriz
    p <- apply(p,1,function(x){x*w})
    weights <- apply(p,2,function(x){sum(abs(x))})
  }
  weights <- weights/sum(weights)
  # los pesos pueden variar entre prcomp y pca_svd, pero el orden es el mismo
  
  indices <- order(weights,decreasing = TRUE) # orden de maximas proyecciones
  
  if (retw){
    nombre <- weights[indices]
  }else{
    nombre <- rownames(results$rotation)[indices]
  }
  return(nombre)
}


# Test para identificar relaciones efectivas gen<->PC
# Inputs:
#  sce: singleCellExperiment con PCA (reducedDim(sca,"PCA")) definido 
#  alpha_min: umbral de significacncia para el filtro de disparidad (Serrano2009PNAS)
#  Si byPC=TRUE : para cada PC se considera la distribucion de pesos de los genes sobre
#                 la misma (columna de la matriz de rotacion) y se filtran contribuciones
#                 que no superan el filtro de disparidad
#     byPC=FALSE: para cada gen (fila de matriz de rotacion) se filtran contribuciones a PC's
#                 que no superen el filtro de disparidad
# Output:
#  Matriz de rotacion con valores que no pasan el filtro puestos a cero
getPCgenes <-function(sce, alpha_min = 0.05, byPC=TRUE){
  rot <-  attr(reducedDim(sce,"PCA"),"rotation") 

  disp<-t(apply(rot,2,function(x){
  xx <- abs(x)/sum(abs(x))
  alphaOK <- (1-xx)^(nrow(rot)-1) < alpha_min
  xx[!alphaOK]<-0
  xx[alphaOK] <- x[alphaOK]
  return(xx)
 }))
  return(t(disp))
}







## singleCellExperiment companion ####

# Calculate feature (row) and cell (column) QC metrics: total_counts, mean_counts
# n_features, pct_dropouts, mol_per_exp
# Inputs:
#  scexp          : a SingleCellExperiment object
#  featureCountMin: threshold for a feature to be considered expressed
# Output:
#   scexp object with updated rowData  and colData
QCmetrics <- function(scexp,featureCountMin=1){
  
  if("mol_per_exp"%in%colnames(colData(scexp))){
    cat("Exiting becasue the singleCellExpression onjects seems to already be QCed.\n")
    return()
  }
  
  mcounts <- counts(scexp)
  fqc<-apply(mcounts,1,function(x){
    total_counts <- sum(x)
    mean_counts  <- mean(x)
    n_cells      <- sum(x>=featureCountMin)
    pct_dropouts <- (length(x)-n_cells)/length(x)
    mol_per_exp  <- total_counts/n_cells
    return(c(total_counts=total_counts,pct_dropouts=pct_dropouts,mean_counts=mean_counts,mol_per_exp=mol_per_exp))
  })
  
  cqc<-apply(mcounts,2,function(x){
    total_counts <- sum(x)
    mean_counts  <- mean(x)
    n_features   <- sum(x>=featureCountMin)
    pct_dropouts <- (length(x)-n_features)/length(x)
    mol_per_exp  <- total_counts/n_features
    return(c(total_counts=total_counts,pct_dropouts=pct_dropouts,mean_counts=mean_counts,mol_per_exp=mol_per_exp))
  })
  
  rowData(scexp) <- cbind(rowData(scexp),DataFrame(t(fqc))) 
  colData(scexp) <- cbind(colData(scexp),DataFrame(t(cqc)))
  return(scexp)
}


## Networks ####

#Devuelve una lista con los primeros vecinos mutuos de cada columna de x
#Para la cuenta se fija dentro de los k primeros vecinos de cada columna y considera
#similaridad tipo pearson o euclidea
# Inputs:
#  xx: matriz de  expresion (i.e. logCounts)
#    k: nro de vecinos usados para identificar vecinos mutuos
#    mode: pearson/euclidea Metrica a considerar
# Outputs:
#  lista de primeros vecinos
getMKNNs <- function(xx,k,mode=c("pearson","euclidean")){
  if(mode=="pearson"){
    if(class(xx)[1]=="dgCMatrix"){
      d100 <- 1- (sparse.cor4(t(xx))$cor+1)/2
    }else{
      d100 <- 1- (cor(t(xx))+1)/2
    }
  }else{
    d100 <- as.matrix(dist(sx))
  }
  
  # KNN
  knn <- apply(d100,1,function(x){
    saux <- names(sort(rank(x),decreasing=FALSE))[2:(k+1)]
    return(saux) #devuelve el nombre de k_primeros_vecinos
  })
  
  #mutualKNN
  lmknn<-list()
  for(i in seq_along(knn[1,])){
    vecinos <- knn[,i]
    lmknn[[colnames(knn)[i]]] <- vecinos[unlist(apply(knn[,vecinos],2,function(x){colnames(knn)[i]%in%x}))]
  }    
  
  return(lmknn)
}

# Funcion para crear redes de K primeros vecinos mutuos
# Se pude utilizar el espacio de logcounts o PCA. En este ultmio caso se puede 
# restringir el nro de componentes a considerar, respecto al total disponible.
# En el espacio elegido es posible utilizar una distancia Euclidea o Correlacion
# El resultado es un grafo, cuyos enlaces tienen como atributo distintas medidas de similaridad
# Inputs:
#   sce: singleCellExperiment.
#   mutualK: nro de vecinos usados para identificar vecinos mutuos
#   use.reducedDim: si NA, se considera logcounts(sce)
#                   si PCA. se considera reducedDim(sce,"PCA")
#   num.PCA: numero de componentes a consdierar para el armado de la red
#   mode: metrica para cuantificar similaridad
# Output:
#   g: grafo igraph
buildMKNN     <-function(sce,mutualK=30,use.reducedDim=c(NA,"PCA")[2],num.PCA=NA,mode=c("pearson","euclidean")){
  if(is.na(use.reducedDim)){
    x <- t(logcounts(sce))
  }else{
    if(!use.reducedDim%in%reducedDimNames(sce)){
      stop("Non recognized reducedDim.\n")
    }else{
      x <- reducedDim(sce,use.reducedDim)
      if(!is.na(num.PCA)){
        x <- x[,1:(min(num.PCA,ncol(x)))]
      }
    }  
  }
  
  if(mode=="pearson"){
    if(class(x)[1]=="dgCMatrix"){
      d100 <- 1- (sparse.cor4(t(x))$cor+1)/2
    }else{
      d100 <- 1- (cor(t(x))+1)/2
    }
  }else{
    #plot(pca$sdev);abline(v=numPCA)
    d100 <- as.matrix(dist(x))
  }
  
  if(is.null(rownames(d100))){
    colnames(d100)<-rownames(d100)<-colData(sce)$cell_id
  }
  
  # KNN
  knn <- apply(d100,1,function(x){
    saux <- names(sort(rank(x),decreasing=FALSE))[2:(mutualK+1)]
    return(saux) #devuelve el nombre de k_primeros_vecinos
  })
  
  #mutualKNN
  lmknn<-list()
  for(i in seq_along(knn[1,])){
    vecinos <- knn[,i]
    lmknn[[colnames(knn)[i]]] <- vecinos[unlist(apply(knn[,vecinos],2,function(x){colnames(knn)[i]%in%x}))]
  }
  
  #Armo edgelist
  el <- c()
  for(i in seq_along(lmknn)){
    a<-lmknn[[i]]
    el<-rbind(el,cbind(rep(names(lmknn)[i],length(a)),a))
  }
  a<-t(apply(el,1,sort))
  el<-a[!duplicated(a),]
  
  g <- igraph::graph_from_edgelist(el,directed=FALSE)
  g <- igraph::simplify(g)
  
  simCellJ <- igraph::similarity(g,method="jaccard")
  simCellD <- igraph::similarity(g,method="dice")
  simCellW <- igraph::similarity(g,method = "invlogweighted")
  colnames(simCellJ)<-colnames(simCellD)<-colnames(simCellW)<-igraph::V(g)$name
  rownames(simCellJ)<-rownames(simCellD)<-rownames(simCellW)<-igraph::V(g)$name
  
  #para asignar las ismilitudes a los enlaces
  eattr<-t(apply(el,1,function(x){
    c(x,simCellJ[x[1],x[2]],simCellD[x[1],x[2]],simCellW[x[1],x[2]])
  }))    
  colnames(eattr) <- c("Source","Target","simJ","simD","simW")
  rownames(eattr)  <- apply(eattr[,1:2],1,function(x){paste0(sort(x),collapse=":")})
  
  #pongo las similutdes calculadas como atrubuto de edges del grafo
  saux <- apply(igraph::ends(g,igraph::E(g)),1,function(x){paste0(sort(x),collapse=":")})
  igraph::E(g)$simJ <- as.numeric(eattr[saux,"simJ"])
  igraph::E(g)$simD <- as.numeric(eattr[saux,"simD"])
  igraph::E(g)$simW <- as.numeric(eattr[saux,"simW"])
  
  g <- igraph::set.graph.attribute(g,"params",paste0("K",mutualK,"_reduced",use.reducedDim,"_Dim",ncol(x),"_mode",mode))
  return(g) 
}


# Calculo de clusterizacion sobre el grafo g. El resultado queda almacenado en colData(sce)
# Inputs:
#   sce: single cell object (usado para calcular previamente el grafo g)
#   g: grafo igraph
#   lAlgorithm: es una lista de funciones igraph de deteccion de comunidades
#       lAlgorithms <- list(louvain=igraph::cluster_louvain, 
#                           infomap=igraph::cluster_infomap,
#                           labProp=igraph::cluster_label_prop,
#                           fgreedy=igraph::cluster_fast_greedy,
#                           walktrap=igraph::cluster_walktrap)
#   prefix: prefijo para armar e nombre de las columnas de colData(sce)
# Outpus:
#   singleCellExperiment object con resultados de clusterizacion en colData(sce)
doClustering<-function(sce,g,lAlgorithms,prefix="label_"){
  lclus<-list()
  #layout(matrix(seq_along(lAlgorithms),2,2))
  for(i in seq_along(lAlgorithms)){
    lclus[[i]]<-lAlgorithms[[i]](g)
    #    (tt<-table(colData(sce)[,"cell_type"],membership(lclus[[i]])[rownames(colData(sce))]))
    #    (inf<-infoPartition(tt))
    #    plot(inf,as.numeric(table(lclus[[i]]$membership)),xlim=c(0,1),xlab="coherence",ylab="size",log="y",main=names(lAlgorithms)[i])
  }
  names(lclus)<-names(lAlgorithms)
  
  laux <- lapply(lclus,igraph::membership)
  aux <- matrix(unlist(laux),ncol=length(lclus))
  colnames(aux)<-names(lclus)
  rownames(aux)<-names(laux[[1]])
  
  labels          <- matrix(0,ncol=length(lclus),nrow=ncol(sce))
  rownames(labels)<- colnames(sce)
  labels[rownames(aux),]<-aux
  colnames(labels)<-paste0(prefix,colnames(aux))
  labels<-apply(labels,2,function(x){   #pongo los ceros como clusters de un elemento
    i0     <- which(x==0)
    newVal <- max(x)+seq_along(i0)
    x[i0]  <- newVal
    return(x)
  })  
  
  colData(sce) <- cbind(colData(sce),labels)
  cat("Clustering done\n")
  return(sce)
}

# Devueve data.frame con las celulas centrales de los clusters 
# cuyas etiquetas estan almacenadas en a columna 'labelcolname' de colData(sce)
getCentroidCells<-function(sce,gg,labelcolname=NULL){
  if(is.null(labelcolname)) warning("You should specified labelcolname.\n")
  if(!labelcolname%in%colnames(colData(sce)))warning(paste(labelcolname),"not found in colData(sce).\n")
  labels   <- colData(sce)[,labelcolname]
  ulabels  <- unique(labels) 
  res<-c()
  for(i in seq_along(ulabels)){
    inodes <- rownames(colData(sce))[which(labels%in%ulabels[i])]
    if(length(inodes)==1){
      cellCentroid<- inodes  
      ssize<-1
    }else{
      subg   <- igraph::induced.subgraph(gg,vids=inodes)
      cellCentroid<-names(which.max(igraph::coreness(subg)+igraph::degree(subg)+igraph::closeness(subg)))[1]
      ssize <- igraph::vcount(subg)
    }
    res<-rbind(res,data.frame(label=ulabels[i],cellID=cellCentroid,clusterSize=ssize))
    
  }
  return(res)
}
 
# (ex grafo-comms de Luz)
# Funcion que toma un grafo y una particion, y genera un nuevo grafo
# donde cada nodo es una comunidad, con enlaces pesados no dirigidos
# input: 
#  g   : igraph graph
#  memb: named membership vector
#  wlabel: name of the igraph edge attribute with weights to be used to merge clusters
#  factorIntra: factor to downweight intra-cluster weights
getClustersGraph <-function(g, memb, wlabel=NULL,factorIntra=1){
  
  if (vcount(g)!=length(memb)){
    interseccion <- intersect(names(memb), vertex_attr(g, name='name'))
    stopifnot(length(interseccion)>0)
    g <- induced_subgraph(g, interseccion)
    memb <- memb[interseccion]
    print('Se utilizo el subconjunto interseccion entre los nodos del grafo y los elementos de la particion')
  }
  
  enlaces <- as_edgelist(g)
  ecoms   <- apply(enlaces, 2, function(x){memb[x]})
  G           <- graph_from_edgelist(ecoms, directed = FALSE)
  if(!is.null(wlabel)){
    E(G)$weight <- get.edge.attribute(g,wlabel)
  }else{
    E(G)$weight <- 1
  }
  G           <- igraph::simplify(G, remove.loops=FALSE, edge.attr.comb=list(weight="sum"))
  
  
  eclusters <- get.edgelist(G)
  isame     <- which(apply(eclusters,1,function(x){x[1]==x[2]}))
  G <- set.edge.attribute(G,"weight",isame,get.edge.attribute(G,"weight",isame) * factorIntra)  
  
  
  return(G)
}



#Exporta el grafo gg a formato Gephi. colData de sce se exporta como atributos de nodos 
# Inputs:
#   sce: singleCellExperiment
#   gg : igraph network to be exported
#   file.name: nombre para el file a exportar
exportGraphToGephi<-function(sce,gg,file.name="foo"){
  aux  <- edge_attr(gg)
  maux <- matrix(unlist(aux),ncol=length(aux))
  colnames(maux)<-names(aux)
  eattr    <- data.frame(ends(gg,E(gg)),maux)
  colnames(eattr)<-c("Source","Target",names(aux))
  
  #edges 
  write.table(eattr,sep=",",row.names=FALSE,col.names=TRUE,
              file=paste0(file.name,".adjList.csv"),quote=FALSE)
  
  el <- get.edgelist(gg)
  el <- t(apply(el,1,sort))
  
  #node attributes
  maux<-colData(sce)
  linkedNodes <- unique(c(el[,1],el[,2]))
  a<-maux[,1]%in%linkedNodes
  if(all(!a)) a<-rownames(maux)%in%linkedNodes
  maux<-maux[a,]
  colnames(maux)[1]<-"Id"
  aux <- grep("percent_top",colnames(maux))
  if(length(aux)>0) maux <- maux[,-aux]
  
  #pongo los atributos de nodo en el grafo
  for(i in 2:(ncol(maux))){
    if(class(maux[1,i])=="character"){
      gg <- set.vertex.attribute(gg,colnames(maux)[i],V(gg),maux[,i])
    }else{
      gg <- set.vertex.attribute(gg,colnames(maux)[i],V(gg),as.numeric(maux[,i]))
    }
  }
  
  write.table(maux,sep=",",row.names=FALSE,col.names=TRUE,file=paste0(file.name,".nodes.csv"),quote=FALSE)
  
  cat("Geaph exported.\n")
}

communicabilityDistance<-function(gr,cells,sce){
  # require(brainGraph)
  # comps<-clusters(gr)
  # for(icomp in 1:comps$no){
  #   ccells <- names(comps$membership)[comps$membership==icomp]
  #   ccomm <- communicability(induced_subgraph(gr,ccells)) 
  #   colnames(ccomm)<-rownames(ccomm)<-ccells
  #   a <- ccomm[intersect(ccells,cells),intersect(ccells,cells)]
  #   for(i in 1:nrow(a)){
  #     for(j in seq(i+1nrow(a)-1L))
  #   }
  # }
}

#Funcion para chequear el contenido de informacion de una clusterizacion
#cuando se la compara con la curacion manual de tipo de celula
infoPartition <- function(tt){
  apply(tt,2,function(x){
    x<-x/sum(x)
    s <- -sum(x[x>0]*log2(x[x>0]))
    s <- 1-s/log2(length(x))
  })
}



## Functional analysis ####
# Estandariza el valor de logcounts de cada gen. 
# Para lidiar con dropouts se consideran medias sobre regiones locales,
#   1)knn.reducedDim : nro de primeros vecinos para buscar primeros vecinos mutuos en reducedDim
#            mode.reducedDim: similaridad tipo pearson o euclidea
#   2)g  :grafo para determinar localia a partir de 1eros vecinos 
#   3)groups: le paso los grupos de localia (en este caso tomo valor y conisdero el resto para mu y sd)
# Inputs:
#     sce: singleCellExperiment
# si bstandarize=TRUE se devuelve el foldchange estandarizado por gene

relativizeExpression<-function(sce,g=NULL,bStandarize=FALSE, remove.zeros=TRUE,
                               use.reducedDim=NULL,reducedDim.d=50,
                               reducedDim.knn=20,reducedDim.mode=c("pearson","euclidean")[1],
                               groups=NULL ){
  
  if(!"logcounts"%in%assayNames(sce)){
    stop("No logcounts assay found for the provided singleCellExperiment object.\n")
  }
  x <- logcounts(sce)
  
  llocal<-NULL
  if(!is.null(use.reducedDim)){
    xx <- reducedDim(sce,use.reducedDim)
    if(!is.null(reducedDim.d)) xx <- xx[,1:min(reducedDim.d,ncol(xx))]
    llocal <- getMKNNs(xx,k=reducedDim.knn,mode=reducedDim.mode)
  }else{
    if(!is.null(g)){
      llocal<-vector(mode="list",length=ncol(x))
      names(llocal)<-colnames(x)
      
      llocal<-neighborhood(g)
      llocal<-lapply(llocal,function(x){x$name})
      names(llocal)<-unlist(lapply(llocal,function(x){x[1]}))
      addCells   <- colnames(x)[!colnames(x)%in%names(llocal)]
      laux       <- as.list(addCells) #vector(mode="list",length=length(addCells))
      names(laux)<- addCells
      llocal <- c(llocal,laux)
    }
  }
  
  #si vengo de use.reducedDim o de graph paso por aca...
  if(!is.null(llocal)){ 
    # uso media de la localia como x y estandarizo con la media y sd del complemento
    Z <- matrix(0,ncol=length(llocal),nrow=nrow(x))
    colnames(Z)<-names(llocal)
    rownames(Z)<-rownames(x)
    
    for(icell in seq_along(llocal)){
      cat("\r", sprintf("%.1f", 100*icell/length(llocal)),"%")
      coi <- llocal[[icell]]  #cells of interest
      
      a<-x[,coi,drop=FALSE]
      if(remove.zeros){
        aux     <- rowMeans_drop0(a)
        mu_cell <- mu <- mu2 <- rep(0,length=nrow(a))   
        names(mu_cell)<-names(mu)<-names(mu2)<-rownames(a)
        mu_cell[as.numeric(names(aux))]<-aux
        
        aux       <- rowMeans_drop0(x[,!colnames(x)%in%coi])
        mu[as.numeric(names(aux))]<-aux
        if(bStandarize){
          aux       <- rowMeans_drop0(x[,!colnames(x)%in%coi]*x[,!colnames(x)%in%coi])
          mu2[as.numeric(names(aux))]<-aux
        }
      }else{
        mu_cell   <- rowMeans(a)
        mu        <- rowMeans(x[,!colnames(x)%in%coi])
        if(bStandarize) mu2       <- rowMeans(x[,!colnames(x)%in%coi]*x[,!colnames(x)%in%coi])
      } 
      Z[,icell]  <- (mu_cell-mu)
      
      if(bStandarize){
        sd <- sqrt(mu2-mu*mu)
        Z[,icell] <- Z[,icell]/sd 
      }
      if(FALSE){
        mu_cell_0 <- rowMeans_drop0(a)
        
        #como es el nivel de dropout local
        rowInd    <- a@i + 1
        localDropOut     <- 1-table(rowInd)/length(coi)
        completeDO       <- rep(1,length=length(mu_cell)-length(localDropOut))
        names(completeDO)<- setdiff(as.character(1:length(mu_cell)),names(localDropOut))
        localDropOut <- c(localDropOut,completeDO)
        hist(localDropOut,main="local dropout distribution")
        
        plot(mu_cell[as.numeric(names(mu_cell_0))],mu_cell_0)
        abline(0,1)
      }
      
    }  
    
  }else{  
    if(remove.zeros){
      mu  <- rowMeans_drop0(x)  
      if(bStandarize) mu2 <- rowMeans_drop0(x*x)
    }else{
      mu  <- rowMeans(x)
      if(bStandarize) mu2 <- rowMeans(x*x)
    }
    Z <- (x - mu)
    if(bStandarize){
      sd<- sqrt(mu2 - mu*mu)
      Z<-Z/sd
    }
  }
  return(Z)
}

# Estandariza el valor de logcounts de cada gen. Lo hace de manera eficiente
# utilizando multiplicacion d ematrices esparsas.
# Para lidiar con dropouts se consideran medias sobre regiones 
# locales definidas a partir de primeros vecinos. 
# Se consideran siempre las contribuciones no nulas en los promedios
# Luego se centran los valores obtenidos a partir de las medias
# (calculadas sin los ceros) sobre los no-vecinos
# Finalmente se estandariza estimando la varianza de los no vecinos
# (nuevamente ignorando contribucioes nulas)
# Inputs:
#     sce: singleCellExperiment
#      g : grafo para estimar vecindades y localias
# Output:
# lista con matrices: mu, centered y Z
standarizeCellExpression<-function(sce,g){ 
  k        <- degree(g)
  neigbs   <- get.adjacency(g,sparse = TRUE)
  Matrix::diag(neigbs) <- 1
  
  noneigbs <- 1-neigbs #get.adjacency(myg40,sparse = FALSE)
  Matrix::diag(noneigbs) <- 0
  
  Ct   <- t(logcounts(sce)[,colnames(neigbs)])
  
  #calculo la media del perfil transcrpcional de cada celula, ignorando contribuciones nulas
  ACt  <- neigbs %*% Ct
  k0   <- neigbs %*% (Ct>0)  #numero de celulas vecinas con contribuciones no-nulas
  mu   <- ACt/k0             #mu es un promedio, ignorando contribuciones nulas
  
  if(FALSE){
    plot(logcounts(sce)[,rownames(mu)[11]],mu[11,],xlab='cpm',ylab='neigb.mean')
    abline(0,1)
  }
  
  #calculo la media del perfil transcripcional fuera del entorno de cada celula
  NNCt  <- noneigbs %*% Ct     # computation time (!)
  k0    <- noneigbs %*% (Ct>0) # computation time (!)
  muOut  <- NNCt/k0
  muOut2 <- NNCt^2/k0
  mySD   <- sqrt(muOut2 - muOut^2) 
  
  Z1 <- mu - muOut
  Z2 <- Z1/mySD
  return(list(mu=t(mu),centered=t(Z1),Z=t(Z2)))
}


# Calcula la media de primeros vecinos para imputar el perfil de cada celula
# La nocion de primeros vecinos se puede obtener: de un grafo, o de calcular internamente
# un KMNN en el espacio reduceddim, con la metrica indicada
#  g: igraph object usded to identify nearest neighbors
#  remove.zeros: logical. If TRUE, disregard zero values in mean estimation
#  imputeOnly0: Only inpute for zero-count genes
#  use.reducedDim: (only if g==NULL) specify the reducedDim name to consider
#   reducedDim.d  : how many components to consider.Makes sense only for ordered components
#   reducedDim.knn: how many first neighbors to consider in reducedDim space.
#   reducedDim.mode: metric in rducedDim space
#  groups: explicitely defined groups for imputation(NOT IMPLEMENTED YET)  
myImputation<-function(sce,g=NULL,remove.zeros=TRUE, imputeOnly0 = TRUE,
                       use.reducedDim=NULL,reducedDim.d=50,
                       reducedDim.knn=20,reducedDim.mode=c("pearson","euclidean")[1],
                       groups=NULL ){
  
  if(!"logcounts"%in%assayNames(sce)){
    stop("No logcounts assay found for the provided singleCellExperiment object.\n")
  }
  x <- logcounts(sce)
  
  llocal<-NULL
  if(!is.null(use.reducedDim)){
    xx <- reducedDim(sce,use.reducedDim)
    if(!is.null(reducedDim.d)) xx <- xx[,1:min(reducedDim.d,ncol(xx))]
    llocal <- getMKNNs(xx,k=reducedDim.knn,mode=reducedDim.mode)
    slocal<-names(llocal)
    llocal <- lapply(slocal,function(x){c(x,llocal[[x]])})
    names(llocal)<-slocal
    cat(paste("Using", use.reducedDim,"neighborhoods.\n"))
  }else{
    if(!is.null(g)){
      llocal<-vector(mode="list",length=ncol(x))
      names(llocal)<-colnames(x)
      
      llocal<-neighborhood(g)
      llocal<-lapply(llocal,function(x){x$name})
      names(llocal)<-unlist(lapply(llocal,function(x){x[1]}))
      addCells   <- colnames(x)[!colnames(x)%in%names(llocal)]
      laux       <- as.list(addCells) #vector(mode="list",length=length(addCells))
      names(laux)<- addCells
      llocal <- c(llocal,laux)
      llocal <- llocal[colnames(x)]
      cat(paste("Using graph neighborhoods.\n"))
      
    }
  }
  
  #si vengo de use.reducedDim o de graph paso por aca...
  if(!is.null(llocal)){ 
    # uso media de la localia como x y estandarizo con la media y sd del complemento
    # Z <- matrix(0,ncol=length(llocal),nrow=nrow(x))
    # colnames(Z)<-names(llocal)
    # rownames(Z)<-rownames(x)
    Z <- as.matrix(logcounts(sce))
        
    for(icell in seq_along(llocal)){
      cat("\r", sprintf("%.1f", 100*icell/length(llocal)),"%")
      coi <- llocal[[icell]]  #cells of interest
      if(length(coi)==0){
        mu_cell <- x[,names(llocal)[icell]]     
      }else{
        a<-x[,coi,drop=FALSE]
        if(remove.zeros){
          aux     <- rowMeans_drop0(a)
          mu_cell <- rep(0,length=nrow(a))   
          names(mu_cell)<-rownames(a)
          mu_cell[as.numeric(names(aux))]<-aux
        }else{
          mu_cell   <- rowMeans(a)
        }
      }
      if(imputeOnly0){
        i0<-logcounts(sce)[,icell]==0
        Z[i0,icell]  <- mu_cell[i0]
      }else{
        Z[,icell]  <- mu_cell
      }
      
    }  
    
  }else{
    warning("No graph, nor reducedDim space specified.\n")
    return(NULL)
  }
  return(Z[,colnames(sce)])
}


# fgsea devuelve tabla de 8 variables: pathway, pval, padj, ES, NES, nMoreExtreme, size, leadingEdge
enrichment.old <- function(standarized=NULL,ssce=NULL,functional_groups=NULL,scoreType=c("std","pos","neg")[1],
                       bparallel=FALSE,ncores=20){
  
  if(is.null(standarized)){
    if("standarized"%in%assayNames(ssce)){
      standarized <- assay(ssce,"standarized")
      nameCols <- colnames(ssce)
    }else{
      stop("No standarized data provided, nor standarized assay available.")
    }
  }
  nameCols <- colnames(standarized)
  Ncells <- ncol(standarized) 
  lenrichment<-list()
  if(bparallel){
    require(foreach)
    require(doParallel)
    registerDoParallel(cores=ncores)

    lenrichment<-foreach(i=1:Ncells) %dopar%{
      gsea_res <- fgseaMultilevel(pathways = functional_groups, stats = standarized[,i],scoreType = scoreType)
      gsea_res
    }
    
  }else{
  for (i in 1:Ncells){
    lenrichment[[i]]<-list()
    cat("\r", sprintf("%.1f", 100*i/Ncells), "%")
    gsea_res <- fgseaMultilevel(pathways = functional_groups, stats = standarized[,i],scoreType = scoreType)
    lenrichment[[i]] <- gsea_res
    if(FALSE) {
      plotEnrichment(functional_groups[[2]],stats=standarized[,i])
      plotGseaTable(functional_groups,stats=standarized[,i],gsea_res)
    }
  }
  }
  names(lenrichment)<-nameCols
  return(lenrichment)
}

enrichment.fgseaML <- function(standarized=NULL,ssce=NULL,functional_groups=NULL,scoreType=c("std","pos","neg")[1],
                                   ncores=1){
  if(is.null(standarized)){
    if("standarized"%in%assayNames(ssce)){
      standarized <- assay(ssce,"standarized")
      nameCols <- colnames(ssce)
    }else{
      stop("No standarized data provided, nor standarized assay available.")
    }
  
  }
  nameCols <- colnames(standarized)
  Ncells <- ncol(standarized) 
  lenrichment<-list()
  
  for (i in 1:Ncells){
    lenrichment[[i]]<-list()
    cat("\r", sprintf("%.1f", 100*i/Ncells), "%")
    aux <- standarized[,i]
    aux <- aux[!is.nan(aux)]

    gsea_res <- fgseaMultilevel(pathways = functional_groups, stats = aux,scoreType = scoreType,
                                 nproc = ncores)
       lenrichment[[i]] <- gsea_res
    if(FALSE) {
      plotEnrichment(functional_groups[[2]],stats=standarized[,i])
      plotGseaTable(functional_groups,stats=standarized[,i],gsea_res)
    }
  }
  
  names(lenrichment)<-nameCols
  return(lenrichment)
}

enrichment.fgseaSimple <- function(standarized=NULL,ssce=NULL,functional_groups=NULL,scoreType=c("std","pos","neg")[1],
                           ncores=1){
  
  if(is.null(standarized)){
    if("standarized"%in%assayNames(ssce)){
      standarized <- assay(ssce,"standarized")
      nameCols <- colnames(ssce)
    }else{
      stop("No standarized data provided, nor standarized assay available.")
    }
    
  }
  nameCols <- colnames(standarized)
  Ncells <- ncol(standarized) 
  lenrichment<-list()

    for (i in 1:Ncells){
      lenrichment[[i]]<-list()
      cat("\r", sprintf("%.1f", 100*i/Ncells), "%")
      # gsea_res <- fgseaMultilevel(pathways = functional_groups, stats = standarized[,i],scoreType = scoreType,
      #                             nproc = ncores)
      aux <- standarized[,i]
      aux <- aux[!is.nan(aux)]
      gsea_res <- fgseaSimple(pathways = functional_groups, 
                              stats    = aux,
                              nperm    = 1000,
                              scoreType = scoreType,
                              nproc = ncores)
      
      lenrichment[[i]] <- gsea_res
      if(FALSE) {
        plotEnrichment(functional_groups[[2]],stats=standarized[,i])
        plotGseaTable(functional_groups,stats=standarized[,i],gsea_res)
      }
    }
  
  names(lenrichment)<-nameCols
  return(lenrichment)
}

#Calculo de enriquecimiento considerando una red Bipartita GO-Genes
# Para cada GO se acumula el score de cada gen (expresion, Zvalue, o Zcentered)
# y se compara con 1000 recableado de la relacion GO<->gene
# devuelve un valor z para cada categoria ((z - mu_0)/sd_0)
# TODO: agregar para que duevla p-value
enrichment.ZGO <- function(inputGOs, sceZ, org="Mm", pvmin=0.05, seed =  123457, 
                           nRND1 = 1000,
                           evGO=c("EXP", "IDA", "IPI", "IMP", "IGI", "IEP", 
                                  "HTP", "HDA", "HMP", "HGI", "HEP")){
  
  db <- paste0("org.",org,".eg.db")
  
  # . Prepraro categorias GO
  require(db)
  #inputGOs <- enrichZpos[[1]]$pathway
  db1 <- paste0("org.",org,".egGO2ALLEGS")
  db2 <- paste0("org.",org,".egSYMBOL")
  a<-mget(unique(inputGOs),db1)
  goCats<-lapply(a,function(x){
    x<-x[names(x)%in%evGO] #mapeos con evidencia apropiada
    x<-unlist(mget(x,db2,ifnotfound = NA))
    x<-x[x%in%rownames(sceZ)]
    return(unique(x))
  })
  
  # . armo la matriz GO/Gene
  genesByGO<-lapply(names(a),function(x){
    cbind(rep(x,length(goCats[[x]])),goCats[[x]])
  })
  names(genesByGO) <- names(a)
  res<-c()
  for(i in seq_along(genesByGO)){
    res<-rbind(res,genesByGO[[i]])
  }
  ugenes <- unique(res[,2])
  i <- match(res[,1],names(genesByGO))
  j <- match(res[,2],ugenes)
  mGOgene<-sparseMatrix(i,j,1)
  colnames(mGOgene)<-ugenes
  rownames(mGOgene)<-names(genesByGO)
  
  observed <- mGOgene %*% sceZ[colnames(mGOgene),]
  
  # . rnd1: GO-gene
  set.seed(seed)
  t1 <- Sys.time()
  
  for(i in 1:nRND1){
    cat("\r", sprintf("%.1f", 100*i/nRND1), "%")
    
    rndGOgene <- t(apply(mGOgene,1,sample,ncol(mGOgene)))
    a <- rndGOgene %*% sceZ[colnames(mGOgene),]
    if(i==1){
      accum <- a
      accum2<- a^2
    }else{
      accum <- accum + a
      accum2<- accum2+ a^2
    }
  }
  zobs <- (observed - accum/nRND1)/sqrt(accum2/nRND1-(accum/nRND1)^2)
  t2 <- Sys.time()
  t2-t1
  return(zobs)
}


# getNESpartition:
#  Devuelve la particion con clusters definidos originalmente en 'lab'
#  mergeados segun intensidad de enlaces alacenadas en wlabel (tipicamente NES-GO)
#  siguiendo prescripcion Louvain. Las intensidades intra-clusters se reescalean segun
#  el factor 'fintra'
# gg   : igraph graph MKNN usado para ingferir particion 'lab'que se dessea 'mejorar'
# lab  : vector de particion (named)
# minSize: comunidades de tamanio menor a minSize no son analizadas
# fintra: (0,1] factor de reescaleo de las interacciones intraclusters
# wlabel: nombre del edgeAttribute donde se encuentran las similaridades utilizadas en 
#         el merging (en el sentido de Louvain)
getNESpartition <- function(gg,lab,fintra=1,minSize=50,wlabel = "simNES"){
  
  if(!wlabel%in%names(get.edge.attribute(gg))){
    stop(paste0("No", wlabel, "edge attributed detected.\n"))
  }
  
  #Comunidades de mas de minSize elementos
  tt     <-table(lab)
  clusOk <- names(tt[tt>minSize])
  
  # armo grafo de comunidades 
  g1 <- getClustersGraph(gg,memb = lab, wlabel = wlabel, factorIntra = fintra)
  
  # . louvain con pesos NES 
  louv <- cluster_louvain(g1,weights=V(g1)$weight)
  tt<-table(louv$membership)
  if(!any(tt>1)){
    warning("No refinement was achieved. Stop\n")
    return(NA)
  }
  
  # . reetiquetado 
  newlab <- louv$membership[lab]
  names(newlab) <- names(lab)
  #se pueden prodcir NAs porque los clusters originales de tamanio 1 no son tenidos en cuenta 
  #
  ina <- which(is.na(newlab))
  if(length(ina)>0){
    start <- max(newlab,na.rm=TRUE)+1
    newlab[ina] <- start:(start+length(ina)-1)
  }
  return(newlab)
}

#
# Z   : matrix de genes x celulas que contiene valor estandarizado de expresion
# lab : vector de particion (named)
# minSize: comunidades de tamanio menor a minSize no son analizadas
# nMarkers: numero de marcadores a identificar por comunidad (si disparity.pv=0)
# output:
#  lista: 
#    markers: nMarkers genes marcadores para cada comunidad 
#    meanMarkersInfo:  mean informativity score para los marcadores de c/comunidad 
#    infoGene: informatividad asociada a cada gen para la particion analizada (1-H)
require(doParallel)
getPartitionMarkers<-function(Z,lab,nMarkers=5,minSize=50,sceData=NULL,numCores=1){
  registerDoParallel(numCores)
  # . comienza el calculo de entropia 
  # para la cuenta de entropia solo analizo clusters de mas de 50 celulas
  tt <- table(lab)
  labOK <- tt > minSize
  labOK <- names(labOK)[labOK]
  
  #calcula entropio/informatividad de cada gen
  resS <- c()
  #for(igen in 1:nrow(Z)){
  aa<-foreach(igen=1:nrow(Z),.combine='c') %dopar% {
    a     <- aggregate(Z[igen,lab%in%labOK],by=list(lab[lab%in%labOK]),mean)
    x <- (a[,2] - min(a[,2]))/(sum(a[,2])-min(a[,2])*nrow(a))
    s <- -sum(x[x>0]*log2(x[x>0]))
    s <- s/log2(length(x))
    s.max <- a[which.max(x),1]
    resS <- c(resS,c(1-s,s.max))
  }
  
  mInfo <- matrix(aa,byrow=TRUE,ncol=2)
  colnames(mInfo) <- c("infoX","clusterMaxX")
  rownames(mInfo) <- rownames(Z)
  
  # . Identifico los top-nMarker mas informativos para cada cluster 
  markers <- by(mInfo,as.factor(mInfo[,"clusterMaxX"]),function(x){
    res <- rownames(x)[order(x$infoX,decreasing=TRUE)[1:nMarkers]]
    return(res[!is.na(res)])
  })
  
  mInfoAux          <- matrix(NA,ncol=2,nrow=nrow(mInfo))
  colnames(mInfoAux)<- c("rankInClust","disparity.pv")
  rownames(mInfoAux)<- rownames(mInfo)
  clusters <- unique(mInfo[,2])
  for(i in seq_along(clusters)){
    aux <- mInfo[mInfo[,2]==clusters[i],]
    disp.pv    <- disparity(aux[,1],alpha_min = 0) 
    rankInClus <- rank(-aux[,1])
    mInfoAux[names(disp.pv),"disparity.pv"]   <- disp.pv
    mInfoAux[names(rankInClus),"rankInClust"] <- rankInClus
  }
  mInfo <- cbind(mInfo,mInfoAux)
  
  if(!is.null(sceData)){
    if(sum(c("mean","detected")%in%colnames(rowData(sceData)))==2){
      aux   <- rowData(sceData)[rownames(mInfoAux),c("mean","detected")]
      mInfo <- data.frame(cbind(mInfo,aux))
    }
  }
  
  
  
  # . Calculo informatividad media de los 5 marcadores
  scores<-unlist(lapply(markers,function(x){
    mean(mInfo[x,"infoX"])
  }))
  
  return(list(markers=markers, meanMarkerInfo=scores, infoGene=mInfo))
}


rankStats <- function(sce,standarized=NULL,functional_groups=NULL){
  
  if(is.null(standarized)){
    if("standarized"%in%assayNames(sce)){
      standarized <- assay(sce,"standarized")
    }else{
      stop("No standarized data provided, nor standarized assay available.")
    }
  }
  
  lrank<-list()
  mrank <- apply(standarized,2,order,decreasing=TRUE)
  rownames(mrank)<-rownames(standarized)
  colnames(mrank)<-colnames(standarized)
  
  a<-lapply(functional_groups,function(fg){
    apply(mrank,2,function(y){
      return(c(min(y[fg]),mean(y[fg])))
    })})
  
  for (i in 1:Ncells){
    cat("\r", sprintf("%.1f", 100*i/Ncells), "%")
    
    gsea_res <- fgseaMultilevel(pathways = functional_groups, stats = standarized[,i])#,scoreType = "pos")
    lenrichment[[i]] <- gsea_res
    
  }  
}


## GO ----
computarIC_ach <- function(dbgo, keytype = "ENTREZID", ont, evCode=NA){
  x<-as.data.frame(get(dbgo))
  x<-x[x$Ontology%in%ont,]
  if(!is.na(evCode[1]))x<-x[x$Evidence%in%evCode,]
  idup <- duplicated(x[,1:2])
  x <- x[!idup,]
  gocounts <- table(x$go_id)
  if(FALSE){
    io<-order(gocounts,decreasing=TRUE)
    gocounts[io][1:4]
  }
  ic <- -log2(gocounts/max(gocounts))
  gos <- names(ic)
  ic  <- as.numeric(ic)
  names(ic)<-gos
  return(ic)  
}

godata_ach <- function(db, keytype = "ENTREZID", ont, evCode=NA){
  x<-as.data.frame(get(paste0(db,"GO2ALLEGS")))
  x<-x[x$Ontology%in%ont,]
  if(!is.na(evCode[1]))x<-x[x$Evidence%in%evCode,]
  idup <- duplicated(x[,1:2])
  x <- x[!idup,]
  gocounts <- table(x$go_id)
  if(FALSE){
    io<-order(gocounts,decreasing=TRUE)
    gocounts[io][1:4]
  }
  ic <- -log2(gocounts/max(gocounts))
  gos <- names(ic)
  ic  <- as.numeric(ic)
  names(ic)<-gos
  
  kk <- unique(x$gene_id)
  res <- new("GOSemSimDATA", keys = kk, ont = ont, IC=ic,geneAnno = x, 
             metadata = metadata(get(paste0(db,".db"))))
}

##' @importFrom AnnotationDbi as.list
##' @importFrom GO.db GOBPOFFSPRING
# fuente: GOSemSim https://rdrr.io/bioc/GOSemSim/src/R/computeIC.R
computarIC <- function(OrgDb, keytype = "ENTREZID", ont, evCode=NA) {
  ont <- toupper(ont)
  ont <- match.arg(ont, c("BP", "CC", "MF"))
  
  OrgDb <- load_OrgDb(OrgDb)
  kk    <- keys(OrgDb, keytype=keytype)
  goAnno <- suppressMessages(select(OrgDb, keys=kk, keytype=keytype,columns=c("GO", "ONTOLOGY")))
  
  goAnno <- goAnno[!is.na(goAnno$GO), ]
  goAnno <- goAnno[goAnno$ONTOLOGY == ont,]
  if(!is.na(evCode[1])) goAnno <- goAnno[goAnno$EVIDENCE%in%evCode,]

  goids <- select(GO.db, keys="BP", columns=c("GOID"), keytype="ONTOLOGY")
  goids <- goids$GOID
  ## all GO terms appearing in an given ontology ###
  goterms=goAnno$GO
  gocount <- table(goterms)
  ## goid of specific organism and selected category.
  goname  <- names(gocount) 
  
  ## ensure goterms not appearing in the specific annotation have 0 frequency
  go.diff        <- setdiff(goids, goname) #esta en goids pero no en goname
  m              <- double(length(go.diff))
  names(m)       <- go.diff
  gocount        <- as.vector(gocount)
  names(gocount) <- goname
  gocount        <- c(gocount, m)
  
  Offsprings <- switch(ont,
                       MF = AnnotationDbi::as.list(GOMFOFFSPRING),
                       BP = AnnotationDbi::as.list(GOBPOFFSPRING),
                       CC = AnnotationDbi::as.list(GOCCOFFSPRING))
  
  cnt <- gocount[goids] + sapply(goids, function(i) sum(gocount[Offsprings[[i]]], na.rm=TRUE))
  names(cnt) <- goids
  # CORRECCION DE ESTAS DOS ULTIMAS LINEAS:
  #cnt <- gocount + sapply(names(gocount), function(i) sum(gocount[Offsprings[[i]]], na.rm=TRUE))
  #names(cnt) <- names(gocount)
  
  ## the probabilities of occurrence of GO terms in a specific corpus.
  p <- cnt/sum(gocount)
  ## IC of GO terms was quantified as the negative log likelihood.
  IC <- -log(p)
  return(IC)
}


indecedTermGraph<-function (r, id, children = TRUE, parents = TRUE) 
{
  if (!children && !parents) 
    stop("children and parents can't both be FALSE")
  goOnt <- testName(r)[2]
  goKidsEnv <- GOenv(paste(goOnt, "CHILDREN", sep = ""))
  goParentsEnv <- GOenv(paste(goOnt, "PARENTS", sep = ""))
  goIds <- character(0)
  wantedNodes <- id
  if (children) {
    wantedNodes <- c(wantedNodes, unlist(edges(goDag(r))[id], 
                                         use.names = FALSE))
  }
  g <- reverseEdgeDirections(goDag(r))
  if (parents) {
    wantedNodes <- c(wantedNodes, unlist(edges(g)[id], use.names = FALSE))
  }
  wantedNodes <- unique(wantedNodes)
  g <- subGraph(wantedNodes, g)
  if (children) {
    for (goid in id) {
      kids <- unique(goKidsEnv[[goid]])
      for (k in kids) {
        if (is.na(k)) 
          next
        if (!(k %in% nodes(g))) {
          g <- addNode(k, g)
          g <- addEdge(k, goid, g)
        }
      }
    }
  }
  if (parents) {
    for (goid in id) {
      elders <- unique(goParentsEnv[[goid]])
      for (p in elders) {
        if (is.na(p)) 
          next
        if (!(p %in% nodes(g))) {
          g <- addNode(p, g)
          g <- addEdge(goid, p, g)
        }
      }
    }
  }
  g
}

# scater functions ----
#Plot expression values for a set of features (e.g. genes or transcripts) in a SingleExperiment object, against a continuous or categorical covariate for all cells.
# . plotExpression ----

# Plot column-level (i.e., cell) metadata in an SingleCellExperiment object
# . plotColData ---- 

#Create a base ggplot object from a SingleCellExperiment, the contents of which can be directly referenced in subsequent layers without prior specification.
# . ggcells ----

# Markers ----

# Cell type annotation ----


