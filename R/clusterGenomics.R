
####################################################################################################################
## Author: Gro Nilsen
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the clusterGenomics package
## Reference: "Identifying clusters in genomics data by recursive partitioning", Nilsen et al. (2013, preprint)
####################################################################################################################

#' @title getDist
#'
#' @description Compute distance matrix for given distance measure; note: distance between rows!
#'
#' @param X Matrix used to calculate the distance
#' @param dist.method The distance method to use in the calculation
#' @param cor.method The correlation method to use. 'pearson' is the default,
#'
#' @importFrom stats as.dist dist
#'
#' @export
#'
getDist <- function(X,dist.method,cor.method="pearson"){
  #Calculate distance matrix:
  if(dist.method=="sq.euclidean"){
    dX <- dist(X,method="euclidean")^2
  }else if(dist.method=="cor"){
    #cor computes correlation between columns
    dX <- as.dist(1-cor(t(X),method=cor.method))
  }else{
    dX <- dist(X,method=dist.method)
  }

  return(dX)

}

#' @title doHclust
#'
#' @description Perform hierarchical clustering with given linkage, and apply horizontal cutting of tree into k clusters
#'
#' @param d a dissimilarity structure as produced by dist.
#' @param k an integer scalar or vector with the desired number of groups
#' @param linkage the agglomeration method to be used.
#'
#' @importFrom stats hclust cutree
#'
#' @export
#'
doHclust <- function(d,k,linkage){
  cl <- hclust(d,method=linkage)
  lab <- cutree(cl,k)
  return(list(cl=cl,lab=lab))
}

#' @title doKmeans
#'
#' @description Perform K-means clustering into k clusters
#'
#' @param X matrix of values
#' @param k value for kmeans
#' @param nstart Where to start
#'
#' @importFrom stats kmeans
#' @export
#'
doKmeans <- function(X,k,nstart){
  cl <- kmeans(X,k,nstart=nstart) #distance is always euclidean
  lab <- cl$cluster
  return(list(cl=cl,lab=lab))
}


#' @title findPartition
#'
#' @description Obtain partitions of the data into k=1,...,Kmax clusters
#' gives the clustering method, distance measure, and other parameters to be used in the clustering
#'
#' @param matrix of values
#' @param maximum number of clusters
#' @param dx distance measure
#'
#' @export
#'
findPartition <- function(X,Kmax,dX=NULL,...){

  arg <- as.list(...)

  #Calculate distances:
  if(is.null(dX) && arg$cl.method!="kmeans"){
    dX <- getDist(X,dist.method=arg$dist.method,cor.method=arg$cor.method)
  }

  #Find the partition for k=1,..,Kmax and return cluster labels stored in a list
  cl.lab <- vector("list",Kmax)
  for(k in 1:Kmax){
    labX <- switch(arg$cl.method,
                   hclust=doHclust(dX,k=k,linkage=arg$linkage)$lab,
                   kmeans=doKmeans(X,k,nstart=arg$nstart)$lab)    #note:kmeans only calculated for euclidean distance!

    cl.lab[[k]] <- labX
  }
  return(cl.lab)
}

#' @title gap
#'
#' @description gap finds the optimal number of clusters based on the "gap statistic" (Tibshirani et al. 2001)
#'
#' @param X matrix of values
#' @param Kmax maximum number of clusters
#' @param B number of recursions
#' @param ref.gen reference gen
#' @param cl.lab cluster labels
#'
#' @export
#'
gap <- function(X,Kmax=10,B=100,ref.gen="PC",cl.lab=NULL,...){

  #Default ... values:
  default.par <- list(dist.method="euclidean",cl.method="hclust",linkage="average",cor.method="pearson",nstart=10)
  #Check for user modifications
  #Note: call could come from part in which case ... is a list called fixed.par, or it could be a independent call in which case ... could contain several parameters which must be converted to list
  if(hasArg(fixed.par)){
    fixed.par <- modifyList(default.par,as.list(list(...)$fixed.par))
  }else{
    fixed.par <- modifyList(default.par,list(...))
  }

  n <- nrow(X)

  if(n<=Kmax){
    if(fixed.par$cl.method=="hclust"){
      Kmax = n
    }else{
      Kmax = n-1   #kmeans does not work if K=n
    }
  }

  #Calculate distances:
  dX <- getDist(X,dist.method=fixed.par$dist.method,cor.method=fixed.par$cor.method)

  if(is.null(cl.lab)){
    cl.lab <- findPartition(X=X,Kmax=Kmax,dX=dX,fixed.par)
  }

  #Find W: vector containing W_K for all choices of K
  W <- findW(dX=dX,K=Kmax,cl.lab=cl.lab)


  #Generate reference data sets and find Wb:
  Wb <- getReferenceW(X,Kmax,B,ref.gen,fixed.par)


  #Calculate gap statistic:
  L <- apply(log(Wb),1,mean)
  gap <- L - log(W)

  sdk <- apply(log(Wb),1,sd)*sqrt((B-1)/B)  #multiply by the last expression to get 1/B (as in original article) instead of 1/(B-1) (normal calculation)
  sk <- sqrt(1+(1/B))*sdk

  #Calculate gap-criterion for each k:
  diff <- gap[1:Kmax-1] - (gap[2:Kmax] - sk[2:Kmax])

  #Select the first k where the criterion is positive:
  kvec <- 1:Kmax
  posDiff <- diff[diff>=0]
  if(length(posDiff)==0){
    hatK <- 1  #no k satisfies criterion, define 1 cluster e.g. to make part run..
  }else{
    hatK <- kvec[diff>=0][1]
  }

  if(all(is.na(diff))){
    hatK <- 1   #when used with part; can happen if sub-cluster consists of only 1 or 2 objects
  }

  #Get labels for best partition:
  lab.hatK <- cl.lab[[hatK]]

  return(list(hatK=hatK,lab.hatK=lab.hatK,gap=gap,sk=sk,W=W))
}

#' @title sim
#'
#' @description helper-functions only used by GAP:
#' Simulate from a uniform distribution according to min and max of the feature vector:
#'
#' @export
#'
sim <- function(Xcol) {
  min <- min(Xcol)
  max <- max(Xcol)
  U <- runif(length(Xcol),min,max)
  return(U)
}

#' @title findW
#'
#' @description Calulates the total within-cluster dispersion (W_K) for a k=1,...,K:
#'
#' @param dx matrix for cluster
#' @param K current cluster
#' @param cl.lab cluster labels
#'
#' @export
#'
findW <- function(dX,K,cl.lab){

  W <- rep(0,K)
  #n <- nrow(X)
  n <- nrow(as.matrix(dX))

  #Calculate W for k=1:
  W[1] <- sum(as.matrix(dX))/(2*n)

  k <- 2
  while(k<=K){
    labX <- cl.lab[[k]]

    for(i in 1:k){
      #Sum of within-cluster dispersion:
      d.k <- as.matrix(dX)[labX==i,labX==i]
      D.k <- sum(d.k)

      nk <- nrow(as.matrix(d.k))
      W[k] <- W[k] + D.k/(2*nk)

    }#endfor
    k <- k+1
  }#endwhile
  return(W)

}#endfunction


#' @title getReferenceW
#'
#' @description Generate reference data sets and find Wb:
#'
#' @param X matrix of values
#' @param Kmax maximum number of clusters
#' @param B Number of recursions
#'
#' @export
#'
getReferenceW <- function(X,Kmax,B,ref.gen,...)	{
  arg <- as.list(...)
  #Generate reference data sets and find Wb:
  Wb <- matrix(0,nrow=Kmax,ncol=B)
  if(ref.gen=="PC"){
    #Transform data using svd:
    m <- apply(X,2,mean,na.rm=TRUE)   #First columncenter X:
    Xc <- sweep(X,2,m)
    #SVD:
    s <- svd(Xc)
    newX <- Xc%*%s$v
  }
  for(b in 1:B){
    if(ref.gen=="PC"){
      U <- apply(newX,2,sim)
      Z1 <- U%*%t(s$v)      #backtransform
      #Add mean
      Z <- sweep(Z1,2,m,FUN="+")
    }else{
      Z <- apply(X,2,sim)
    }
    #Calculate distances
    dZ <- getDist(Z,dist.method=arg$dist.method,cor.method=arg$cor.method)
    #Cluster reference data
    clW.lab <- findPartition(X=Z,Kmax=Kmax,dX=dZ,arg)
    #Calculate Wb_K for all values of K
    Wb[,b] <- findW(dX=dZ,K=Kmax,cl.lab=clW.lab)
  }#endfor

  return(Wb)

}

#' @title part
#'
#' @description Main function for PART:
#'
#' @param X matrix of values
#' @param kmax maximum number of clusters
#' @param minSize minimum cluster size
#' @param minDist Minimum distance between clusters
#' @param cl.lab cluster labels
#'
#' @export
#'
part <- function(X,Kmax=10,minSize=8,minDist=NULL,cl.lab=NULL,...){

  #Default ... values:
  default.par <- list(q=0.25,Kmax.rec=5,B=100,ref.gen="PC",dist.method="euclidean",cl.method="hclust",linkage="average",cor.method="pearson",nstart=10)
  #Check for user modifications:
  fixed.par <- c(minDist=minDist,minSize=minSize,modifyList(default.par,list(...)))


  #Find stopping threshold if minDist is NULL
  if(is.null(minDist)){
    minDist <- get.threshold(X,q=fixed.par$q,fixed.par)
    fixed.par$minDist <- minDist
  }

  #Start recursive runs:
  clusters = PartRec(X,Kmax=Kmax,ind=rep(1,nrow(X)),cl.lab=cl.lab,fixed.par)


  #Check for possible outliers and assign cluster labels
  label <- getPARTlabels(clusters,minSize)
  outliers <- which(label==0)
  if(length(outliers)==0){
    outliers <- NULL
    hatK <- length(unique(label))
  }else{
    hatK <- length(unique(label[-outliers]))
  }


  return(list(hatK=hatK,lab.hatK=label,outliers=outliers))
}


#' @title PartRec
#'
#' @description The recursive function used by PART clustering
#'
#' @param X Matrix of values
#' @param Kmax max number of cluster
#' @param ind number of clusters
#' @param cl.lab cluster label
#'
#' @export
#'
PartRec <- function(X,Kmax,ind,cl.lab=NULL,...){

  fixed.par <- as.list(...)

  #STEP 1: Make sure it is feasible to split X into two clusters each of size >= minSize, otherwise return this cluster
  if(sum(ind)<(2*fixed.par$minSize)){
    return(ind)
  }


  #STEP 2: Use a clustering algorithm to partition the objects in X into K=1,..,Kmax clusters:
  #First Make sure Kmax does not exceed the number of objects in X:
  n <- sum(ind)
  if(n<=Kmax){
    if(fixed.par$cl.method=="hclust"){
      Kmax = n
    }else{
      Kmax = n-1   #kmeans does not work if K=n
    }
  }
  if(is.null(cl.lab)){
    cl.lab <- findPartition(X=X,Kmax=Kmax,dX=NULL,fixed.par)
  }



  #STEP 3: #Use an objective function (Gap) to decide on the optimal number of clusters, hatK, for this set X
  gap.res <- gap(X=X,Kmax=Kmax,cl.lab=cl.lab,B=fixed.par$B,ref.gen=fixed.par$ref.gen,fixed.par=fixed.par)
  hatK <- gap.res$hatK
  lab.hatK <- gap.res$lab.hatK


  #STEP 4.

  #In case hatK > 1, make sure at least two of them are >= minSize:
  if(sum(table(lab.hatK) >= fixed.par$minSize)<2){
    hatK <- 1
  }


  #a) hatK == 1:
  if(hatK==1 && length(cl.lab)==1){
    #Special case if minSize=1 and only 2 obs in X; kmeans cannot return 2 clusters (see findPartition) and Kmax is therefore set to 1 -> cannot divide into two tentative clusters
    return(ind)
  }

  if(hatK==1){
    #Divide set into two tentative clusters:
    obs1 <- cl.lab[[2]]==1
    obs2 <- cl.lab[[2]]==2
    T1 <- X[obs1,,drop=FALSE]   #drop is necessary in case obs1 is of length 1!
    T2 <- X[obs2,,drop=FALSE]


    # Check if stopping threshold has been reached:
    #Use hierarchical clustering to determine the distance between T1 and T2:
    hc.res <- doHclust(getDist(X,dist.method=fixed.par$dist.method,cor.method=fixed.par$cor.method),k=1,linkage=fixed.par$linkage)$cl    #(k is irrelevant here, only specified because doHclust needs it)
    T.height <- max(hc.res$height)


    if(T.height > fixed.par$minDist){
      #Create new index-vectors corresponding to the two tentative clusters:
      t1.ind <- ind
      t2.ind <- ind
      t1.ind[which(ind==1)[!obs1]] <- 0
      t2.ind[which(ind==1)[!obs2]] <- 0
      #Tentative runs:
      t1 = PartRec(X=T1,Kmax=fixed.par$Kmax.rec,ind=t1.ind,cl.lab=NULL,fixed.par)
      t2 = PartRec(X=T2,Kmax=fixed.par$Kmax.rec,ind=t2.ind,cl.lab=NULL,fixed.par)
      #If no clusters are found in recursive runs the return value will be a vector, otherwise it will be a matrix:
      if(!is.matrix(t1) && !is.matrix(t2)){
        return(ind)
      }else{
        return(cbind(t1,t2))
      }
    }else{
      #Return the current set if stopping criterion has been reached
      return(ind)
    }

  }else{
    #b) hatK > 1

    res <- matrix(NA,nrow=length(ind),ncol=0)
    for(k in 1:hatK){
      #Pick out subset in this cluster
      obs = cl.lab[[hatK]]==k
      S = X[obs,,drop=FALSE]
      ind.S = ind
      ind.S[which(ind==1)[!obs]] <- 0
      res = cbind(res,PartRec(X=S,Kmax=fixed.par$Kmax.rec,ind=ind.S,cl.lab=NULL,fixed.par))
    }
    return(res)

  }

}


#' @title get.threshold
#'
#' @description Help-functions only used by part:
#' Find a stopping threshold given the percentage of heights to be used in the dendrogram
#'
#' @importFrom stats hclust
#'
#' @export
#'
get.threshold <- function(X,q,...){
  arg <- as.list(...)

  #Get distance matrix
  dX <- getDist(X,dist.method=arg$dist.method,cor.method=arg$cor.method)

  #Get hierarchical clustering result:
  cl <- hclust(dX,method=arg$linkage)

  #The total set of cluster heights in dendrogram:
  h <- cl$height

  use.h <- quantile(h,probs=1-q)

  return(use.h)

}


#' @title getPARTlabels
#'
#' @description Retrives the labels for the clusters found by PART
#'
#' @param clusters The found clusters
#' @param minSize The parameter for minimum cluster size
#'
#' @export
#'
getPARTlabels <- function(clusters,minSize){
  if(!is.matrix(clusters)){
    clusters <- as.matrix(clusters)
  }
  label <- rep(NA,nrow(clusters))
  ncl <- ncol(clusters)
  id=1
  for(j in 1:ncl){
    if(sum(clusters[,j])<minSize){
      label[clusters[,j]==1] <- 0
    }else{
      label[clusters[,j]==1] <- id
      id=id+1
    }

  }
  return(label)
}




