
## Author: Gro Nilsen
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the clusterGenomics package
## Reference: "Identifying clusters in genomics data by recursive partitioning", Nilsen et al. (2013, preprint)


#' Compute a Distance Matrix for PART
#'
#' @description
#' `getDist()` computes a distance matrix from a data matrix `X` using either
#' standard Euclidean-based distances or correlation-based distances.
#' The function supports:
#' * squared Euclidean distance,
#' * standard distances from `stats::dist()`, and
#' * correlation-based distances, computed as `1 - cor(X)` with user-specified
#'   correlation methods.
#'
#' This helper function is used throughout the PART algorithm to construct
#' distance matrices for hierarchical clustering and threshold computations.
#'
#' @param X
#' A numeric data matrix where rows represent observations and columns represent
#' features.
#'
#' @param dist.method
#' A character string specifying which distance metric to use. Supported options:
#' * `"sq.euclidean"` — squared Euclidean distance,
#' * `"cor"` — correlation-based distance,
#' * any other string accepted by `stats::dist()` (e.g., `"euclidean"`,
#'   `"manhattan"`, `"maximum"`, etc.).
#'
#' @param cor.method
#' Method for computing correlations when `dist.method = "cor"`.
#' Passed to `stats::cor()`. Options include `"pearson"` (default),
#' `"spearman"`, or `"kendall"`.
#' @export
#'
getDist <- function(X,dist.method,cor.method="pearson"){
  #Calculate distance matrix:
  if(dist.method=="sq.euclidean"){
    # Squared Euclidean distance
    dX <- dist(X,method="euclidean")^2
  }else if(dist.method=="cor"){
    #cor computes correlation between columns
    # Correlation-based distance: 1 - cor(row_i, row_j)
    # cor() computes column correlations to transpose X
    dX <- as.dist(1-cor(t(X),method=cor.method))
  }else{
    # Any distance supported by stats::dist()
    dX <- dist(X,method=dist.method)
  }

  return(dX)

}

#' Wrapper for Hierarchical Clustering
#'
#' @description
#' `doHclust()` performs hierarchical clustering using a specified linkage method
#' and returns both the full `hclust` object and the cluster labels obtained
#' by cutting the dendrogram into `k` clusters.
#'
#' This wrapper is used throughout PART to standardize calls to hierarchical
#' clustering and ensure consistent linkage and cutting behavior.
#'
#' @param d
#' A distance object (typically produced by `dist()` or `getDist()`) representing
#' pairwise distances between observations.
#'
#' @param k
#' Integer specifying the number of clusters to extract from the dendrogram
#' using `stats::cutree()`.
#'
#' @param linkage
#' The hierarchical clustering linkage method. Must be a valid method for
#' `stats::hclust()`, such as `"average"` (default in PART), `"ward.D2"`,
#' `"complete"`, `"single"`, etc.
#' @importFrom stats hclust cutree
#' @export
#'
doHclust <- function(d,k,linkage){
  cl <- hclust(d,method=linkage)
  lab <- cutree(cl,k)
  return(list(cl=cl,lab=lab))
}

#' Wrapper for K-means Clustering
#'
#' @description
#' `doKmeans()` performs k-means clustering on a data matrix `X` and returns both
#' the full `kmeans` object along with the resulting cluster labels.
#'
#' This function standardizes how PART calls k-means and ensures consistent use
#' of the `nstart` parameter. K-means in this context always uses Euclidean
#' distance.
#'
#' @param X
#' A numeric data matrix where rows represent observations and columns represent
#' features.
#'
#' @param k
#' The number of clusters to form. Must be at least 1 and strictly less than the
#' number of observations.
#'
#' @param nstart
#' The number of random initializations used in the k-means algorithm.
#' Passed directly to `stats::kmeans()` to improve robustness of results.
#'
#' @importFrom stats kmeans
#' @export
#'
doKmeans <- function(X,k,nstart){
  cl <- kmeans(X,k,nstart=nstart) #distance is always euclidean
  lab <- cl$cluster
  return(list(cl=cl,lab=lab))
}


#' Generate Cluster Partitions for K = 1, …, Kmax
#'
#' @description
#' `findPartition()` computes cluster assignments for a sequence of cluster
#' counts `k = 1, …, Kmax` using either hierarchical clustering or k-means,
#' depending on the method specified in the `...` argument list.
#'
#' For hierarchical clustering, a distance matrix is computed (unless provided),
#' and `cutree()` is applied at each level `k`.
#' For k-means, the algorithm is run independently for each `k`, using Euclidean
#' distance.
#'
#' The function returns a list where each element `cl.lab[[k]]` contains a vector
#' of cluster labels corresponding to the partition with `k` clusters.
#'
#' @param X
#' A numeric data matrix where rows correspond to observations and columns to
#' features. Used directly for k-means and for computing distances when needed.
#'
#' @param Kmax
#' Maximum number of clusters to generate. Partitions for all `k` from `1` to
#' `Kmax` are returned.
#'
#' @param dX
#' Optional precomputed distance object. If `NULL` and the clustering method is
#' hierarchical (`"hclust"`), the distance matrix is computed internally using
#' `getDist()`. Ignored when `cl.method = "kmeans"`.
#'
#' @param ...
#' Additional parameters passed from higher-level PART functions, provided as a
#' list (typically `fixed.par`). Relevant values include `cl.method`, `linkage`,
#' `nstart`, `dist.method`, and `cor.method`.
#'
#' @return list of cluster labels
#'
#' @export
#'
findPartition <- function(X,Kmax,dX=NULL,...){

  # Collect algorithm parameters passed from PART (cl.method, linkage, etc.)
  arg <- as.list(...)

  #Calculate distances:
  if(is.null(dX) && arg$cl.method!="kmeans"){
    dX <- getDist(X,dist.method=arg$dist.method,cor.method=arg$cor.method)
  }

  #Find the partition for k=1,..,Kmax and return cluster labels stored in a list
  # Prepare output: list where cl.lab[[k]] holds labels for k clusters
  cl.lab <- vector("list",Kmax)
  for(k in 1:Kmax){
    labX <- switch(arg$cl.method,
                   hclust=doHclust(dX,k=k,linkage=arg$linkage)$lab,
                   kmeans=doKmeans(X,k,nstart=arg$nstart)$lab)    #note:kmeans only calculated for euclidean distance!

    cl.lab[[k]] <- labX
  }
  return(cl.lab)
}

#' Compute Gap Statistic for Selecting the Number of Clusters
#'
#' @description
#' `gap()` computes the Gap statistic of Tibshirani et al. (2001) for a dataset
#' `X` over candidate cluster numbers `K = 1, …, Kmax`.
#'
#' The procedure evaluates:
#' * `W_k` — within-cluster dispersion for the observed dataset,
#' * `W^*_k` — reference dispersions from `B` Monte Carlo samples generated
#'    according to the reference distribution (`ref.gen`),
#' * `gap_k = mean(log W^*_k) − log W_k`,
#' * the standard error `sk`, and
#' * the Gap selection criterion:
#'     select the smallest `k` such that
#'      `gap_k ≥ gap_{k+1} − sk_{k+1}`.
#'
#' The function returns both the selected number of clusters `hatK` and the
#' corresponding cluster labels.
#'
#' The function can be called internally by PART—where parameters are passed via
#' `fixed.par`—or used independently, in which case arguments in `...` override
#' default values directly.
#'
#'
#' @param X
#' A numeric data matrix where rows are observations and columns are features.
#'
#' @param Kmax
#' Maximum number of clusters to evaluate. If `Kmax` exceeds the number of
#' observations, it is reduced appropriately.
#'
#' @param B
#' Number of Monte Carlo reference datasets to generate for estimating the
#' reference dispersion `W^*_k`.
#'
#' @param ref.gen
#' Method used to generate reference datasets. Typically `"PC"` (principal
#' components) but may include other options supported by `getReferenceW()`.
#'
#' @param cl.lab
#' Optional list of cluster label vectors for `k = 1, …, Kmax`. If `NULL`,
#' labels are computed by `findPartition()`.
#'
#' @param ...
#' Additional tuning parameters for clustering and distance computation.
#' These may be supplied as a `fixed.par` list (when called from PART) or as
#' individual named arguments (for stand-alone usage). Relevant fields include
#' `cl.method`, `dist.method`, `linkage`, `cor.method`, and others.
#'
#' @return list of gaps found
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
  # Adjust Kmax if too large for the current dataset
  if(n<=Kmax){
    if(fixed.par$cl.method=="hclust"){
      Kmax = n
    }else{
      Kmax = n-1   #kmeans does not work if K=n
    }
  }

  # Compute the distance matrix used for hierarchical clustering
  dX <- getDist(X,dist.method=fixed.par$dist.method,cor.method=fixed.par$cor.method)

  # Compute cluster labels for all K if none were provided
  if(is.null(cl.lab)){
    cl.lab <- findPartition(X=X,Kmax=Kmax,dX=dX,fixed.par)
  }

  #Find W: vector containing W_K for all choices of K
  W <- findW(dX=dX,K=Kmax,cl.lab=cl.lab)


  # Generate B reference datasets and compute W*_k for all K
  # Wb is a B × Kmax matrix of reference dispersions
  Wb <- getReferenceW(X,Kmax,B,ref.gen,fixed.par)


  #Calculate gap statistic:
  #gap_k = mean(log W*_k) − log(W_k)
  L <- apply(log(Wb),1,mean)
  gap <- L - log(W)

  # Compute standard error for each k
  sdk <- apply(log(Wb),1,sd)*sqrt((B-1)/B)  #multiply by the last expression to get 1/B (as in original article) instead of 1/(B-1) (normal calculation)
  sk <- sqrt(1+(1/B))*sdk

  #Calculate gap-criterion for each k:
  diff <- gap[1:Kmax-1] - (gap[2:Kmax] - sk[2:Kmax])

  # Select smallest k satisfying the criterion
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

#' Generate Uniform Reference Data for a Single Variable
#'
#' @description
#' `sim()` generates a univariate reference sample by drawing values uniformly
#' between the minimum and maximum of the input vector `Xcol`.
#' This is used when creating reference datasets for computing the Gap statistic
#' and other resampling-based clustering criteria.
#'
#' @param Xcol
#' A numeric vector. The simulated values are drawn uniformly from the interval
#' spanning `min(Xcol)` to `max(Xcol)`.
#' @return Simulated uniform distribution
#' @export
#'
sim <- function(Xcol) {
  min <- min(Xcol)
  max <- max(Xcol)
  U <- runif(length(Xcol),min,max)
  return(U)
}

#' Compute Within-Cluster Dispersion for Multiple Clusterings
#'
#' @description
#' `findW()` calculates the within-cluster dispersion \(W_k\) for a distance
#' matrix `dX` and cluster labels for different numbers of clusters.
#' *****
#' Within-cluster dispersion tells us how ‘tight’ each cluster is.
#' Tighter clusters are better. We use it to decide if we should split
#' clusters further, and to figure out the optimal number of clusters.
#' *****
#'
#' For each `k = 1, …, K`, the function sums the pairwise distances between
#' points within each cluster and scales by the cluster size:
#'
#' When `k = 1`, `W_1` is simply the sum of all pairwise distances divided by
#' twice the number of points.
#'
#' This function is typically used within Gap statistic computations.
#'
#' @param dX
#' A distance object (`dist`) or a symmetric matrix of pairwise distances
#' between observations.
#'
#' @param K
#' Maximum number of clusters to compute dispersions for.
#'
#' @param cl.lab
#' A list of integer vectors, where `cl.lab[[k]]` contains cluster labels
#' for `k` clusters.
#'
#' @return The within cluster distribution
#'
#' @export
#'
findW <- function(dX,K,cl.lab){

  # Initialize vector to store within-cluster dispersion W_k
  W <- rep(0,K)

  # Convert dX to a full matrix if it is a 'dist' object
  n <- nrow(as.matrix(dX))

  #Calculate W for k=1:
  #Compute W for k = 1 (entire dataset as a single cluster)
  W[1] <- sum(as.matrix(dX))/(2*n)

  # Step 2: Compute W for k = 2, …, K
  k <- 2
  while(k<=K){
    labX <- cl.lab[[k]]
    # Compute within-cluster dispersion for each cluster
    for(i in 1:k){
      # Submatrix of pairwise distances within cluster i
      d.k <- as.matrix(dX)[labX==i,labX==i]
      # Sum of distances within cluster i
      D.k <- sum(d.k)
      # Number of observations in cluster i
      nk <- nrow(as.matrix(d.k))
      # Add scaled dispersion to W[k]
      W[k] <- W[k] + D.k/(2*nk)

    }#endfor
    k <- k+1
  }#endwhile
  return(W)

}#endfunction


#' Generate Reference Within-Cluster Dispersion Matrices
#'
#' @description
#' `getReferenceW()` generates multiple reference datasets and calculates the
#' within-cluster dispersion \(W_b\) for each, across 1 to Kmax clusters.
#'
#' Reference datasets can be generated either:
#' - Uniformly within the original data range for each variable (`ref.gen="uniform"`)
#' - Using a principal components transformation to preserve variance structure (`ref.gen="PC"`)
#'
#' This function is used in the computation of the Gap statistic to compare
#' the observed clustering dispersion to what is expected under a null reference.
#'
#' @param X
#' Numeric data matrix (observations × variables) to generate reference datasets from.
#'
#' @param Kmax
#' Maximum number of clusters to compute W for.
#'
#' @param B
#' Number of reference datasets to generate.
#'
#' @param ref.gen
#' Reference generation method:
#' - `"PC"`: uses PCA-based transformation for preserving variance structure
#' - any other value: generates uniform reference data per variable.
#'
#' @param ...
#' Additional arguments passed on to other functions (e.g., `dist.method` for distance calculation, `cl.method` for clustering method, `linkage`, `cor.method`, `nstart`).
#'
#' @return calculated Wb
#'
#' @export
#'
getReferenceW <- function(X,Kmax,B,ref.gen,...)	{
  arg <- as.list(...)

  # Initialize matrix to store within-cluster dispersion for reference datasets
  # Each column = a reference dataset, each row = W for k = 1..Kmax
  Wb <- matrix(0,nrow=Kmax,ncol=B)

  # If using PCA-based reference, transform X to principal component space
  if(ref.gen=="PC"){
    #center columns
    m <- apply(X,2,mean,na.rm=TRUE)   #First columncenter X:
    Xc <- sweep(X,2,m)
    #SVD:
    s <- svd(Xc)
    # Rotate data into principal component space
    newX <- Xc%*%s$v
  }

  # Loop over B reference datasets
  for(b in 1:B){
    if(ref.gen=="PC"){
      # Generate uniform reference values in PC space
      U <- apply(newX,2,sim)
      Z1 <- U%*%t(s$v)      #backtransform
      # Re-add original column means
      Z <- sweep(Z1,2,m,FUN="+")
    }else{
      # Generate uniform reference data per variable
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

#' Perform Recursive Partitioning Using the PART Algorithm
#'
#' @description
#' The `part()` function performs clustering using the PART algorithm
#' (Partitioning Algorithm based on Recursive Thresholding). The method begins
#' with a global clustering of the data and then recursively evaluates potential
#' subclusters. At each recursive step, it applies the Gap statistic and a
#' separation threshold (`minDist`) to determine whether a cluster should be
#' further split. The procedure is designed to detect both major cluster
#' structure and finer substructure while guarding against spurious splits.
#'
#' The function returns the estimated number of clusters, cluster labels for all
#' observations, and a list of samples flagged as potential outliers.
#'
#'
#' @param X
#' A numeric data matrix with observations in rows and variables in columns.
#' This is the dataset to be clustered.
#'
#' @param Kmax
#' The maximum number of clusters to consider during the *initial*, global PART
#' run. Recursive subdivision uses its own upper limit (`Kmax.rec`) provided via
#' `...`. Default is 10.
#'
#' @param minSize
#' The minimum number of observations allowed in a cluster in order for that
#' cluster to be considered for further recursive splitting. Clusters smaller
#' than `minSize` are treated as terminal. Default is 8.
#'
#' @param minDist
#' Minimum required separation (in dendrogram height units) between subclusters
#' for a recursive split to be accepted. If `NULL`, a separation threshold is
#' computed automatically using `get.threshold()` based on the user-specified
#' value of `q` (passed via `...`). Default is `NULL`.
#'
#' @param cl.lab
#' Optional list of precomputed cluster label vectors. If supplied, the list
#' should contain the cluster assignments for `k = 1, …, Kmax`, typically
#' produced by hierarchical clustering. If `NULL`, cluster labels are generated
#' internally. Default is `NULL`.
#'
#' @param ...
#' Additional arguments controlling the PART algorithm. These override the
#' internal defaults defined in the function. Common parameters include:
#' \itemize{
#'   \item `q` — proportion of dendrogram height used to compute the stopping
#'   threshold when `minDist` is not supplied. Default: `0.25`.
#'   \item `Kmax.rec` — maximum number of clusters evaluated at each recursive
#'   splitting step. Default: `5`.
#'   \item `B` — number of bootstrap samples used in the Gap statistic. Default:
#'   `100`.
#'   \item `ref.gen` — method used to generate reference datasets for the Gap
#'   statistic (“PC” by default).
#'   \item `dist.method` — distance metric used for hierarchical clustering.
#'   \item `cl.method` — clustering method (“hclust” or “kmeans”).
#'   \item `linkage` — linkage method if hierarchical clustering is used.
#'   \item `cor.method` — correlation type if correlation distance is selected.
#'   \item `nstart` — number of random initializations if k-means is used.
#' }
#' @export
part <- function(X,Kmax=10,minSize=8,minDist=NULL,cl.lab=NULL,...){

  #Default ... values:
  default.par <- list(q=0.25,Kmax.rec=5,B=100,ref.gen="PC",dist.method="euclidean",cl.method="hclust",linkage="average",cor.method="pearson",nstart=10)
  #Check for user modifications:
  fixed.par <- c(minDist=minDist,minSize=minSize,modifyList(default.par,list(...)))

  #If the user did NOT provide minDist, compute it adaptively.
  if(is.null(minDist)){
    minDist <- get.threshold(X,q=fixed.par$q,fixed.par)
    fixed.par$minDist <- minDist
  }
  print(minDist)

  # PartRec applies:
  #   - global clustering (up to Kmax)
  #   - recursive splitting (up to Kmax.rec)
  #   - Gap statistic at each step
  #   - stopping rules using minDist and minSize
  #
  # 'ind' starts as all ones, meaning the entire dataset is one group initially.
  clusters = PartRec(X,Kmax=Kmax,ind=rep(1,nrow(X)),cl.lab=cl.lab,fixed.par)

  # getPARTlabels merges recursive output into a single label vector,
  # marking small unstable clusters as outliers (label = 0).
  label <- getPARTlabels(clusters,minSize)
  outliers <- which(label==0)
  if(length(outliers)==0){
    outliers <- NULL
    hatK <- length(unique(label))
  }else{
    hatK <- length(unique(label[-outliers]))
  }

  # hatK      : number of non-outlier clusters
  # lab.hatK  : cluster labels for each observation (0 = outlier)
  # outliers  : index vector identifying outlier samples
  return(list(hatK=hatK,lab.hatK=label,outliers=outliers))
}


#' Recursive Splitting Step of the PART Algorithm
#'
#' @description
#' `PartRec()` performs the recursive splitting step of the PART clustering
#' algorithm. Given a subset of observations and a maximum number of clusters to
#' evaluate, the function:
#'
#' 1. Checks whether the current subset is large enough to split at all.
#' 2. Computes candidate partitions for `K = 1, …, Kmax` using a user-specified
#'    clustering method.
#' 3. Selects the optimal number of clusters (`hatK`) using the Gap statistic.
#' 4. Applies size constraints (`minSize`) to ensure meaningful splits.
#' 5. When a split is accepted, recursively applies PART to each subcluster.
#'
#' The output is either:
#' * a vector representing an unsplit terminal cluster, or
#' * a matrix whose columns represent recursively generated subclusters.
#'
#'
#' @param X
#' A numeric data matrix containing the observations in the cluster currently
#' under consideration.
#'
#' @param Kmax
#' Maximum number of clusters to consider for the *current* (non-recursive)
#' partitioning step.
#'
#' @param ind
#' A binary index vector of the same length as the number of rows in the
#' original dataset. Values of `1` indicate membership in the current subset and
#' `0` indicate exclusion. This index is updated recursively to track cluster
#' identities across all levels of the dendrogram.
#'
#' @param cl.lab
#' Optional list containing pre-computed cluster label vectors for
#' `K = 1, …, Kmax`. If supplied, these labels are used instead of recomputing
#' partitions. If `NULL`, labels are obtained through `findPartition()`. Default
#' is `NULL`.
#'
#' @param ...
#' Additional parameters passed from the top-level `part()` function, provided as
#' a list in `fixed.par`. Important values include `minSize`, `minDist`,
#' `Kmax.rec`, `dist.method`, `cl.method`, `linkage`, `cor.method`, and others
#' controlling the clustering algorithm and recursive stopping rules.
PartRec <- function(X,Kmax,ind,cl.lab=NULL,...){

  fixed.par <- as.list(...)

  #STEP 1: Make sure it is feasible to split X into two clusters each of size >= minSize,
  #otherwise return this cluster
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


  # --
  # Case A — hatK == 1: No valid split accepted
  #
  # Attempt:
  #   - Form two tentative clusters from cl.lab[[2]] (K=2 solution)
  #   - Test whether dendrogram height between them exceeds minDist
  #   - If yes, recursively evaluate both tentative clusters
  #   - If no, treat as terminal cluster
  #
  # Special case: if only one cluster level exists (length(cl.lab) == 1),
  # splitting is impossible → return terminal cluster.
  # --
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
    # --
    # Case B — hatK > 1:
    #
    # Accept the K = hatK split, and recursively split each resulting cluster.
    # --

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


#' Compute Stopping Threshold for Recursive Clustering
#'
#' @description
#' `get.threshold()` calculates a threshold distance used to decide whether
#' a cluster should be split further in recursive clustering (e.g., PART).
#' The threshold is based on the hierarchical clustering dendrogram of the
#' data and a quantile parameter `q`.
#'
#' @param X
#' Numeric data matrix (observations × variables) to calculate distances from.
#'
#' @param q
#' Quantile parameter (0 < q < 1). The threshold is set as the (1-q) quantile
#' of all dendrogram heights. Higher q → lower threshold → more splitting.
#'
#' @param ...
#' Additional parameters passed on to distance and clustering functions (e.g.,
#' `dist.method`, `linkage`, `cor.method`).
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


#' Assign Final Cluster Labels in PART
#'
#' @description
#' `getPARTlabels()` takes the output of the recursive clustering (PART)
#' and assigns final cluster labels to each observation. Clusters smaller
#' than `minSize` are labeled as 0 (treated as outliers).
#'
#' @param clusters
#' A matrix (or vector) indicating the tentative clusters from PART. Each column
#' corresponds to a recursive split, with 1 indicating membership in that subcluster.
#'
#' @param minSize
#' Minimum number of observations required for a cluster to be considered valid.
#' Clusters smaller than this are treated as outliers.
#'
#' @return
#' An integer vector of cluster labels for each observation. Outliers are labeled 0.
#'
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




