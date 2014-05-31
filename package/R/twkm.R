twkm <- function(x, k, strGroup, lambda, eta, maxiter=100, delta=0.000001, maxrestart=10,seed=-1) 
{
  if (missing(k))
    stop("the number of clusters 'k' must be provided")
   
  if(seed<=0){
    seed <-runif(1,0,10000000)[1]
  }
    
  vars <- colnames(x)
  
  nr <-nrow(x) # nrow() return a integer type
  nc <-ncol(x) # integer
  
 # get the setting of feature group
  G <- .C("parseGroup",as.character(strGroup),numGroups=integer(1), groupInfo=integer(nc),PACKAGE="wskm")
  
  
  Z <- .C("twkm",
          x = as.double(as.matrix(x)),
          nr,
          nc,
          k = as.integer(k),
          lambda = as.double(lambda),
          eta = as.double(eta),
          G$numGroups,
          G$groupInfo,
          delta = as.double(delta),
          maxIterations = as.integer(maxiter),
          maxRestarts = as.integer(maxrestart),
          seed,
          cluster = integer(nr),
          centers = double(k * nc),
          featureWeight = double( nc),
          groupWeight = double( G$numGroups),
          iterations = integer(1),
          restarts = integer(1),
          totiters = integer(1),
          totalCost = double(1),
          totss = double(1),
		  withiness = double(k),
          PACKAGE="wskm"
          )
       
  centers <- matrix( Z$centers)
  dim(centers) <- c(k, nc)
  colnames(centers) <- vars
  
  featureWeight <- Z$featureWeight
  
  groupWeight <- Z$groupWeight
  
  ignore <- which(rowSums(centers==0) == ncol(centers))
  if (length(ignore)) {
    centers <- centers[-ignore,, drop=FALSE]
    featureWeight <- featureWeight[-ignore,, drop=FALSE]
  }
  
  rownames(centers) <- 1:nrow(centers)
  
  cluster <- Z$cluster + 1
  
  size <- aggregate(cluster, list(cluster=cluster), length)[[2]]
  
  result <- list(cluster = cluster,
                 centers = Z$centers,
                 totss = Z$totss, 
                 withinss = Z$withinss, 
                 tot.withinss = sum(Z$withiness), 
                 betweenss = Z$totss-sum(Z$withinss),
                 size = size,
                 iterations = Z$iterations,
                 restarts = Z$restarts,
		 totiters=Z$totiters,
                 featureWeight = Z$featureWeight,
                 groupWeight = Z$groupWeight)
  
  dim(result$centers) <- c(k, nc)
  
  class(result) <- c("kmeans", "twkm")
  return(result)
}
