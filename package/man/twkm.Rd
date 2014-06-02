\name{twkm}
\alias{twkm}
\title{
  Two-level variable weighting clustering
}
\description{
  Two-level variable weighting clustering.
}
\usage{
  twkm(x, k, strGroup, lambda, eta, maxiter=100, delta=0.000001, maxrestart=10,seed=-1)
}
\arguments{
  \item{x}{numeric matrix of observations and features.}
  \item{k}{target number of clusters.}
  \item{strGroup}{a string give the group information, formatted as "0,1,2,4;3,5;6,7,8" or "0-2,4;3,4,8;6-8", where ";" defines a group}
  \item{lambda}{parameter of feature weight distribution.}
  \item{eta}{parameter of group weight distribution.}
  \item{delta}{maximum change allowed between iterations for convergence.}
  \item{maxiter}{maximum number of iterations.}

  \item{maxrestart}{maximum number of restarts. Default is 10 so that we
  stand a good chance of getting a full set of clusters. Normally, any
  empty clusters that result are removed from the result, and so we may
  obtain fewer than k clusters if we don't allow restarts(i.e.,
  maxrestart=0). If < 0, then there is no limit on the number of
  restarts and we are much likely to get a full set of k clusters.}
  \item{seed}{random seed. If it was set below 0, then a randomly generated number will be assigned.}
}
\details{
  
  The two-level variable weighting clustering algorithm is a
  extension to \link{ewkm}, which itself is a soft subspace clustering
  method.
  
  The algorithm weights subspaces in both feature groups and individual
  features.
  
  Always check the number of iterations, the number of restarts, and the
  total number of iterations as they give a good indication of whether
  the algorithm converged.

  As with any distance based algorithm, be sure to rescale your numeric
  data so that large values do not bias the clustering. A quick
  rescaling method to use is \code{\link{scale}}.

}
\value{

  Return an object of class "kmeans" and "twkm", compatible with other
  function that work with kmeans objects, such as the 'print'
  method. The object is a list with the following components in addition
  to the components of the kmeans object:
  
  \item{cluster}{A vector of integer (from 1:k) indicating the cluster
  to which each point is allocated.}
  
  \item{centers}{A matrix of cluster centers.}
  
  \item{featureWeight}{A vector of weights recording the relative
  importance of each feature.}
  
  \item{groupWeight}{A vector of group weights recording the relative
  importance of each feature group.}
  
  \item{iterations}{This report on the number of iterations before
  termination. Check this to see whether the maxiters was reached. If so
  then teh algorithm may not be converging, and thus the resulting
  clustering may not be particularly good.}
  
  \item{restarts}{The number of times the clustering restarted because
  of a disappearing cluster resulting from one or more k-means having no
  observations associated with it. An number here greater than zero
  indicates that the algorithm is not converging on a clustering for the
  given k. It is recommeded that k be reduced.}
  
  \item{totalIterations}{The total number of iterations over all
  restarts.}
  
  \item{totolCost}{The total cost calculated in the cost function.}

}
\references{

 Xiaojun Chen, Xiaofei Xu, Yunming Ye and Joshua Zhexue Huang: 
 TW-k-means: A Two-level Variable Weighting Clustering Algorithm for Multi-view Data. 
 IEEE Transactions on Knowledge and Data Engineering, 2012.
}
\author{

  Xiaojun Chen \email{xjchen.hitsz@gmail.com}

}
\seealso{
  \code{\link{kmeans}}
  \code{\link{ewkm}}
  \code{\link{fgkm}}
}
\examples{

# The data twkm.sample has 2000 objects and 410 variables.
# Scale the data before clustering
x <- scale(twkm.sample[,1:409])

# Group information is formated as below.
# Each group is separated by ';'.
strGroup <- "0-75;76-291;292-355;356-402;403-408"

\donttest{
# Use the twkm algorithm.
mytwkm <- twkm(x, 10, strGroup, 3, 1)

# You can print the clustering result now.
mytwkm$cluster
mytwkm$featureWeight
mytwkm$groupWeight
mytwkm$iterations
mytwkm$restarts
mytwkm$totiters
mytwkm$totss

# Use a cluster validation method from package 'clv'.

# real.cluster is the real class label of the data 'twkm.sample'.
real.cluster <- twkm.sample[,410]

# std.ext() returns four values SS, SD, DS, DD.
std <- std.ext(as.integer(mytwkm$cluster), real.cluster)

# Rand index 
clv.Rand(std)

# Jaccard coefficient
clv.Jaccard(std)
}

}
\keyword{
  variable weighting
}

