\name{fgkm}
\alias{fgkm}
\title{
  Feature Group Weighting K-Means for Subspace clustering
}
\description{
  Perform an feature group weighting subspace k-means.
}
\usage{
  fgkm(x, centers, groups, lambda, eta, maxiter=100, delta=0.000001,
       maxrestart=10,seed=-1)
}
\arguments{
  \item{x}{numeric matrix of observations and features.}
  \item{centers}{target number of clusters or the initial centers for
    clustering.}
  \item{groups}{a string give the group information, formatted as
    "0,1,2,4;3,5;6,7,8" or "0-2,4;3,5;6-8", where ";" defines a group;
    or a vector of length of features, each element of the vector
    indicates the group of the feature. For example, }
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
  \item{seed}{random seed. If it was set below 0, then a randomly
  generated number will be assigned.}
}
\details{
  
  The feature group weighting k-means clustering algorithm is a
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

  Return an object of class "kmeans" and "fgkm", compatible with other
  function that work with kmeans objects, such as the 'print'
  method. The object is a list with the following components in addition
  to the components of the kmeans object:
  
  \item{cluster}{A vector of integer (from 1:k) indicating the cluster
  to which each point is allocated.}
  
  \item{centers}{A matrix of cluster centers.}
  
  \item{featureWeight}{A matrix of weights recording the relative
  importance of each feature for each cluster.}
  
  \item{groupWeight}{A matrix of group weights recording the relative
  importance of each feature group for each cluster.}
  
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

  X. Chen, et al., A feature group weighting method for subspace
  clustering of high-dimensional data, Pattern Recognition(2011),
  doi:10.1016/j.patcog.2011.06.004

}
\author{

  Longfei Xiao \email{lf.xiao@siat.ac.cn}

}
\seealso{
  \code{\link{kmeans}}
  \code{\link{ewkm}}
  \code{\link{twkm}}
}
\examples{

# The data fgkm.sample has 600 objects and 50 dimensions.
# Scale the data before clustering
x <- scale(fgkm.sample)

# Group information is formated as below.
# Each group is separated by ';'.
strGroup <- "0-9;10-19;20-49"

# Use the fgkm algorithm.
myfgkm <- fgkm(x, 3, strGroup, 3, 1)

# You can print the clustering result now.
myfgkm$cluster
myfgkm$featureWeight
myfgkm$groupWeight
myfgkm$iterations
myfgkm$restarts
myfgkm$totiters
myfgkm$totss

# Use a cluster validation method from package 'clv'.

# real.cluster is the real class label of the data 'fgkm.sample'.
real.cluster <- rep(1:3, each=200)

# std.ext() returns four values SS, SD, DS, DD.
std <- std.ext(as.integer(myfgkm$cluster), real.cluster)

# Rand index 
clv.Rand(std)

# Jaccard coefficient
clv.Jaccard(std)

}
\keyword{
  soft subspace clustering, feature group
}

