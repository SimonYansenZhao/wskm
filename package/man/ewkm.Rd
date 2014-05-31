\name{ewkm}
\alias{ewkm}
\title{
  
  Entropy Weighted K-Means

}
\description{

  Perform an entropy weighted subspace k-means.

}
\usage{

 ewkm(x, k, lambda=1, maxiter=100, delta=0.00001, maxrestart=10)

}
\arguments{

  \item{x}{numeric matrix of observations and variables.}

  \item{k}{target number of clusters.}

  \item{lambda}{parameter for variable weight distribution.}

  \item{maxiter}{maximum number of iterations.}
    
  \item{delta}{maximum change allowed between iterations for convergence.}

  \item{maxrestart}{maximum number of restarts. Default is 10 so that we
  stand a good chance of getting a full set of clusters. Normally, any
  empty clusters that result are removed from the result, and so we may
  obtain fewer than k clusters if we don't allow restarts (i.e.,
  maxrestart=0). If < 0 then there is no limit on the number of restarts
  and we are much more likely to get a full set of k clusters.}

}
\details{
  
  The entopy weighted k-means clustering algorithm is a subspace
  clusterer ideal for high dimensional data. Along with each cluster we
  also obtain variable weights that provide a relative measure of the
  importance of each variable to that cluster.

  The algorithm is based on the k-means approach to clustering. An
  initial set of k means are identified as the starting
  centroids. Observartions are clustered to the nearest centroid
  according to a distance measure. This defines the initial
  clustrering. New centroids are then identified based on these
  clusters.
  
  Weights are then calculated for each variable within each cluster,
  based on the current clustering.  The weights are a measure of the
  relative importance of each variable with regard to the membership of
  the observations to that cluster. These weights are then incorporated
  into the distance function, typically reducing the distance for the
  more important variables.

  New centroids are then calculated, and using the weighted distance
  measure each observation is once again clustered to its nearest
  centroid.

  The process continues until convergence (using a measure of dispersion
  and stopping when the change becomes less than delta) or until a
  specified number of iterations has been reached (maxiter).

  Large lambda (e,g, > 3) lead to a relatively even distribution of
  weights across the variables. Small lambda (e.g., < 1) lead to a more
  uneven distribution of weights, giving more discrimintation between
  features. Recommended values are between 1 and 3.

  Always check the number of iterations, the number of restarts, and the
  total number of iterations as they give a good indication of whether
  the algorithm converged.

  As with any distance based algorithm, be sure to rescale your numeric
  data so that large values do not bias the clustering. A quick
  rescaling method to use is \code{\link{scale}}.
  
}
\value{
  
  Returns an object of class "kmeans" and "ewkm", compatible with other
  functions that work with kmeans objects, such as the 'print'
  method. The object is a list with the following components in addition
  to the components of the kmeans object:

  weights: A matrix of weights recording the relative importance of each
  variable for each cluster.

  iterations: This reports on the number of iterations before
  termination. Check this to see whether the maxiters was reached. If so
  then the algroithm may not be converging,and thus the resulting
  clustering may not be particularly good.

  restarts: The number of times the clustering restarted because of a
  disappearing cluster resulting from one or more k-means having no
  observations associated with it. An number here greater than 0
  indicates that the algorithm is not converging on a clustering for the
  given k. It is recommended that k be reduced.

  total.iterations: The total number of iterations over all restarts.

}
\references{

  Jing Liping, M. K. Ng, J. Z. Huang, An Entropy Weighting k-Means
  Algorithm for Subspace Clustering of High-Dimensional Sparse Data,
  IEEE Transactions on Knowledge and Data Engineering, 19 (8), Aug 2007
  pp 1026--1041.
  
}
\author{

  Qiang Wang, Xiaojun Chen, Graham J Williams, Joshua Z Huang
  
}
\seealso{

  \code{\link{plot.ewkm}}.

}
\examples{

myewkm <- ewkm(iris[1:4], k=3, lambda=0.5, maxiter=100)

plot(iris[1:4], col=myewkm$cluster)

# For comparative testing

mykm <- kmeans(iris[1:4], 3)

plot(iris[1:4], col=mykm$cluster)

}
