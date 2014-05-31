\name{predict.ewkm}

\alias{predict}
\alias{predict.ewkm}

\title{Predict method for \code{ewkm} model.}

\description{Return the nearest cluster to each observation based on a
  Euclidean distance with each variable weighted differently per cluster.}

\usage{\method{predict}{ewkm}(object, data, \dots)}

\arguments{
  \item{object}{object of class \code{ewkm}.}
  \item{data}{the data that needs to be predicted.  Variables should
  have the same names and order as used in building the
  \code{\link{ewkm}} model.}
  \item{\dots}{other arguments.}
}

\value{a vector of cluster numbers of length \code{nrow(data)}.}

\author{
  Graham Williams (Togaware)
}

\seealso{
  \code{\link{ewkm}}
}
