\name{weights}
\title{Weights added to Meta Data}
\alias{weights}
\description{Adds weights to the meta-analysis data set.}
\usage{
weights(g, var.g, data)
}
\arguments{
\item{g}{Hedges g (unbiased standardized mean difference estimate). 
}
\item{var.g}{Variance of g. 
}
\item{data}{\code{data.frame} with values above. 
}
}
\value{Adds fixed and random-effects weights and confidence intervals to meta data set.
}
\author{ AC Del Re & William T. Hoyt

Maintainer: AC Del Re \email{acdelre@gmail.com}
}
\keyword{data}