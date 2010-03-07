\name{wifun}
\title{Weights added to Meta Data}
\alias{Wifun}
\description{Adds weights to the meta-analysis data set.}
\usage{
Wifun(meta)
}
\arguments{
\item{meta}{\code{data.frame} with id, g (unbiased standardized mean difference), var.g (variance of g), n.1 (sample size of treatment group), and n.2 (sample size for comparison group)for each study. The function requires that the column name for these values be exactly as stated here.
}
}
\value{Adds fixed and random-effects weights and confidence intervals to meta data set.
}
\author{ AC Del Re & William T. Hoyt

Maintainer: AC Del Re \email{acdelre@gmail.com}
}
\keyword{data}