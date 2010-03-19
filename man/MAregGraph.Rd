\name{MAregGraph}
\alias{MAregGraph}
\title{Meta Regression Scatterplot
}
\description{Outputs a scatterplot from a fixed or random effects meta regression (continuous and/or categorical). 
}
\usage{
MAregGraph(meta, mod, method = "random", modname = NULL, 
title = NULL, ylim = c(0, 1))
}
\arguments{
 \item{meta}{Aggregated \code{data.frame} (see \code{ComplData} function) with id, g (unbiased standardized mean difference), var.g (variance of g), for each study. The function requires that the column name for these values be exactly as stated here.
}
  \item{mod}{ Moderator variable used for meta regression.  
}
  \item{method}{ Model used, either \code{random} or \code{fixed} effects. Default is \code{random}.
}
  \item{modname}{Name of moderator to appear on x axis of plot. Default is NULL.
}
  \item{title}{Plot title. Default is NULL.
}
  \item{ylim}{Limits of y-axis with the first argument for the minimum y-value and the second for the maximum y value. Default is \code{c(0, 1)}.
}
}
\value{ Scatterplot with fixed or random effects regression line with size of visual points based on study weights, where the more precise studies have larger points. The ggplot2 package outputs the rich graphics.
}
\references{ Cooper, H., Hedges, L.V., & Valentine, J.C. (2009). \emph{The handbook of research synthesis and meta analysis} (2nd edition). New York: Russell Sage Foundation. 
}
\author{AC Del Re & William T. Hoyt

Maintainer: AC Del Re \email{acdelre@gmail.com}
}
\seealso{
\code{\link{MAreg1}},
\code{\link{MAreg2}}
}
\examples{
id<-c(1, 1:19)
n.1<-c(10,20,13,22,28,12,12,36,19,12,36,75,33,121,37,14,40,16,14,20)
n.2 <- c(11,22,10,20,25,12,12,36,19,11,34,75,33,120,37,14,40,16,10,21)
g <- c(.68,.56,.23,.64,.49,-.04,1.49,1.33,.58,1.18,-.11,1.27,.26,.40,.49,
.51,.40,.34,.42,1.16)
var.g <- c(.08,.06,.03,.04,.09,.04,.009,.033,.0058,.018,.011,.027,.026,.0040,
.049,.0051,.040,.034,.0042,.016)
mod<-factor(c(rep(c(1,1,2,3),5)))
df<-data.frame(id, n.1,n.2, g, var.g,mod)

# Example

\dontrun{MAregGraph(df, df$mod)   }
}
\keyword{ aplot }

