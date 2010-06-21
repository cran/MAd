\name{plotcon}
\alias{plotcon}
\title{Meta Regression Scatterplot
}
\description{Outputs a scatterplot from a fixed or random effects meta regression (continuous and/or categorical). 
}
\usage{
plotcon(g, var, mod, data, method= "random", modname=NULL, 
  title=NULL, ylim=c(0, 1), ...)
}
\arguments{
 \item{g}{Hedges g (unbiased estimate of d) effect size.
}
\item{var}{Vaiance of g.
}
  \item{mod}{Categorical moderator variable used for moderator analysis.
} 
 \item{method}{ Default is \code{random} (Restricted-Maximal Likelihood), which is the standard random effects method. For fixed effects, use \code{fixed}. 
}
\item{data}{\code{data.frame} with values above.
}
  \item{modname}{Name of moderator to appear on x axis of plot. Default is NULL.
}
  \item{title}{Plot title. Default is NULL.
}
  \item{ylim}{Limits of y-axis with the first argument for the minimum y-value and the second for the maximum y value. Default is \code{c(0, 1)}.
}
\item{...}{ Additional arguments to be passed to ggplot.
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
\code{\link{mareg}},
\code{\link{plotcat}}
}
\examples{
id<-c(1, 1:19)
n.1<-c(10,20,13,22,28,12,12,36,19,12,36,75,33,121,37,14,40,16,14,20)
n.2 <- c(11,22,10,20,25,12,12,36,19,11,34,75,33,120,37,14,40,16,10,21)
g <- c(.68,.56,.23,.64,.49,-.04,1.49,1.33,.58,1.18,-.11,1.27,.26,.40,.49,
.51,.40,.34,.42,1.16)
var.g <- c(.08,.06,.03,.04,.09,.04,.009,.033,.0058,.018,.011,.027,.026,.0040,
.049,.0051,.040,.034,.0042,.016)
mod<-rep(c(1:4),5)
df<-data.frame(id, n.1,n.2, g, var.g,mod)

# Example

\dontrun{plotcon(g = g, var = var.g, mod = mod, data = df, method= "random", 
modname= "Moderator") }
}
\keyword{ aplot }

