\name{plotcon}
\alias{plotcon}
\title{Meta Regression Scatterplot
}
\description{Outputs a scatterplot from a fixed or random effects meta regression (continuous and/or categorical). 
}
\usage{
plotcon(g, var, mod, data, method= "random", modname=NULL, 
  title=NULL, ...)
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
# SAMPLE DATA
MA2 <-read.table(textConnection("
id       es1       var1 n.1 n.2 mod1 mod2
1   1 0.5695938 0.04906967  26  30    a   20
2   2 0.4123667 0.04362541  28  34    b   30
3   3 0.4084333 0.04458363  34  28    a   25
4   4 0.5014756 0.04186354  37  29    b   35
5   5 0.5540745 0.04339382  31  32    b   40
6   1 0.5695938 0.04906967  26  30    a   20
7   2 0.4123667 0.04362541  28  34    b   30
8   3 0.4084333 0.04458363  34  28    a   25
9   4 0.5014756 0.04186354  37  29    b   35
10  5 0.5540745 0.04339382  31  32    b   40"))


# EXAMPLE
plotcon(es1, var1, mod2, data=MA2, method= "fixed", modname="NULL",title="NULL")
}
\keyword{ aplot }

