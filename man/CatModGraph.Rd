\name{CatModGraph}
\alias{CatModGraph}
\title{Categorical Moderator Graph
}
\description{ Outputs a rich and detailed boxplot graphic for each level of the specified moderator (under a fixed or random effects model). 
}
\usage{
CatModGraph(meta, mod, method = "random", modname = NULL,title = NULL)
}
\arguments{
 \item{meta}{Aggregated \code{data.frame} (see \code{ComplData} function) with id, g (unbiased standardized mean difference), var.g (variance of g), for each study. The function requires that the column name for these values be exactly as stated here.
}
  \item{mod}{ Categorical moderator variable used for analysis. 
}
  \item{method}{ Model used, either \code{random} or \code{fixed} effects model. Default is \code{random}. 
}
  \item{modname}{ Name of moderator to appear on x-axis of plot. Default is NULL. 
}
  \item{title}{ Plot title. Default is NULL. 
}
}
\value{ Boxplot graph with median, mean (denoted by '+' inter-quartile range, max, min, and outliers from a fixed or random effects categorical moderator analysis. Places jitter points (for each study) on the boxplots. The size of each point (representing a study in the analysis) are based on study weights where more precise studies have larger points. The ggplot2 package outputs the rich graphics.
}
\references{ Cooper, H., Hedges, L.V., & Valentine, J.C. (2009). \emph{The handbook of research synthesis and meta-analysis} (2nd edition). New York: Russell Sage Foundation. 
}
\author{AC Del Re & William T. Hoyt

Maintainer: AC Del Re \email{acdelre@gmail.com}
}
\seealso{
\code{\link{CatComp}},
\code{\link{CatMod}},
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

\dontrun{CatModGraph(df,df$mod)}
}   
\keyword{aplot}

