\name{FunnelPlot}
\alias{FunnelPlot}
\title{Funnel Plot
}
\description{Outputs a funnel plot from a fixed or random effects omnibus analysis to visually assess meta-analytic data for publication bias.
}
\usage{
FunnelPlot(meta, method = "random", title = NULL)
}
\arguments{
 \item{meta}{\code{data.frame} with id, g (unbiased standardized mean difference), var.g (variance of g), n.1 (sample size of treatment group), and n.2 (sample size for comparison group)for each study. The function requires that the column name for these values be exactly as stated here.
}

  \item{method}{Plot Model used is either \code{random} or \code{fixed} effects. Default is \code{random}.
}
  \item{title}{Plot title. Default is NULL.
}
}
\value{Outputs a funnel plot with omnibus effect size (fixed or random), visual point for each study, where size of point is based on the study's precision (based primarily on sample size), and predicted standard error lines drawn to assist in assessing for publication bias in the meta-analysis. The ggplot2 package outputs the rich graphics. 
}
\references{ Cooper, H., Hedges, L.V., & Valentine, J.C. (2009). \emph{The handbook of research synthesis and meta-analysis} (2nd edition). New York: Russell Sage Foundation. 

Sutton, A. J. (2009). Publication bias. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta analysis} (pp. 435-452). New York: Russell Sage Foundation.
}
\author{AC Del Re & William T. Hoyt

Maintainer: AC Del Re \email{acdelre@gmail.com}
}
\examples{
id<-c(1, rep(1:19))
n.1<-c(10,20,13,22,28,12,12,36,19,12,36,75,33,121,37,14,40,16,14,20)
n.2 <- c(11,22,10,20,25,12,12,36,19,11,34,75,33,120,37,14,40,16,10,21)
g <- c(.68,.56,.23,.64,.49,-.04,1.49,1.33,.58,1.18,-.11,1.27,.26,.40,.49,
.51,.40,.34,.42,1.16)
var.g <- c(.08,.06,.03,.04,.09,.04,.009,.033,.0058,.018,.011,.027,.026,.0040,
.049,.0051,.040,.034,.0042,.016)
mod<-factor(c(rep(c(1,1,2,3),5)))
df<-data.frame(id, n.1,n.2, g, var.g,mod)

# Example

FunnelPlot(df)
}
\seealso{
\code{\link{PubBias}},
}
\keyword{aplot}
