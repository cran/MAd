\name{CatComp}
\alias{CatComp}
\title{Direct Categorical Moderator Comparison 
}
\description{Function for a planned comparison between two levels of a moderator under a fixed and random effects model.
}
\usage{
CatComp(meta, mod, x1, x2, method = "post.hoc1")
}
\arguments{
  \item{meta}{\code{data.frame} with id (study), g (unbiased standardized mean difference), var.g (variance of g), for each study. The function requires that the column name for these values be exactly as stated here.
}
  \item{mod}{Categorical moderator variable used for moderator analysis.
}
  \item{x1}{One level of categorical moderator.
}
  \item{x2}{Comparison level of same categorical moderator.
}
  \item{method}{\code{post.hoc1} assumes the comparison was not planned prior to conducting the meta analysis. The a priori option, \code{planned}, assumes the researcher planned to conduct the analysis a priori. Default is \code{post.hoc1} using Tukey's Honest Significant Differences method. \code{post.hoc2} uses the Scheffe post hoc statistical method. 
}
}
\details{See Konstantopoulos & Hedges (2009; pp. 280-288) for the computations used in this function.
}
\value{
\item{d}{Standardized mean difference (d).
}
\item{Var(d)}{Variance of d (standardized mean difference).
}
\item{CI.lower}{ Lower 95\% confidence interval.
}
\item{CI.lower}{ Upper 95\% confidence interval.
}
}
\references{ Konstantopoulos & Hedges (2009). Analyzing effect sizes: Fixed-effects models. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta analysis} (pp. 279-293). New York: Russell Sage Foundation. 
} 
\author{AC Del Re & William T. Hoyt

Maintainer: AC Del Re \email{acdelre@gmail.com}
}
\seealso{
\code{\link{CatMod}},
\code{\link{CatModGraph}}
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

CatComp(df,df$mod,1,2, method = "post.hoc2") # compares level 1 vs level 2 of moderator. 
}
\keyword{models}

