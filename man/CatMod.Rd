\name{CatMod}
\alias{CatMod}
\title{Categorical Moderator Analysis 
}
\description{Computes single predictor categorical moderator analysis under a fixed and Random effects model.
}
\usage{
CatMod(meta, mod)
}
\arguments{
  \item{meta}{Aggregated \code{data.frame} (see \code{ComplData} function) with id, g (unbiased standardized mean difference), var.g (variance of g), for each study. The function requires that the column name for these values be exactly as stated here.
}
  \item{mod}{Categorical moderator variable used for moderator analysis.
}
}
\details{See Konstantopoulos & Hedges (2009; pp. 280-288) for the computations used in this function.
}
\value{

\item{Mod}{ Level of the categorical moderator.
} 
\item{k}{ Number of studies for each level of the moderator.
}
\item{ES}{ Mean effect size of each level of the moderator.
}
\item{CI.lower}{ Lower 95\% confidence interval.
}
\item{CI.upper}{ Upper 95\% confidence interval.
}
\item{Z}{ z-score (standardized value).
}
\item{p}{ Significance level.
}
\item{var.ES}{ Variance of effect size.
}
\item{SE}{ Square root of variance.
}
\item{Q}{ Q-statistic (measure of homogeneity).
}
\item{df}{ Degrees of freedom for Q-statistic.
}
\item{p.h}{p-value for homogeneity within that level of the moderator. 
}
\item{I2}{ Proportion of total variation in effect size that is due to heterogeneity rather than chance (see Shadish & Haddock, 2009; pp. 263).
}
\item{Q}{ Q-statistic (measure of homogeneity).
} 
\item{Qw}{ Q-within (measure of within-group heterogeneity).
}
\item{df.w}{ Degrees of freedom for Q-within.
}
\item{p.w}{ Q-within p-value (for homogeneity).
}
\item{Qb}{ Q-between (measure of between-group heterogeneity)
}
\item{df.b}{ Degrees of freedom for Q-between.
}
\item{p.b}{ Q-between p-value (for homogeneity)
} 
}
\references{Konstantopoulos & Hedges (2009). Analyzing effect sizes: Fixed-effects models. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta analysis} (pp. 279-293). New York: Russell Sage Foundation.  

Shadish & Haddock (2009). Analyzing effect sizes: Fixed-effects models. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta analysis} (pp. 257-278). New York: Russell Sage Foundation. 
}
\author{AC Del Re & William T. Hoyt

Maintainer: AC Del Re \email{acdelre@gmail.com}
}
\seealso{
\code{\link{CatComp}},
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

CatMod(df,df$mod)   
}
\keyword{models}
