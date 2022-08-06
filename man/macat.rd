\name{macat}
\alias{macat}
\title{Categorical Moderator Analysis 
}
\description{Computes single predictor categorical moderator analysis under a fixed or random effects model.
}
\usage{
macat(g, var, mod, data, method= "random")
}
\arguments{
 \item{g}{Hedges g (unbiased estimate of d) effect size.
}
\item{var}{Vaiance of g.
}
  \item{mod}{Categorical moderator variable used for moderator analysis.
} 
 \item{method}{ Default is \code{random}. For fixed effects, use \code{fixed}. 
}

 \item{data}{\code{data.frame} with values above.
}
}
\details{See Konstantopoulos & Hedges (2009; pp. 280-288) for the computations used in this function.
}
\value{

\item{mod}{ Level of the categorical moderator.
} 
\item{k}{ Number of studies for each level of the moderator.
}
\item{estimate}{ Mean effect size of each level of the moderator.
}
\item{ci.l}{ Lower 95\% confidence interval.
}
\item{ci.u}{ Upper 95\% confidence interval.
}
\item{z}{ z-score (standardized value).
}
\item{p}{ Significance level.
}
\item{var}{ Variance of effect size.
}
\item{se}{ Square root of variance.
}
\item{Q}{ Q-statistic (measure of homogeneity).
}
\item{df}{ Degrees of freedom for Q-statistic.
}
\item{p.h}{p-value for homogeneity within that level of the moderator. 
}
\item{I2}{ Proportion of total variation in effect size that is due to heterogeneity rather than chance (see Shadish & Haddock, 2009; pp. 263).
}
\item{Q}{ Q-statistic overall. Note: Whether fixed or random effects analyses are conducted, the Q statistic reported is for the fixed effect model. Therefore, Qb + Qw != Q in the random effects output.
}  
\item{Qw}{ Q-within (or error). Measure of within-group heterogeneity.
}
\item{Qw.df}{ Degrees of freedom for Q-within.
}
\item{Qw.p}{ Q-within p-value (for homogeneity).
}
\item{Qb}{ Q-between (or model). Measure of model fit.
}
\item{Qb.df}{ Degrees of freedom for Q-between.
}
\item{Qb.p}{ Q-between p-value (for homogeneity). Qb and Qb.p provide the test of whether the moderator variable(s) account for significant variance among effect sizes.
} 
}
\references{Konstantopoulos & Hedges (2009). Analyzing effect sizes: Fixed-effects models. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta analysis} (pp. 279-293). New York: Russell Sage Foundation.  

Shadish & Haddock (2009). Analyzing effect sizes: Fixed-effects models. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta analysis} (pp. 257-278). New York: Russell Sage Foundation. 
}
\author{AC Del Re & William T. Hoyt

Maintainer: AC Del Re \email{acdelre@gmail.com}
}
\seealso{
\code{\link{plotcat}}
}
\examples{
id<-c(1:20)
n.1<-c(10,20,13,22,28,12,12,36,19,12,36,75,33,121,37,14,40,16,14,20)
n.2 <- c(11,22,10,20,25,12,12,36,19,11,34,75,33,120,37,14,40,16,10,21)
g <- c(.68,.56,.23,.64,.49,-.04,1.49,1.33,.58,1.18,-.11,1.27,.26,.40,.49,
.51,.40,.34,.42,1.16)
var.g <- c(.08,.06,.03,.04,.09,.04,.009,.033,.0058,.018,.011,.027,.026,.0040,
.049,.0051,.040,.034,.0042,.016)
mod<-factor(c(rep(c(1,1,2,3),5)))
df<-data.frame(id, n.1,n.2, g, var.g,mod)

# Example

# Random effects
macat(g = g, var= var.g, mod = mod, data = df, method= "random")
}
\keyword{models}
