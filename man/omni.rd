\name{omni}
\alias{omni}
\title{Omnibus Effect Size (Fixed and Random Effects) 
}
\description{Computes fixed and random effects omnibus effect size for correlations. 
}
\usage{
omni(g, var, data, type="weighted", method = "random")
}
\arguments{
\item{g}{Hedges g (unbiased estimate of d) effect size.
}
\item{var}{Vaiance of g.
}
 
  \item{type}{\code{weighted} or \code{unweighted}. Default is \code{weighted}. Use the \code{unweighted} variance method only if Q is rejected and is very large relative to the number of studies in the meta-analysis. 
}
 \item{method}{ Default is \code{random}. For fixed effects, use \code{fixed}. 
}
\item{data}{\code{data.frame} with above values.
}
}
\value{
Fixed and random effects:
 
\item{k}{ Number of studies in the meta-analysis.
}
\item{estimate}{ Unstandardized regression coefficient estimate.
} 
\item{se}{ Standard error of the estimate coefficient.
}
\item{z}{ z-value.
}
\item{ci.l}{ Lower 95\% confidence interval.
}
\item{ci.u}{ Upper 95\% confidence interval.
}
\item{p}{ Significance level.
}
\item{Q}{ Q-statistic (measure of homogeneity).
}
\item{df.Q}{ Degrees of freedom for Q-statistic.
}
\item{Qp}{ Q-statistic p-value (assesses overall homogeneity between studies).
}
\item{I2}{ Proportion of total variation in effect size that is due to heterogeneity rather than chance (see Shadish & Haddock, 2009; pp. 263).
}
}
\references{ Shadish & Haddock (2009). Analyzing effect sizes: Fixed-effects models. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta analysis} (pp. 257-278). New York: Russell Sage Foundation. 
}
\author{AC Del Re & William T. Hoyt

Maintainer: AC Del Re \email{acdelre@gmail.com}
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

omni(g = g, var = var.g, data = df, type="weighted", method = "random")
}
\keyword{ models }

