\name{OmnibusES}
\alias{OmnibusES}
\title{Omnibus Effect Size (Fixed and Random Effects) 
}
\description{Computes fixed and random effects omnibus effect size for correlations. 
}
\usage{
OmnibusES(meta, var = "weighted", cor = .50)
}
\arguments{
 \item{meta}{\code{data.frame} with id, g (unbiased standardized mean difference), var.g (variance of g), n.1 (sample size of treatment group), and n.2 (sample size for comparison group)for each study. The function requires that the column name for these values be exactly as stated here.
}
  \item{var}{\code{weighted} or \code{unweighted}. Default is \code{weighted}. Use the \code{unweighted} variance method only if Q is rejected and is very large relative to the number of studies in the meta-analysis. 
}
\item{cor}{Estimated correlation among within-study outcome variables. Default is .50 based on Wampold et al. recommended procedures (1997).
}
}
\value{
Fixed and random effects:
 
\item{k}{ Number of studies in the meta-analysis.
}
\item{ES}{ Overall Omnibus effect size for the meta-analysis.
}
\item{var.ES}{ Variance of Omnibus effect size.
}
\item{SE}{ Square root of variance.
}
\item{CI.lower}{ Lower 95\% confidence interval.
}
\item{CI.upper}{ Upper 95\% confidence interval.
}
\item{Z}{ Z-score (standardized value).
}
\item{p}{ Significance level.
}
\item{Q}{ Q-statistic (measure of homogeneity).
}
\item{df}{ Degrees of freedom for Q-statistic.
}
\item{p.hetero}{ Q-statistic p-value (assesses overall heterogeneity between studies).
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

OmnibusES(df)
}
\keyword{ models }

