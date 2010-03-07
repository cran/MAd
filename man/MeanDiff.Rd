\name{MeanDiff}
\alias{MeanDiff}
\title{Automated Conversion to Unbiased Standardized Mean Difference
}
\description{Adds d (standardized mean differnce), var.d (variance of d), g (unbiased standardized mean difference), var.g (variance of g) to a \code{data.frame}. Required inputs are: n.1 (sample size of group one), m.1 (raw post-mean value of group one), sd.1 (standard deviation of group one), n.2 (sample size of group two), m.2 (raw post-mean value of group two), sd.2 (standard deviation of group two). 
}
\usage{
MeanDiff(meta, n.1 = meta$n.1, m.1 = meta$m.1, sd.1 = meta$sd.1, n.2 = meta$n.2, m.2 = meta$m.2, sd.2 = meta$sd.2, denom = "pooled.sd")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{meta}{\code{data.frame} with n.1 (sample size of group one), m.1 (raw post-mean value of group one), sd.1 (standard deviation of group one), n.2 (sample size of group two), m.2 (raw post-mean value of group two), sd.2 (standard deviation of group two). 
}
  \item{n.1}{sample size of group one.
}
  \item{m.1}{raw post-mean value of group one.
}
  \item{sd.1}{standard deviation of group one.
}
  \item{n.2}{sample size of group two.
}
  \item{m.2}{raw post-mean value of group two.
}
  \item{sd.2}{standard deviation of group two.
}
  \item{denom}{Value in the denominator to standardize mean difference. \code{pooled.sd} pools both groups in deriving d. \code{control.sd} uses the standard deviation of group two to calculate d. 
}
}
\value{
  \item{d }{Standardized mean difference}
  \item{var.d }{Variance of d}
  \item{se.d }{Standard error of d}
  \item{dg}{Unbiased standardized mean difference}
  \item{var.g }{Variance of g}
  \item{se.g }{Standard error of g}
}
\author{AC Del Re & William T. Hoyt

Maintainer: AC Del Re \email{acdelre@gmail.com}
}
\references{Borenstein (2009). Effect sizes for continuous data. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta analysis} (pp. 279-293). New York: Russell Sage Foundation.
}
\seealso{
\code{\link{MeanDiffd}},
\code{\link{MeanDiffg}}
}
\examples{
id<-rep(1:20)
n.1<-c(10,20,13,22,28,12,12,36,19,12,36,75,33,121,37,14,40,16,14,20)
n.2 <- c(11,22,10,20,25,12,12,36,19,11,34,75,33,120,37,14,40,16,10,21)
m.1 <- c(1.68,.56,.23,.64,.49,.4,1.49,.53,.58,1.18,.11,1.27,.26,.40,.49,
.51,.40,.34,.42,.66)
m.2 <- c(1.38,.36,.23,.34,.29,.4,1.9,.33,.28,1.1,.111,.27,.21,.140,.149,
.51,.140,.134,.42,.16)
sd.1 <- c(.28,.26,.23,.44,.49,.34,.39,.33,.58,.38,.31,.27,.26,.40,
.49,.51,.140,.134,.42,.46)
sd.2 <- c(.28,.26,.23,.44,.49,.44,.39,.33,.58,.38,.51,.27,.26,.40,
.49,.51,.140,.134,.142,.36)
mod1<-c(1,2,3,4,1,2,8,7,5,3,9,7,5,4,3,2,3,5,7,1)
mod2<-factor(c(rep(c(1,2,3,4),5)))
df<-data.frame(id, n.1,m.1, sd.1, n.2, m.2, sd.2,mod1, mod2)


MeanDiff(df)
}
\keyword{data}