\name{CorAtten}
\alias{CorAtten}
\title{Correction for Attenuation
}
\description{ Used to correct for attenuated effect sizes due to measurement unreliability.
}
\usage{
CorAtten(meta, xx, yy)
}
\arguments{
  \item{meta}{ \code{data.frame} with g (unbiased standardized mean difference statistic) and var.g (variance of g) for each study.
}
  \item{xx}{ Column for reliability of predictor variable ("independent variable"). 
}
  \item{yy}{ Column for reliability of outcome variable ("dependent variable").
}
}
\value{ \code{data.frame} with a column for standardized mean differences corrected for measurement unreliability (\code{g}), updated standard errors, variance, and study weights based on these corrected values (see Hunter & Schmidt, 2004; pp. 97-98). Studys without reliability information will remain unchanged, as will their standard errors, variances and weights.
}
\author{AC Del Re & William T. Hoyt

Maintainer: AC Del Re \email{acdelre@gmail.com}
}
\references{Hunter, J. E., Schmidt, F. L. (2004). \emph{Methods of meta-analysis} (2nd edition). Thousand Oaks, CA: Sage. 
}
\examples{
# Sample data:

id<-c(1, 1:19)
n<-c(10,20,13,22,28,12,12,36,19,12,36,75,33,121,37,14,40,16,14,20)
g<-c(.68,.56,.23,.64,.49,-.04,.49,.33,.58,.18,-.11,.27,.26,.40,.49,
 .51,.40,.34,.42,.16)
var.g <- c(.08,.06,.03,.04,.09,.04,.009,.033,.0058,.018,.011,.027,.026,.0040,
  .049,.0051,.040,.034,.0042,.016)
xx<-c(.88,.86,.83,.64,.89,.84,.89,.83,.99,.88,.81,.77,.86,.70,.79,
 .71,.80,.74,.82,.86)  # Reliability of "independent variable"
yy<-c(.99,.86,.83,.94,.89,.94,.89,.93,.99,.88,.81,.77,.86,.70,.79,
 .71,.80,.94,.92,.96)  # Reliability of "dependent variable"
   
df<-data.frame(id,n,g, var.g, xx,yy)

# Example        
CorAtten(df,df$xx,df$yy) 
}
\keyword{data}

