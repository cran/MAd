\name{atten}
\alias{atten}
\title{Correction for Attenuation
}
\description{ Used to correct for attenuated effect sizes due to measurement unreliability.
}
\usage{
atten(g, xx, yy, data)
}
\arguments{
 \item{g}{Hedges g (unbiased estimate of d) effect size.
} 
  \item{xx}{ Column for reliability of predictor variable ("independent variable"). 
}
  \item{yy}{ Column for reliability of outcome variable ("dependent variable").
  }
  \item{data}{ \code{data.frame} with the above values.
}
}
\value{A new column for g corrected for attenuation (\code{g.corrected}) will be added to the data, for those xx & yy columns with complete data. 
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
atten(g= g, xx = xx, yy = yy, data= df) 
}
\keyword{data}

