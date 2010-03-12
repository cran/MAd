\name{MeanDiffg}
\alias{MeanDiffg}
\title{Automated Conversion to Standardized Mean Difference
}
\description{Adds d (standardized mean difference) to a \code{data.frame}. Required inputs are: n.1 (sample size of group one),  sd.1 (standard deviation of group one), n.2 (sample size of group two). 
}
\usage{
MeanDiffg(meta, d = meta$d, var.d = meta$var.d, n.1= meta$n.1, 
          n.2 =  meta$n.2)
}
\arguments{
  \item{meta}{\code{data.frame} with d (Standardized mean difference), var.d (Variance of d), n.1 (sample size of group one),n.2 (sample size of group two). 
}
\item{d}{Standardized mean difference (biased).
}
  \item{var.d}{Variance of d.
}
  \item{n.1}{sample size of group one.
}
 \item{n.2}{sample size of group two.
}
}
\value{
  \item{g }{Unbiased standardized mean difference.
  }
  \item{var.g }{Variance of g.
  }
  \item{se.g}{Standard error of g.
  }
}
\author{AC Del Re & William T. Hoyt

Maintainer: AC Del Re \email{acdelre@gmail.com}
}
\references{Borenstein (2009). Effect sizes for continuous data. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta analysis} (pp. 279-293). New York: Russell Sage Foundation.
}
\seealso{
\code{\link{MeanDiff}},
\code{\link{MeanDiffd}}
}
\examples{
id<-c(1, 1:19)
n.1<-c(10,20,13,22,28,12,12,36,19,12,36,75,33,121,37,14,40,16,14,20)
n.2 <- c(11,22,10,20,25,12,12,36,19,11,34,75,33,120,37,14,40,16,10,21)
d <- c(.68,.56,.23,.64,.49,-.04,1.49,1.33,.58,1.18,-.11,1.27,.26,.40,.49,
.51,.40,.34,.42,1.16)
var.d <- c(.08,.06,.03,.04,.09,.04,.009,.033,.0058,.018,.011,.027,.026,.0040,
.049,.0051,.040,.034,.0042,.016)
mod<-factor(c(rep(c(1,1,2,3),5)))
df<-data.frame(id, n.1,n.2, d, var.d,mod)

MeanDiffg(df)
}
\keyword{data}

