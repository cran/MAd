\name{MAreg1}
\alias{MAreg1}
\title{ Meta-Regression (Single Predictor Model)
}
\description{ Meta-regression function for a single predictor model. This function corrects for inaccuracies in the standard errors produced from an OLS regression model (Konstantopoulos & Hedges, 2009; pp. 289-290).
}
\usage{
MAreg1(meta, mod, method = "random")
}
\arguments{
\item{meta}{\code{data.frame} with id, g (unbiased standardized mean difference), var.g (variance of g), wi (fixed effects study weight), wi.tau (random-effects study weight) for each study. The function requires that the column name for these values be exactly as stated here.
}

  \item{mod}{ Moderator variable used for meta-regression.  
}
  \item{method}{ Model used, either \code{random} or \code{fixed effects}. Default is \code{random}.
}
}
\details{See Konstantopoulos & Hedges (2009; pp. 289-290) for the computations used in this function.
}
\value{
\item{b}{ Unstandardized regression coefficient.
} 
\item{SE}{ Standard error of b coefficient (adjusted SE value).
}
\item{t}{ t-value (adjusted value).
}
\item{CI.Lower}{ Lower 95\% confidence interval.
}
\item{CI.Upper}{ Upper 95\% confidence interval.
}
\item{p}{ Significance level.
}
}
\references{Konstantopoulos & Hedges (2009). Analyzing effect sizes: Fixed-effects models. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta analysis} (pp. 279-293). New York: Russell Sage Foundation.  
}
\author{AC Del Re & William T. Hoyt

Maintainer: AC Del Re \email{acdelre@gmail.com}
}
\seealso{
\code{\link{MAreg1}},
\code{\link{MAreg2}},
\code{\link{MAregGraph}}
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

MAreg1(df, df$mod)   
}
\keyword{ models }
\keyword{ regression }


