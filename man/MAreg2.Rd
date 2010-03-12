\name{MAreg2}
\alias{MAreg2}
\title{Meta-Regression (Multiple Predictor Model) 
}
\description{ Computes multiple predictor fixed or random effects meta-regression (continuous and/or categorical).
}
\usage{
MAreg2(reg)
}
\arguments{
  \item{reg}{ Weighted linear regression (saved as an object, see details below).
}
}
\details{ Running this function requires the researcher to conduct a weighted linear regression while saving it as an object, e.g.:
  
\code{reg <- lm(df$g ~ df$mod8 + df$mod1, weights= df$wi.tau)}

where lm = linear model, df$g = Hedges g, df$mod8 = predictor/moderator, df$mod1 = second predictor/moderator, and df$wi.tau = weighting variable. The outcome variable used in the initial regression should be a Hedges g (Cooper et al., 2009). The predictors (moderator(s)) can be either continuous or categorical. Weight the regression by either the fixed or random effect weight (e.g., fixed equals \code{data$wi} and random equals \code{data$wi.tau}). See \code{ComplData} function to assist in reducing the dataset, if applicable. See Konstantopoulos & Hedges (2009; pp. 289-290) for the computations used in this function.
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
\code{\link{MAregGraph}},
\code{\link{ComplData}}

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

regr <- lm(df$g ~ df$mod, weights= df$wi)
 
MAreg2(regr) # Function will correct inaccurate OLS standard errors
}  
\keyword{ model }
\keyword{ regression }

