\name{MRfit}
\alias{MRfit}
\title{Meta-Regression Model Fit 
}
\description{ Computes the model fit of a moderator analysis comparing an empty model to one with x number of predictors.
}
\usage{
MRfit(...)
}
\arguments{
  \item{...}{ Any number of weighted linear regressions (saved as objects, see details in example section below).
}
}
\details{ Running this function requires the researcher to conduct a weighted hierarchical linear regression while saving each block as an object (see example below).
The outcome variable is Hedges g (unbiased standardized mean difference) and the predictor moderator(s) can be either continuous or categorical. Weight the regression by either the fixed or random effect weight (e.g., fixed equals \code{df$wi} and random equals \code{df$wi.tau}). See \code{ComplData} function to assist in reducing the dataset, if applicable. See Konstantopoulos & Hedges (2009; pp. 289-290) for the computations used in this function.
}
\value{
\item{df.Q}{ Q degrees of freedom.
} 
\item{Qe}{ Q-error, weighted residual sum of squares.
}
\item{predictors}{ Number of predictors in the model for each block.
}
\item{Qb}{ Q-between, indicating the between-group heterogeneity.
}
\item{F}{ F-value.
}
\item{Pvalue}{ Significance level.
}
\item{R^2}{ R-squared, the coefficient of determination, is the proportion of variability in a data set that is accounted for by the statistical model.
}
\item{R^2.change}{ Change in R-squared from one model to the next (from one block to the next).

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

reg0 <- lm(df$g ~ 1, weights= df$wi.tau)  # empty model
reg1 <- lm(df$g ~ df$mod, weights= df$wi.tau)  # 1 predictor model
MRfit(reg0, reg1) 
}  
\keyword{ model }
\keyword{ regression }

