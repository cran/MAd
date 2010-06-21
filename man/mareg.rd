\name{mareg}
\alias{mareg}
\title{ Meta-Regression 
}
\description{ Meta-regression function for a single or multiple predictor model. This function is a wrapper for the \code{rma()} function in the metafor package (Viechtbauer, W, 2010). Please see http://CRAN.R-project.org/package=metafor for details or for more advanced functionality with the \code{rma()} function. 
}
\usage{
mareg(formula, var, data, method = "REML", subset,  ...)
}
\arguments{
  \item{formula}{ This is a formula based function, similar to other functions in R (e.g., lm), where the criterion variable (e.g., Hedges g in this case) is dependent on ('~') the predictor variables (e.g., moderators). The formula for two moderators would take this form: mareg(g ~ mod1 + mod2, var.g, data), where g is the criterion variable predicted by mod1 and mod2. The variance (var) of each g is var.g in this case.   
}
  \item{var}{ Variance of g.
  }
\item{data}{Aggregated \code{data.frame} (see \code{ComplData} function for setting up the dataset for these analyses) with id, g (unbiased standardized mean difference), var.g (variance of g) for each study. 
}
  \item{method}{ Default is \code{REML} (Restricted-Maximal Likelihood), which is the standard random effects method. For fixed effects, use \code{FE}. Other options are specified in the \code{metafor} package manual ('rma' function).
}
 \item{subset}{ an optional vector specifying a subset of observations to be used in the fitting process.
  }
 \item{...}{ Additional arguments to be passed to rma().
  }  
}
\details{See Wolfgang Viechtbauer (2010). metafor: Meta-Analysis Package for
  R. R package version 1.1-0. for the details  of the \code{rma()}  function. http://CRAN.R-project.org/package=metafor
}
\value{
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
}
\references{Wolfgang Viechtbauer (2010). metafor: Meta-Analysis Package for
  R. R package version 1.1-0. http://CRAN.R-project.org/package=metafor
}
\seealso{
\code{\link{wd}}
\code{\link{plotcon}}
}
\examples{
# install metafor
# install.packages('metafor', dependencies = TRUE)

# Sample data
id<-c(1:20)
n.1<-c(10,20,13,22,28,12,12,36,19,12,36,75,33,121,37,14,40,16,14,20)
n.2 <- c(11,22,10,20,25,12,12,36,19,11,34,75,33,120,37,14,40,16,10,21)
g <- c(.68,.56,.23,.64,.49,-.04,1.49,1.33,.58,1.18,-.11,1.27,.26,.40,.49,
.51,.40,.34,.42,1.16)
var.g <- c(.08,.06,.03,.04,.09,.04,.009,.033,.0058,.018,.011,.027,.026,.0040,
.049,.0051,.040,.034,.0042,.016)
mod<-factor(c(rep(c(1,1,2,3),5)))
mods2<-c(rep(1:5,4))
df<-data.frame(id, n.1,n.2, g, var.g,mod, mods2)


# Examples

# Random Effects
mareg(g~ mod + mods2, var = var.g, method = "REML", data = df)

# Fixed Effects
mareg(g~ mod + mods2, var = var.g, method = "FE", data = df)   
}



