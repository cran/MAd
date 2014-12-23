\name{mareg}
\alias{mareg}
\title{ Meta-Regression 
}
\description{ Meta-regression function for a single or multiple predictor model. This function is a wrapper for the \code{rma()} function in the metafor package (Viechtbauer, W, 2010). Please see http://CRAN.R-project.org/package=metafor for details or for more advanced functionality with the \code{rma()} function. 
}
\usage{
mareg(formula, var, data, method = "REML", subset,  digits = 3, ...)
}
\arguments{
  \item{formula}{ This is a formula based function, similar to other functions in R (e.g., lm), where the criterion variable (e.g., Hedges g in this case) is dependent on ('~') the predictor variables (e.g., moderators). The formula for two moderators would take this form: mareg(g ~ mod1 + mod2, var.g, data), where g is the criterion variable predicted by mod1 and mod2. The variance (var) of each g is var.g in this case.   
}
  \item{var}{ Variance of g.
  }
\item{data}{Aggregated \code{data.frame} (see \code{agg} function for setting up the dataset for these analyses) with id, g (unbiased standardized mean difference), var.g (variance of g) for each study. 
}
  \item{method}{ Default is \code{REML} (Restricted-Maximal Likelihood), which is the standard random effects method. For fixed effects, use \code{FE}. Other options are specified in the \code{metafor} package manual ('rma' function).
}
 \item{subset}{ an optional vector specifying a subset of observations to be used in the fitting process.
  }
   \item{digits}{ Number of digits to output. Default is 3.
  } 
 \item{...}{ Additional arguments to be passed to rma().
  }  
}
\details{See Wolfgang Viechtbauer (2010). metafor: Meta-Analysis Package for
  R. R package version 1.1-0. for the details  of the \code{rma()}  function. http://CRAN.R-project.org/package=metafor
}
\value{
\item{estimate}{Meta-regression coefficient estimate.
} 
\item{se}{ Standard error of the estimate coefficient.
}
\item{z}{ z-value.
}
\item{ci.l}{ Lower 95\% confidence interval.
}
\item{ci.u}{ Upper 95\% confidence interval.
}
\item{p}{ p-value (significance level).
}
\item{QE}{ Q-error. Measure of error in the model.
}
\item{QE.df}{ Degrees of freedom for Q-error.
}
\item{QEp}{ Q-error p-value (for homogeneity).
}
\item{QM}{ Q-model. Measure of model fit.
}
\item{QM.df}{ Degrees of freedom for Q-model.
}
\item{QMp}{ Q-between p-value (for homogeneity). QM and QMp provide the test of whether the moderator variable(s) account for significant variance among effect sizes.
} 
}
\references{Wolfgang Viechtbauer (2010). metafor: Meta-Analysis Package for
  R. R package version 1.1-0. http://CRAN.R-project.org/package=metafor
}
\seealso{
\code{\link{wd}},
\code{\link{plotcon}}
}
\examples{
# install metafor
# install.packages('metafor', dependencies = TRUE)

# Sample data
data(dat.sim.final)


# Examples

# OMNIBUS
m0 <- mareg(es~1, var=var, data=dat.sim.final)
summary(m0)

# META-REGRESSION
m1 <- mareg(es~dose, var=var, data=dat.sim.final)
summary(m1)  # SINGLE MODERATOR

m2 <- mareg(es~stress, var=var, data=dat.sim.final)
summary(m2)  # SINGLE MODERATOR

m3 <- mareg(es~dose + stress, var=var, data=dat.sim.final)
summary(m3)  # MULTIPLE MODERATOR  
}



