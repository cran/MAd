\name{ancova_to_d2}
\alias{ancova_to_d2}
\title{ANCOVA F-statistic to d II
}
\description{
Converts an ANCOVA F-statistic with pooled standard deviation to d (standardized mean difference)
}
\usage{
ancova_to_d2(m.1.adj, m.2.adj, s.pooled, n.1, n.2, R, q)
}
\arguments{
  \item{m.1.adj}{Adjusted mean of treatment group from ANCOVA.
}
  \item{m.2.adj}{Adjusted mean of comparison group from ANCOVA.
}
  \item{s.pooled}{Pooled standard deviation.
}
  \item{n.1}{Treatment group sample size.
}
  \item{n.2}{Comparison group sample size.
}
  \item{R}{Covariate outcome correlation or multiple correlation.
}
  \item{q}{number of covariates.
}
}
\value{
\item{d}{Standardized mean difference (d).}
\item{var_d }{Variance of d.}

}
\author{ AC Del Re & William T. Hoyt

Maintainer: AC Del Re \email{acdelre@gmail.com}
}
\references{Borenstein (2009). Effect sizes for continuous data. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta analysis} (pp. 279-293). New York: Russell Sage Foundation.
}
\seealso{
\code{\link{d_to_g}},
\code{\link{mean_to_d}},
\code{\link{mean_to_d2}},
\code{\link{t_to_d}},
\code{\link{f_to_d}},
\code{\link{p_to_d1}},
\code{\link{p_to_d2}}, 
\code{\link{ancova_to_d1}}, 
\code{\link{ancova_to_d2}}, 
\code{\link{tt.ancova_to_d}}, 
\code{\link{f.ancova_to_d}}, 
\code{\link{r_to_d}}, 
\code{\link{p.ancova_to_d1}},
\code{\link{p.ancova_to_d2}},
\code{\link{lor_to_d}},
\code{\link{prop_to_or}}, 
\code{\link{prop_to_d}}, 
\code{\link{r_from_chi}}, 
\code{\link{r_from_d}},
\code{\link{r_from_d1}},
\code{\link{r_from_t}}
}
\keyword{ arith }