\name{robustSE}
\alias{robustSE}
\title{Robust standard error
}
\author{Mike Cheung with modifications by AC Del Re
}

\description{ When the correlation between dependent effect sizes are unknown, one approach is to conduct the meta-analysis by assuming that the effect sizes are independent. Then, Hedges et al. (2010) robust standard error procedure can be calculated to adjust for dependence. 
}
\usage{
robustSE(model, cluster=NULL, CI=.95, digits=3)
}
\arguments{
  \item{model}{omnibus or moderator model object fitted from mareg() function.  
}
  \item{cluster}{Name of variable where the dependencies are present. This will typically be the variable for \emph{study id} where \code{length(unique(study_id))>1}.
  }
 \item{CI}{Confidence interval. Defaults to .95.  
}
 \item{digits}{Number of digits to output. Defaults to 3.  
}
}
\value{
\item{estimate}{ Meta-regression coefficient estimate.
} 
\item{se}{Adjusted Standard error of the estimate coefficient.
}
\item{t}{ t-value.
}
\item{ci.l}{ Adjusted Lower 95\% confidence interval.
}
\item{ci.u}{ Adjusted Upper 95\% confidence interval.
}
\item{p}{ p-value.
}
}

\references{Hedges, L. V., Tipton, E., & Johnson, M. C. (2010). Robust variance estimation in meta-regression with dependent effect size estimates. Research Synthesis Methods, 1(1), 39-65. doi:10.1002/jrsm.5
 
 Cheung, M.W.L. (2012). metaSEM: An R package for meta-analysis using structural equation modeling. Manuscript submitted for publication.
}
\seealso{
\code{\link{mareg}},

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

# Adjusted SE
robustSE(mareg(g~ mod + mods2, var = var.g, method = "REML", data = df))

   
}



