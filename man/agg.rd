\name{agg}
\alias{agg}
\title{Meta-Analysis Aggregation
}
\description{This fuction has automated (i.e., will compute for all studies simultaneously) the process of aggregating within-study effect sizes while taking into account the correlations among the within-study outcome measures (Gleser & Olkin 2009; Gleser & Olkin 2009; Hedges & Olkin, 1985; Rosenthal et al., 2006). These functions default the correlation between within-study effect sizes at .50 (Wampold et al., 1997) and will compute the correct aggregated effect size for all studies. This default of .50 is adjustable.  \code{MAd} aggregation functions implement Gleser & Olkin (1994; 2009) recommendations for aggregating dependent correlations.}
\usage{
agg(id, g, var, n.1, n.2, cor = .50, mod=NULL, data)
}
\arguments{
 \item{id}{Study id with multiple rows of same id.
}
\item{g}{Hedges g (unbiased estimate of d) effect size.
}
\item{var}{Vaiance of g.
}
\item{n.1}{Sample size of group one.
}
\item{n.2}{Sample size of group two.
}
\item{cor}{Estimated correlation among within-study outcome variables. Default is .50 based on Wampold et al. recommended procedures (1997).
}
\item{mod}{Default is NULL. To aggregate by id and one moderator. If there are multiple levels of a categorical moderator within study and one can in derive seperate effect size estimates for each level within and between studies. However, there will be dependency issues and one way to resolve is shown below in the examples. 
}
\item{data}{\code{data.frame} with above values.
}
}
\value{Outputs a \code{data.frame} with aggregated effect sizes where each study is reduced to one row per study (unless aggregated by a moderator) by a weighted average formula. This formula is based on Gleser & Olkin's (1994) approach to aggregation of dependent effect sizes (see chapter 22, pp. 339-356). 
}
\references{ Cooper, H., Hedges, L.V., & Valentine, J.C. (2009). \emph{The handbook of research synthesis and meta-analysis} (2nd edition). New York: Russell Sage Foundation. 

Gleser & Olkin (1994). Stochastically dependent effect sizes. In H. Cooper, & L. V. Hedges, & J. C.(Eds.), \emph{The handbook of research synthesis}(pp. 339-356). New York: Russell Sage Foundation.

Gleser & Olkin (2009). Stochastically dependent effect sizes. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta analysis} (pp. 357-377). New York: Russell Sage Foundation.

Shadish & Haddock (2009). Analyzing effect sizes: Fixed-effects models. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta analysis} (pp. 257-278). New York: Russell Sage Foundation. 
}
\author{AC Del Re & William T. Hoyt

Maintainer: AC Del Re \email{acdelre@gmail.com}
}
\examples{
id <- c(1, 1:19)
n.1 <- c(10,20,13,22,28,12,12,36,19,12,36,75,33,121,37,14,40,16,14,20)
n.2 <- c(11,22,10,20,25,12,12,36,19,11,34,75,33,120,37,14,40,16,10,21)
g <- c(.68,.56,.23,.64,.49,-.04,1.49,1.33,.58,1.18,-.11,1.27,.26,.40,.49,
       .51,.40,.34,.42,1.16)
var.g <- c(.08,.06,.03,.04,.09,.04,.009,.033,.0058,.018,.011,.027,.026,.0040,
.049,.0051,.040,.034,.0042,.016)
tmt <- factor(c(rep(c(1,2,2,3),5)))
df <- data.frame(id, n.1,n.2, g, var.g, tmt)

# Examples

# aggregate to 1 id per study (independent sample)
agg(id = id, g = g, var = var.g, n.1 = n.1, n.2 = n.2, data=df)

# aggregate by id & a moderator (non-independent sample)
temp <- agg(id = id, g = g, var = var.g, n.1 = n.1, n.2 = n.2, mod = tmt, data=df)  

temp

# This function below will randomly select one within
# study level of the moderator (if there are more than one) and output an
# independent sample. Replace temp with the name of your data.
do.call(rbind, lapply(split(temp, temp$id), 
          function(.data) .data[sample(nrow(.data), 1),]))


}
\keyword{ data }

