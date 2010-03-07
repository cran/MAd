\name{ComplData}
\alias{ComplData}
\title{Complete Dataset
}
\description{ Function to reduce data set with complete data for 1 predictor. Will accept a \code{data.frame} with multiple rows per study and reduce the \code{data.frame} to one effect size per study based on Hedges & Olkin's recommended aggregation procedures (1994; 2009).
}
\usage{
ComplData(meta, mod, type = "independent", cor = .50)
}
\arguments{
\item{meta}{\code{data.frame} with id, g (unbiased standardized mean difference), var.g (variance of g), n.1 (sample size of treatment group), and n.2 (sample size for comparison group)for each study. The function requires that the column name for these values be exactly as stated here.
}
  \item{mod}{ Moderator variable wanting to be retained for further analysis. 
}
  \item{type}{ Type of data aggregation. Options are \code{independent} or \code{dependent}. Choose \code{independent} if a fully independent data set is desired, for which the function will randomly select one factor level of the within-study moderator (if there are > 1 coded from the primary studies, e.g., the primary study reports early, mid, and late symptom change from CBT treatment). Choose \code{dependent} if all factor levels of the within-study moderator are to be kept in the data set. Default is \code{independent}.
}
\item{cor}{Estimated correlation among within-study outcome variables. Default is .50 based on Wampold et al. recommended procedures (1997).
}
}
\value{ Reduced \code{data.frame} (with complete data) for the moderator entered into the function. 
}
\author{AC Del Re & William T. Hoyt

Maintainer: AC Del Re \email{acdelre@gmail.com}
}
\references{
Gleser & Olkin (1994). Stochastically dependent effect sizes. In H. Cooper, & L. V. Hedges, & J. C.(Eds.), \emph{The handbook of research synthesis} (pp. 339-356). New York: Russell Sage Foundation.

Gleser & Olkin (2009). Stochastically dependent effect sizes. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta analysis} (pp. 357-376). New York: Russell Sage Foundation.

Raudenbush (2009). Analyzing effect sizes: Random-effects models. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta analysis} (pp. 295-316). New York: Russell Sage Foundation.  

Wampold, B. E., Mondin, G. W., Moody, M., Stich, F., Benson, K., & Ahn, H. (1997). A meta-analysis of outcome studies comparing bona fide psychotherapies: Empiricially, 'all must have prizes.' \emph{Psychological Bulletin, 122(3)}, 203-215.
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

ComplData(df, df$mod)
}
\keyword{data}

