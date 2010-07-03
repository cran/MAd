\name{MAd-package}
\alias{MAd-package}
\alias{MAd}
\docType{package}
\title{Meta-Analysis with Mean Differences
}
\description{ This package contains a variety of functions relevant to conducting a mean differences meta-analysis using recommended procedures as described in \emph{The Handbook of Research Synthesis and Meta-Analysis} (Cooper, Hedges, and Valentine, 2009). The goal in creating this package was to provide user-friendly functions to assist researchers in the process of conducting a meta-analysis, from the initial to final stages of their analytic endeavor. The meta-analyst can begin their project by using \code{MAd} functions to derive d (standardized mean differences) and g (unbiased d) from a variety of statistics/values reported in the primary studies (e.g., raw means and sd, t-test). Then, the analyst can aggregate all within-study effect sizes (while accounting for within-study correlations among outcome measures and eliminating any dependencies in the dataset) based on recommended procedures by Gleser & Olkin (1994; 2009), calculate omnibus effect sizes under a fixed and random effects model, and assess for significant moderators (categorical and continuous, single and multi-predictor models) in the dataset. Finally, the meta-analyst can use one of several user-friendly graphics functions to visually represent their data. 
}
\details{
\tabular{ll}{
Package: \tab MAd\cr
Type: \tab Package\cr
Version: \tab 0.6\cr
Date: \tab 2010-06-21\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
The \code{MAd} package has integrated functions to facilitate the meta-analytic process at nearly every analytical stage. There are five broad areas of analysis that the \code{MAd} package targets:
 
1. Computations to Calculate Mean Differences: 

There are a variety of functions to compute d (standardized mean difference) and g (unbiased d) from various designs reported in the primary studies. Most functions were derived from Borenstein's chapter in \emph{The Handbook of Research Synthesis and Meta-Analysis} (Cooper, Hedges, & Valentine, 2009; pp. 228-234).  For additional conversion formulas see the \bold{compute.es} package:  \url{http://CRAN.R-project.org/package=compute.es}


2. Within-Study Aggregation of Effect Sizes:

This fuction will simultaneously aggregate all within-study effect sizes while taking into account the correlations among the within-study outcome measures (Gleser & Olkin 2009; Gleser & Olkin 2009; Hedges & Olkin, 1985; Rosenthal et al., 2006). The default imputed correlation between within-study effect sizes is set at .50 (Wampold et al., 1997) and will compute an aggregated effect size for each study. This default of .50 is adjustable.  This \code{MAd} aggregation function implements Gleser & Olkin's (1994; 2009) recommended procedures for aggregating dependent effect sizes. To our knowledge, this is the first statistical package/program to explicitly utilize and automate this type of aggregation procedure, which has a dual effect of saving the researcher \bold{substantial} time while improving the accuracy of their analyses. 

3. Fixed and Random Effects Omnibus Analysis: 

This package contains all the relevant functions to calculate fixed and random effects omnibus effect sizes, outputting the omnibus (i.e., overall) effect size, variance, standard error, upper and lower confidence intervals, and the Q-statistic (heterogeneity test). 

4. Moderator Analyses:

There are user-friendly functions to compute fixed and random effects moderator analyses. These include single and multiple predictor models for both categorical and continuous moderator data. 

5. Graphics:

This package has a variety of functions visually representing data. This includes boxplots and meta-regression scatterplots.

6. Sample of Additional Functions:

Export MA output to nicely formatted Word tables.

}
\author{AC Del Re & William T. Hoyt

Maintainer: AC Del Re \email{acdelre@gmail.com}
}
\references{ Cooper, H., Hedges, L.V., & Valentine, J.C. (2009). \emph{The handbook of research synthesis and meta-analysis} (2nd edition). New York: Russell Sage Foundation.

Gleser & Olkin (1994). Stochastically dependent effect sizes. In H. Cooper, & L. V. Hedges, & J. C.(Eds.), \emph{The handbook of research synthesis} (pp. 339-356). New York: Russell Sage Foundation.
  
Gleser & Olkin (2009). Stochastically dependent effect sizes. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta analysis} (pp. 357-376). New York: Russell Sage Foundation.

 Viechtbauer, W (2010). metafor: Meta-Analysis Package for
  R. R package version 1.1-0. http://CRAN.R-project.org/package=metafor

Wampold, B. E., Mondin, G. W., Moody, M., Stich, F., Benson, K., & Ahn, H. (1997). A meta-analysis of outcome studies comparing bona fide psychotherapies: Empiricially, 'all must have prizes.' \emph{Psychological Bulletin, 122(3)}, 203-215.
}
\keyword{ package }
\examples{ 
# Examples for each broad area:

id <- c(1, 1:19)
n.1 <- c(10,20,13,22,28,12,12,36,19,12,36,75,33,121,37,14,40,16,14,20)
n.2 <- c(11,22,10,20,25,12,12,36,19,11,34,75,33,120,37,14,40,16,10,21)
g <- c(.68,.56,.23,.64,.49,-.04,1.49,1.33,.58,1.18,-.11,1.27,.26,.40,.49,
       .51,.40,.34,.42,1.16)
var.g <- c(.08,.06,.03,.04,.09,.04,.009,.033,.0058,.018,.011,.027,.026,.0040,
           .049,.0051,.040,.034,.0042,.016)
mod <- factor(c(rep(c(1,1,2,3),5)))
mods2 <- c(rep(1:5,4))
df <- data.frame(id, n.1, n.2, g, var.g, mod, mods2)
  

# 1. Computations to Mean Differences:
 
# For example, suppose the primary study reported an log odds ratio for different 
# proportions between 2 groups. Then, running:

lor_to_d(.9070,.0676)
           
# reported log odds ratio (lor = .9070) and variance (.0676) will output the 
# standardized mean difference (d) and its variance (var.d) for use in 
# the meta-analysis.

# 2. Within-Study Aggregation of Effect Sizes: 

agg(id = id, g = g, var = var.g, n.1 = n.1, n.2 = n.2, data=df) 

# where df = data.frame with columns for id, g (unbiased standardized
# mean difference), var.g (variance of g), n.1 (sample size of group
# one), and n.2 (sample size of comparison group) with multiple rows per
# study. Outputs an aggregated data.frame with 1 effect size per study. 

# 3. Fixed and Random Effects Omnibus Analysis

omni(g = g, var = var.g, data = df, type="weighted", method = "random")

# where df = data.frame with columns for id, g (unbiased standardized
# mean difference), var.g (variance of g), n.1 (sample size of group
# one), and n.2 (sample size of comparison group).
 
# 4. Moderator Analyses:

# Random effects
mareg(g~ mod + mods2, var = var.g, method = "REML", data = df) 
 
# where df = data.frame with columns  g (unbiased standardized
# mean difference), var (variance of g), and moderators.

# 5. Graphics:

\dontrun{plotcon(g = g, var = var.g, mod = mod, data = df, method= "random", 
modname= "Moderator") }

# Additional Functions

# Export MA output to nicely formatted Word tables.

# install R2wd
# install.packages('R2wd', dependencies = TRUE)

# mareg fuction

temp <- mareg(g~ mod + mods2, var = var.g, method = "REML", data = df)

# Export data to Word in formatted table

#  wd(temp, get = TRUE, new = TRUE)
}
\seealso{
\bold{RcmdrPlugin.MAd} package:  \url{http://CRAN.R-project.org/package=RcmdrPlugin.MAd};

\bold{MAc} package:  \url{http://CRAN.R-project.org/package=MAc};

\bold{RcmdrPlugin.MAc} package:  \url{http://CRAN.R-project.org/package=RcmdrPlugin.MAc};

\bold{compute.es} package:  \url{http://CRAN.R-project.org/package=compute.es};

\bold{metafor} package:  \url{http://CRAN.R-project.org/package=metafor}
}