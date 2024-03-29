\name{MAd-package}
\alias{MAd-package}
\alias{MAd}
\docType{package}
\title{Meta-Analysis with Mean Differences
}
\description{ The Mad package contains a variety of functions for conducting a mean differences meta-analysis using recommended procedures as described in \emph{The Handbook of Research Synthesis and Meta-Analysis} (Cooper, Hedges, and Valentine, 2009). The goal in creating this package was to provide user-friendly functions to assist researchers in the process of conducting a meta-analysis, from the initial to final stages of their analytic endeavor. The meta-analyst can begin their project by using \code{MAd} functions to derive d (standardized mean differences) and g (unbiased d) from a variety of statistics/values reported in the primary studies (e.g., raw means and sd, t-test). Then, the analyst can aggregate all within-study effect sizes (while accounting for within-study correlations among outcome measures and eliminating any dependencies in the dataset), calculate omnibus effect sizes under a fixed and random effects model, and assess for significant moderators (categorical and continuous, single and multi-predictor models) in the dataset. Finally, the meta-analyst can use one of several user-friendly graphics functions to visually represent their data. 
}
\details{
\tabular{ll}{
Package: \tab MAd\cr
Type: \tab Package\cr
Version: \tab 0.8-2\cr
Date: \tab 2014-12-23\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
The \code{MAd} package has integrated functions to facilitate the meta-analytic process at nearly every analytical stage. There are five broad areas of analysis that the \code{MAd} package targets:
 
1. Computations to Calculate Mean Differences: 

There are a variety of functions to compute d (standardized mean difference) and g (unbiased d) from various designs reported in the primary studies. Most functions were derived from Borenstein's chapter in \emph{The Handbook of Research Synthesis and Meta-Analysis} (Cooper, Hedges, & Valentine, 2009; pp. 228-234).  For additional conversion formulas see the \bold{compute.es} package:  \url{https://CRAN.R-project.org/package=compute.es}


2. Within-Study Aggregation of Effect Sizes:

This fuction will simultaneously aggregate all within-study effect sizes while taking into account the correlations among the within-study outcome measures (Gleser & Olkin 2009; Gleser & Olkin 2009; Hedges & Olkin, 1985; Rosenthal et al., 2006). The default imputed correlation between within-study effect sizes is set at .50 (Wampold et al., 1997) and will compute an aggregated effect size for each study. This default of .50 is adjustable and can vary between outcome types.  This \code{MAd} aggregation function implements Gleser & Olkin (1994; 2009) and Borenstein et al. (2009) procedures for aggregating dependent effect sizes. To our knowledge, this is the first statistical package/program to explicitly utilize and automate this type of aggregation procedure, which has a dual effect of saving the researcher \bold{substantial} time while improving the accuracy of their analyses. 

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
\references{ 
Borenstein (2009). Effect sizes for continuous data. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta analysis} (pp. 279-293). New York: Russell Sage Foundation.

Cooper, H., Hedges, L.V., & Valentine, J.C. (2009). \emph{The handbook of research synthesis and meta-analysis} (2nd edition). New York: Russell Sage Foundation.

Gleser & Olkin (1994). Stochastically dependent effect sizes. In H. Cooper, & L. V. Hedges, & J. C.(Eds.), \emph{The handbook of research synthesis} (pp. 339-356). New York: Russell Sage Foundation.
  
Gleser & Olkin (2009). Stochastically dependent effect sizes. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta analysis} (pp. 357-376). New York: Russell Sage Foundation.

 Viechtbauer, W (2010). metafor: Meta-Analysis Package for
  R. R package version 1.1-0. https://CRAN.R-project.org/package=metafor

Wampold, B. E., Mondin, G. W., Moody, M., Stich, F., Benson, K., & Ahn, H. (1997). A meta-analysis of outcome studies comparing bona fide psychotherapies: Empiricially, 'all must have prizes.' \emph{Psychological Bulletin, 122(3)}, 203-215.
}
\keyword{ package }
\examples{ 
## EXAMPLES FOR EACH BROAD AREA

# SAMPLE DATA:
MA <- data.frame(id=factor(rep(1:5, 3)),
                 measure=c(rep("dep",5), rep("anx",5), rep("shy",5)),
                 es=c(rnorm(5, 0.8, .2), rnorm(5, 0.5, .1), rnorm(5, 0.4, .1)),
                 var.es=abs(rnorm(5*3,0.05, .03)),
                 nT=round(rnorm(5*3, 30, 5),0),
                 nC=round(rnorm(5*3, 30, 5),0),
                 mod1=factor(rep(c("a","b","c","d","e"),3)),
                 mod2=rep(seq(20, 60, 10), 3))
                 
# 1. COMPUTE MEAN DIFFERENCE STATISTIC FROM 
# REPORTED STATS (GENERALLY FROM A PRIMARY STUDY):
 
# suppose the primary study reported an log odds ratio for different 
# proportions between 2 groups. Then, running:

lor_to_d(.9070,.0676)
           
# reported log odds ratio (lor = .9070) and variance (.0676) will output the 
# standardized mean difference (d) and its variance (var.d) that can be used in 
# a meta-analysis.

## 2. ACCOUNT FOR DEPENDENCIES: WITHIN-STUDY EFFECT SIZES (ES): 

## 2 EXAMPLES:  
# EXAMPLE 1:  AGGREGATING EFFECT SIZES FOR A DATA FRAME 
# (MULTIPLE STUDIES AT LEAST SOME OF WHICH HAVE MULTIPLE DEPENDENT EFFECT SIZES)   
# EXAMPLE 2:  AGGREGATING EFFECT SIZES FOR SINGLE STUDY WITH THREE OR MORE
# EFFECT SIZES WHEN PAIRS OF DVS HAVE DIFFERENT CORRELATIONS 

## EXAMPLE 1:  MA IS A DATA FRAME CONTAINING MULTIPLE STUDIES (id),
## EACH WITH MULTIPLE EFFECT SIZES (CORRELATIONS BETWEEN ALL PAIRS OF DVS ARE r=.5.)


# AGGREGATION PROCEDURE:  
# method="GO1"; GLESER AND OLKIN (1994) PROCEDURE WHEN d IS COMPUTED 
# USING POOLED SD IN THE DENOMINATOR

MA1 <- agg(id=id, es=es, var=var.es, n.1=nT, n.2=nC, cor = .5, method="GO1", data=MA)   

MA1

## EXAMPLE 2: STUDY 1 COMPARES A TREATMENT AND CONTROL GROUP ON 
## THREE OUTCOME MEASURES (DEPRESSION, ANXIETY, and SHYNESS).
# THE CORRELATION AMONG THE THREE PAIRS OF DVS ARE r12=.5, r13=.2, and r23=.3.  

study1 <- data.frame( id=factor(rep(1, 3)), 
                      measure=c("dep", "anx", "shy"),
                      es=c(0.8, 0.5, 0.4), 
                      var.es=c(0.01, 0.02, 0.1), 
                      nT=rep(30, 3), 
                      nC=rep(30, 3))  

# ONE WOULD CONSTRUCT THE CORRELATION MATRIX AS FOLLOWS:  

cors <- matrix(c(1,.5,.2, 
                 .5,1,.3, 
                 .2,.3,1), nrow=3)  

# AGGREGATION PROCEDURE:
# method="GO1"; GLESER AND OLKIN (1994) PROCEDURE WHEN d 
# IS COMPUTED USING POOLED SD IN THE DENOMINATOR 

agg(id=id, es=es, var=var.es, n.1=nT, n.2=nC, cor=cors, method="GO1", mod = NULL, data=study1)   


# where MA = data.frame with columns for id, es (standardized
# mean difference), var.es (variance of es), n.1 (sample size of group
# one), and n.2 (sample size of comparison group) with multiple rows per
# study. Outputs an aggregated data.frame with 1 effect size per study. 

## 3.OMNIBUS ANALYSIS

# FIRST ADD MODERATORS TO THE AGGREGATED DATASET:

MODS <- data.frame(id=1:5,
                   mod1=factor(c("a","b","a","b","b")),
                   mod2=as.numeric(c(20, 30, 25, 35, 40)))
                   
MA2 <- merge(MA1, MODS, by='id')

# Random Effects
m0 <- mareg(es1~ 1, var = var1, method = "REML", data = MA2)

# where MA = data.frame with columns for id, es (standardized
# mean difference), var.es (variance of es), n.1 (sample size of group
# one), and n.2 (sample size of comparison group).
 
# view output:
summary(m0)


# 4. MODERATOR ANALYSIS:

# Random Effects
m1 <- mareg(es1~ mod1 + mod2 , var = var1, method = "REML", data = MA2) 
 
# view output:
summary(m1)

# 5. Graphics:

\dontrun{plotcon(g = es1, var = var1, mod = mod1, data = MA2, method= "random", 
modname= "Moderator") }

# Additional Functions

# Export MA output to nicely formatted Word tables.

# install R2wd
# install.packages('R2wd', dependencies = TRUE)

# Export data to Word in formatted table

#  wd(m1, get = TRUE, new = TRUE)
}
\seealso{

\bold{compute.es} package:  \url{https://CRAN.R-project.org/package=compute.es};

\bold{metafor} package:  \url{https://CRAN.R-project.org/package=metafor}
}