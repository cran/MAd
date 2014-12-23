\name{agg}
\alias{agg}
\title{Aggregate Dependent Effect Sizes
}
\description{This function will simultaneously aggregate all within-study effect sizes while taking into account the correlations among the within-study outcomes. The default correlation between outcome measures is set at .50 (Wampold et al., 1997), and can be adjusted as needed.  An aggregate effect size and its variance is computed for each study. This \code{MAd} aggregation function implements Gleser & Olkin (1994) and Borenstein et al. (2009) procedures for aggregating dependent effect sizes..}
\usage{
agg(id, es, var, n.1=NULL, n.2=NULL, method = "BHHR", cor = .50,  mod=NULL, data)
}
\arguments{
 \item{id}{Study id with multiple rows of same id.
}
\item{es}{Effect size. Use Cohen's d for GO1 and GO2 method and use Hedges g for BHHR method.
}
\item{var}{Variance of g.
}
\item{n.1}{Sample size of group one. Only required if using method='GO1' or 'GO2'.
}
\item{n.2}{Sample size of group two. Only required if using method='GO1' or 'GO2'.
}
\item{method}{'BHHR'= Borenstein, et al. (2009) procedure.'GO1'= Gleser and Olkin (1994) procedure using pooled standard deviation (SD). 'GO2'= Gleser and Olkin (1994) procedure using group 2 SD  (typically control group SD).  Default is 'BHHR'.
}
\item{cor}{Estimated correlation among within-study outcome variables. Input should either be a fixed correlation (default is r=.5) or a correlation matrix to allow different correlations between outcome types (in which case aggregation should typically be on one study at a time rather than the entire dataset--see examples below). 
}
\item{mod}{To aggregate by id and one moderator. If there are multiple levels of a categorical moderator within study and one can in derive seperate effect size estimates for each level within and between studies.  Default is NULL.
}
\item{data}{\code{data.frame} with above values.
}
}
\value{Outputs a \code{data.frame} with aggregated effect sizes and variances of effect sizes where each study is reduced to one row per study (unless aggregated by a moderator) by a weighted average formula.  
}
\references{ 
Borenstein (2009). Effect sizes for continuous data. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta analysis} (pp. 279-293). New York: Russell Sage Foundation.

Cooper, H., Hedges, L.V., & Valentine, J.C. (2009). \emph{The handbook of research synthesis and meta-analysis} (2nd edition). New York: Russell Sage Foundation. 

Gleser & Olkin (1994). Stochastically dependent effect sizes. In H. Cooper, & L. V. Hedges, & J. C.(Eds.), \emph{The handbook of research synthesis}(pp. 339-356). New York: Russell Sage Foundation.

Gleser & Olkin (2009). Stochastically dependent effect sizes. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta analysis} (pp. 357-377). New York: Russell Sage Foundation.

Shadish & Haddock (2009). Analyzing effect sizes: Fixed-effects models. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta analysis} (pp. 257-278). New York: Russell Sage Foundation. 
}
\author{AC Del Re & William T. Hoyt

Maintainer: AC Del Re \email{acdelre@gmail.com}
}
\examples{
## 2 EXAMPLES:  

# EXAMPLE 1:  Aggregating effect sizes for a data frame
# (multiple studies at least some of which have multiple
# effect sizes), assuming equal correlations (r=.5) between
# pairs of DVs.
# EXAMPLE 2:  Aggregating effect sizes for a single study
# with 3 or more effect sizes when pairs of DVs have 
# different correlations.

# LOAD DATA (EXAMPLE DATA FROM HOYT & DEL RE, 2015 SIMULATION):
data(dat.hoyt)

## EXAMPLE 1:  dat.hoyt is a data frame with multiple studies identified
## by variable 'id'.  Each study has multiple effect sizes based on 
## multiple DVs.  Correlations between all pairs of DVs are r=.5.

# NOTE: Based on a simulation study by Hoyt & Del Re (2015), it is
# recommended that methods "G01" and "G02" (Gleser and Olkin) 
# should aggregate Cohen's d, without using Hedges & Olkin's 
# recommended bias correction.  (Studies providing only a single
# effect size should still be corrected for bias, after aggregation.)

# Method "BHHR" should aggregate Hedges' g, after bias correction.


# Option 1: method="BHHR"; Borenstein et al. (2009) procedure. 
# Use with Hedges' g; can also be used with any other effect 
# size (e.g., z', LOR).

agg(id=id, es=g, var=vg,  cor=.5,
    method="BHHR", mod=NULL, data=dat.hoyt) 

#  Option 2: method="GO1"; Gleser & Olkin (1994) procedure when
#  d is computed using pooled sd in denominator. 

agg(id=id, es=d, var=vd, n.1=n.T, n.2=n.C, cor = .5, 
    method="GO1", mod = NULL, data=dat.hoyt)   

# Option 3: method="GO2"; Gleser & Olkin (1994) procedure when
# d is computed using sd.2 (typically control group sd)
# in denominator

agg(id=id, es=d, var=vd, n.1=n.T, n.2=n.C, cor = .5, 
    method="GO2", mod = NULL, data=dat.hoyt)   


## EXAMPLE 2: Single study comparing T and C group
## on three DVs:  depression, anxiety, and shyness
## r12=.5; r13=.2; r23=.3

data <- dat.hoyt[20:22,]

# Step 1:  Create the correlation matrix, based on r12, r13, and r23:  

cors <- matrix(c(1,.5,.2, 
                 .5,1,.3, 
                 .2,.3,1), nrow=3)  

# Step 2:  Aggregate using agg() function.

# Option 1: method="BHHR"; Borenstein et al. (2009) procedure. 
# Use with Hedges' g; can also be used with any other effect 
# size (e.g., z', LOR).

agg(id=id, es=g, var=vg,  cor=cors,
    method="BHHR", mod=NULL, data=data) 

#  Option 2: method="GO1"; Gleser & Olkin (1994) procedure when
#  d is computed using pooled sd in denominator. 

agg(id=id, es=d, var=vd, n.1=n.T, n.2=n.C, cor = cors, 
    method="GO1", mod = NULL, data=data)   

# Option 3: method="GO2"; Gleser & Olkin (1994) procedure when
# d is computed using sd.2 (typically control group sd)
# in denominator

agg(id=id, es=d, var=vd, n.1=n.T, n.2=n.C, cor = cors, 
    method="GO2", mod = NULL, data=data)   


## Citation ##
# Hoyt, W. T., & Del Re, A. C. (2013).  Comparison of methods for 
# aggregating dependent effect sizes in meta-analysis.  
# Manuscript submitted for publication.

  


}
\keyword{ aggregation }

