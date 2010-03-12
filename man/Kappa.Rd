\name{Kappa}
\alias{Kappa}
\title{Inter-Rater Agreement
}
\description{Kappa coefficients for assessing inter-rater agreement between two coders (for categorical \cr
variables/moderators).
}
\usage{
Kappa(rater1, rater2)
}
\arguments{
  \item{rater1}{First rater of categorical variable to be analyzed. This vector needs to be in a separate column if in a \code{data.frame}.
}
  \item{rater2}{Second rater on same categorical variable to be analyzed. This vector needs to be in a separate column if in a \code{data.frame}.
}
}
\value{Kappa coefficients for inter-rater reliability (categorical variables).
}
\author{AC Del Re & William T. Hoyt

Maintainer: AC Del Re \email{acdelre@gmail.com}
}
\keyword{ models }
