% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/project.R
\name{project}
\alias{project}
\title{Project population dynamics}
\usage{
project(
  A,
  vector = "n",
  time = 100,
  standard.A = FALSE,
  standard.vec = FALSE,
  return.vec = TRUE,
  Aseq = "unif",
  Astart = NULL,
  draws = 1000,
  alpha.draws = "unif",
  PREcheck = TRUE
)
}
\arguments{
\item{A}{a matrix, or list of matrices. If \code{A} is a matrix, then 
\code{project} performs a 'deterministic' projection, where the matrix
does not change with each timestep. If \code{A} is a list of matrices, then 
\code{project} performs a 'stochastic' projection where the matrix varies 
with each timestep. The sequence of matrices is determined using \code{Aseq}. 
Matrices must be square, non-negative and numeric. If \code{A} is a list, 
all matrices must have the same dimension. 'Projection' objects inherit
names from \code{A}: if \code{A} is a matrix, stage names (in mat and 
vec slots) are inherited from its column names..}

\item{vector}{(optional) a numeric vector or matrix describing 
the age/stage distribution(s) used to calculate the projection. Single
population vectors can be given either as a numeric vector or 
one-column matrix. Multiple vectors are specified as a matrix, where 
each column describes a single population vector. Therefore the number
of rows of the matrix should be equal to the matrix dimension, whilst the 
number of columns gives the number of vectors to project. \code{vector} may
also take either "n" (default) to calculate the set of stage-biased projections 
(see details), or "diri" to project random population vectors drawn from a 
dirichlet distribution (see details).}

\item{time}{the number of projection intervals.}

\item{standard.A}{(optional) if \code{TRUE}, scales each matrix in \code{A}
by dividing all elements by the dominant eigenvalue. This standardises 
asymptotic dynamics: the dominant eigenvalue of the scaled matrix is 1. 
Useful for assessing transient dynamics.}

\item{standard.vec}{(optional) if \code{TRUE}, standardises each \code{vector} 
to sum to 1, by dividing each vector by its sum. Useful for assessing projection
relative to initial population size.}

\item{return.vec}{(optional) if \code{TRUE}, returns the time series of 
demographic (st)age vectors as well as overall population size.}

\item{Aseq}{(optional, for stochastic projections only) the sequence of 
matrices in a stochastic projection. 
\code{Aseq} may be either:
\itemize{
 \item "unif" (default), which results in every matrix in \code{A} having an 
 equal, random chance of being chosen at each timestep.
 \item a square, nonnegative left-stochastic matrix describing a first-order 
 Markov chain used to choose the matrices. The transitions are defined COLUMNWISE: 
 each column j describes the probability of choosing stage (row) i at time t+1, 
 given that stage (column) j was chosen at time t. \code{Aseq}  should have the 
 same dimension as the number of matrices in \code{A}. 
 \item a numeric vector giving a specific sequence which corresponds to the
 matrices in \code{A}.
 \item a character vector giving a specific sequence which corresponds to the
 names of the matrices in \code{A}.
}}

\item{Astart}{(optional) in a stochastic projection, the matrix with which to
initialise the projection (either numeric, corresponding to the matrices in 
\code{A}, or character, corresponding to the names of matrices in \code{A}). 
When \code{Astart = NULL} (the default), a random initial matrix is chosen.}

\item{draws}{if \code{vector="diri"}, the number of population vectors drawn
from dirichlet.}

\item{alpha.draws}{if \code{vector="diri"}, the alpha values passed to 
\code{rdirichlet}: used to bias draws towards or away from a certain population
structure.}

\item{PREcheck}{many functions in \code{popdemo} first check Primitivity, 
Reducibility and/or Ergodicity of matrices, with associated warnings and/or 
errors if a matrix breaks any assumptions. Set \code{PREcheck=FALSE} if you
want to bypass these checks.}
}
\value{
A \code{\link{Projection-class}} item. 
'Projection' objects inherit from a standard array, and can be treated as 
such. Therefore, if if \code{vector} is specified, the 'Projection' object will 
behave as: 
\itemize{
 \item if a single \code{vector} is given, a numeric vector of population sizes 
 of length \code{time+1}
 \item if multiple \code{vector}s are given, a numeric matrix of population 
 projections where each column represents a single population projection and 
 is of length \code{time+1}
 \item if \code{vector="n"}, a numeric matrix of population projections where each column 
 represents a single stage-biased projection and is of length \code{time+1}.
 \item if \code{vector="diri"}, a numeric matrix of population projections where each 
 column represents projection of a single vector draw and each column is of 
 length \code{time+1}
}

See documentation on \code{\link{Projection-class}} objects to understand how 
to access other slots (e.g. (st)age vectors through the population projection) 
and for S4 methods (e.g. plotting projections).
Some examples for understanding the structure of 3D arrays returned when 
\code{return.vec=TRUE}: when projecting a 3 by 3 matrix for >10 time intervals 
(see examples), element [11,3,2] represents the density of stage 3 at time 10 
for either vector 2 (multiple vectors), stage-bias 2 (\code{vector="n"}) or draw 2 
(\code{vector="diri"}); note that because element 1 represents t=0, then t=10 
is found at element 11. The vector [,3,2] represents the time series of densities 
of stage 3 in the projection of vector 2 / stage-bias 2 / draw 2. The matrix [,,2] 
represents the time series of all stages in the projection of vector 2 / stage-bias 
2 / draw 2.

Note that the projections inherit the labelling from \code{A} and \code{vector}, if
it exists. Both stage and vector names are taken from the COLUMN names of \code{A} 
and \code{vector} respectively. These may be useful for selecting from the
\code{projection} object, and for labelling graphs when plotting Projection
objects.
}
\description{
Project dynamics of a specified population matrix projection model.
}
\details{
If \code{vector} is specified, \code{project} will calculate population 
dynamics through time by projecting this vector / these vectors through 
\code{A}. If multiple vectors are specified, a separate population projection
is calculated for each.
 
If \code{vector="n"}, \code{project} will automatically project the set of 
'stage-biased' vectors of \code{A}. Effectively, each vector is a population
consisting of all individuals in one stage. These projections are achieved using a 
set of standard basis vectors equal in number to the dimension of \code{A}.
The vectors have every element equal to 0, except for a single element equal to 1,  
i.e. for a matrix of dimension 3, the set of stage-biased vectors are: 
\code{c(1,0,0)}, \code{c(0,1,0)} and \code{c(0,0,1)}. Stage-biased projections are 
useful for seeing how extreme transient dynamics can be.

If \code{vector="diri"}, \code{project} draws random population vectors from 
the dirichlet distribution. \code{draws} gives the number of population vectors
to draw. \code{alpha.draws} gives the parameters for the dirichlet and can be
used to bias the draws towards or away from certain population structures.
The default is \code{alpha.draws="unif"}, which passes \code{rep(1,dim)} (where
dim is the dimension of the matrix), resulting in an equal probability of 
any random population vector. Relative values in the vector give the population
structure to focus the distribution on, and the absolute value of the vector
entries (and their sum) gives the strength of the distribution: values greater
than 1 make it more likely to draw from nearby that population structure, 
whilst values less than 1 make it less likely to draw from nearby that population
structure.

Projections returned are of length \code{time+1}, as the first element 
represents the population at \code{t=0}.

Projections have their own plotting method (see \code{\link{Projection-plots}})
to enable easy graphing.

In addition to the examples below, see the "Deterministic population dynamics" 
and "Stochastic population dynamics" vignettes for worked examples that use 
the \code{project} function.
}
\examples{
  ### USING PROJECTION OBJECTS

  # Create a 3x3 PPM
  ( A <- matrix(c(0,1,2,0.5,0.1,0,0,0.6,0.6), byrow=TRUE, ncol=3) )

  # Project stage-biased dynamics of A over 70 intervals
  ( pr <- project(A, vector="n", time=70) )
  plot(pr)

  # Access other slots
  vec(pr)  #time sequence of population vectors
  bounds(pr)  #bounds on population dynamics
  mat(pr)  #matrix used to create projection
  Aseq(pr)  #sequence of matrices (more useful for stochastic projections)
  projtype(pr)  #type of projection
  vectype(pr)  #type of vector(s) initiating projection

  # Extra information on the projection
  nproj(pr)  #number of projections
  nmat(pr)  #number of matrices (more usefulk for stochastic projections)
  ntime(pr)  #number of time intervals
  
  # Select the projection of stage 2 bias
  pr[,2]

  # Project stage-biased dynamics of standardised A over 30 intervals
  ( pr2 <- project(A, vector="n", time=30, standard.A=TRUE) )
  plot(pr2)

  #Select the projection of stage 2 bias
  pr2[,2]

  # Select the density of stage 3 in bias 2 at time 10
  vec(pr2)[11,3,2]

  # Select the time series of densities of stage 2 in bias 1
  vec(pr2)[,2,1]

  #Select the matrix of population vectors for bias 2
  vec(pr2)[,,2]

  # Create an initial stage structure
  ( initial <- c(1,3,2) )

  # Project A over 50 intervals using a specified population structure
  ( pr3 <- project(A, vector=initial, time=50) )
  plot(pr3)

  # Project standardised dynamics of A over 10 intervals using 
  # standardised initial structure and return demographic vectors
  ( pr4 <- project(A, vector=initial, time=10, standard.vec=TRUE, 
                   standard.A=TRUE, return.vec=TRUE) )
  plot(pr4)

  # Select the time series for stage 1
  vec(pr4)[,1]

  ### DETERMINISTIC PROJECTIONS

  # Load the desert Tortoise matrix
  data(Tort)

  # Create an initial stage structure
  Tortvec1 <- c(8, 7, 6, 5, 4, 3, 2, 1)
  
  # Create a projection over 30 time intervals
  ( Tortp1 <- project(Tort, vector = Tortvec1, time = 10) )

  # plot p1
  plot(Tortp1)
  plot(Tortp1, bounds = TRUE) #with bounds
 
  # new display parameters
  plot(Tortp1, bounds = TRUE, col = "red", bty = "n", log = "y", 
       ylab = "Number of individuals (log scale)",
       bounds.args = list(lty = 2, lwd = 2) )

  # multiple vectors
  Tortvec2 <- cbind(Tortvec1, c(1, 2, 3, 4, 5, 6, 7, 8))
  plot(project(Tort, vector = Tortvec2), log = "y")
  plot(project(Tort, vector = Tortvec2), log = "y", labs = FALSE) #no labels

  # dirichlet distribution 
  # darker shading indicates more likely population size
  Tortshade <- project(Tort, time = 30, vector = "diri", standard.A = TRUE,
               draws = 500, alpha.draws = "unif")
  plot(Tortshade, plottype = "shady", bounds = TRUE)
  
  ### STOCHASTIC PROJECTIONS
  # load polar bear data
  data(Pbear)
  
  # project over 50 years with uniform matrix selection
  Pbearvec <- c(0.106, 0.068, 0.106, 0.461, 0.151, 0.108)
  p2 <- project(Pbear, Pbearvec, time = 50, Aseq = "unif")

  # stochastic projection information
  Aseq(p2)
  projtype(p2)
  nmat(p2)
  
  # plot
  plot(p2, log = "y")
  
}
\seealso{
\code{\link{Projection-class}} \code{\link{Projection-plots}}
}
\concept{projection project population}
