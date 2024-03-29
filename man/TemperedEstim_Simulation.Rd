% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MonteCarloFunction.R
\name{TemperedEstim_Simulation}
\alias{TemperedEstim_Simulation}
\title{Monte Carlo Simulation}
\usage{
TemperedEstim_Simulation(
  ParameterMatrix,
  SampleSizes = c(200, 1600),
  MCparam = 100,
  TemperedType = c("CTS", "TSS", "NTS", "MTS", "GTS", "KRTS", "RDTS"),
  Estimfct = c("ML", "GMM", "Cgmm", "GMC"),
  HandleError = TRUE,
  saveOutput = FALSE,
  SeedOptions = NULL,
  eps = 1e-06,
  algo = NULL,
  regularization = NULL,
  WeightingMatrix = NULL,
  t_scheme = NULL,
  alphaReg = NULL,
  t_free = NULL,
  nb_t = NULL,
  subdivisions = NULL,
  IntegrationMethod = NULL,
  randomIntegrationLaw = NULL,
  s_min = NULL,
  s_max = NULL,
  ncond = NULL,
  IterationControl = NULL,
  methodR = "TM",
  ...
)
}
\arguments{
\item{ParameterMatrix}{The matrix is to be composed of vectors, row by row.
Each vector must fit the pattern of theta of the \code{TemperedType}.}

\item{SampleSizes}{Sample sizes to be used to simulate the data. By default,
we use 200 (small sample size) and 1600 (large sample size);
vector of integer.}

\item{MCparam}{Number of Monte Carlo simulation for each couple of parameter,
default=100; integer}

\item{TemperedType}{A String. Either "CTS", "TSS", "NTS", "MTS", "GTS",
"KRTS", "RDTS".}

\item{Estimfct}{The estimation function to be used. A String.
Either "ML", "GMM", "Cgmm", or "GMC".}

\item{HandleError}{Logical flag: if set to TRUE, the simulation doesn't stop
when an error in the estimation function is encountered. A vector of
(size 4) NA is saved and the the simulation carries on. See details.}

\item{saveOutput}{Logical flag: if set to TRUE, a csv file (for each couple
of parameter) with the the estimation
information is saved in the current directory. See details.}

\item{SeedOptions}{List to control the seed generation. See details.}

\item{eps}{Numerical error tolerance. \code{1e-06} by default.}

\item{algo}{algorithm: For GMM: \code{"2SGMM"} is the two step GMM proposed
by Hansen (1982). \code{"CueGMM"} and \code{"ITGMM"} are respectively the
continuous updated and the iterative GMM proposed by Hansen, Eaton et Yaron
(1996) and adapted to the continuum case. For GMC: \code{"2SGMC", "CueGMC"}.
For Cgmm: \code{"2SCgmm", "CueCgmm", ...}.}

\item{regularization}{regularization scheme to be used for moment methods,
one of \code{"Tikhonov"} (Tikhonov), \code{"LF"} (Landweber-Fridmann) and
\code{"cut-off"} (spectral cut-off).}

\item{WeightingMatrix}{type of weighting matrix used to compute the
objective function for the GMM and GMC methods, one of \code{"OptAsym"} (the
optimal asymptotic), \code{"DataVar"} (the data driven, only for GMM) and
\code{"Id"} (the identity matrix).}

\item{t_scheme}{scheme used to select the points for the GMM method where the
moment conditions are evaluated, one of \code{"equally"} (equally placed),
\code{"NonOptAr"} (non optimal arithmetic placement), \code{"uniformOpt"}
(uniform optimal placement), \code{"ArithOpt"} (arithmetic optimal placement)
, \code{"Var Opt"} (optimal variance placement) and \code{"free"} (users need
to pass their own set of points in ...).}

\item{alphaReg}{value of the regularisation parameter; numeric. Example Value
could be ==0.01.}

\item{t_free}{sequence, if \code{t_scheme=="free"}.}

\item{nb_t}{integer, if you set \code{t_scheme <- "equally"}. nb_t could be
== 20 for example.}

\item{subdivisions}{Number of subdivisions used to compute the different
integrals involved in the computation of the objective function for the Cgmm
method (to minimise); numeric.}

\item{IntegrationMethod}{Numerical integration method to be used to
approximate the (vectorial) integrals for the Cgmm method. Users can choose
between "Uniform" discretization or the "Simpson"'s rule (the 3-point
Newton-Cotes quadrature rule).}

\item{randomIntegrationLaw}{Probability measure associated to the Hilbert
space spanned by the moment conditions for the Cgmm method.}

\item{s_min, s_max}{Lower and Upper bounds of the interval where the moment
conditions are considered for the Cgmm method; numeric.}

\item{ncond}{Integer. Number of moment conditions (until order \code{ncond})
for the GMC method. Must not be less than 3 for TSS, 6 for CTS, 5 for NTS.}

\item{IterationControl}{only used if algo = "IT..." or algo = "Cue..."
to control the iterations. See Details.}

\item{methodR}{A string. Method generates random variates of TS distribution.
"TM" by default. Switches automatically if the method is not applicable in
this way.}

\item{...}{Other arguments to be passed to the estimation function.}
}
\value{
If \code{saveOutput == FALSE}, the return object is a list of 2.
Results of the simulation are listed in \code{$outputMat}. If \code{
saveOutput == TRUE}, only a csv file is saved and nothing is returned.
}
\description{
Runs Monte Carlo simulation for a selected estimation method. The function
can save results in a file.
}
\details{
\strong{TemperedTyp} With the parameter 'TemperedTyp' you can choose the
tempered stable distribution you want to use. Here is a list of distribution
you can choose from:
\describe{
\item{TSS}{Tempered stabel subordinator: See \code{\link[=charTSS]{charTSS()}} for details.}
\item{CTS}{Classical tempered stable distribution: See \code{\link[=charCTS]{charCTS()}} for
details.}
\item{GTS}{Generalized classical tempered stable distribution: See
\code{\link[=charGTS]{charGTS()}} for details.}
\item{NTS}{Normal tempered stable distribution: See \code{\link[=charNTS]{charNTS()}} for
details.}
\item{MTS}{Modified tempered stable distribution: See \code{\link[=charMTS]{charMTS()}} for
details.}
\item{RDTS}{Rapid decreasing tempered stable distribution: See \code{\link[=charRDTS]{charRDTS()}}
for details.}
\item{KRTS}{Kim-Rachev tempered stable distribution: See \code{\link[=charKRTS]{charKRTS()}} for
details.}
}

\strong{Error Handling} It is advisable to set it to TRUE when user is
planning to launch long simulations as it will prevent the procedure to stop
if an error occurs for one sample data. The estimation function will produce
a vector of NA as estimated parameters related to this (error generating)
sample data and move on to the next Monte Carlo step.

\strong{Output file} Setting \code{saveOutput} to \code{TRUE} will have the
side effect of saving a csv file in the working directory. This file will
have \code{MCparam*length(SampleSizes)} lines and its columns will be:
\describe{
\item{alphaT, ...:}{the true value of the parameters.}
\item{data size:}{the sample size used to generate the simulated data.}
\item{seed:}{the seed value used to generate the simulated data.}
\item{alphaE, ...:}{the estimate of the parameters.}
\item{failure:}{binary: 0 for success, 1 for failure.}
\item{time:}{estimation running time in seconds.}
}
The file name is informative to let the user identify the value of the true
parameters, the MC parameters as well as the options selected for the
estimation method. The csv file is updated after each MC estimation which is
useful when the simulation stops before it finishes.

\strong{SeedOptions} If users does not want to control the seed generation,
they could ignore this argument (default value NULL). This argument can be
more useful when they wants to cut the simulation (even for one parameter
value) into pieces. In that case, they can control which part of the seed
vector they want to use.
\describe{
\item{MCtot:}{total values of MC simulations in the entire process.}
\item{seedStart:}{starting index in the seed vector. The vector extracted
will be of size MCparam.}
}

\strong{Estimfct} Additional parameters are needed for different estimation
functions. These are listed below for each function. The list of additional
parameters starts after the parameter \code{eps} in the parameter list.
\describe{
\item{For ML:}{ See usage of Maximum likelihood estimation in Kim et al.
(2008).No additional parameters are needed.}
\item{For GMM:}{Generalized Method of Moments by Feuerverger (1981).
The parameters \code{algo, alphaReg, regularization, WeightingMatrix, and
  t_scheme} must be specified.

Parameter \code{t_scheme}: One of the most important features of this
method is that it allows the user to choose how to place the points where
the moment conditions are evaluated. One can choose among 6 different
options. Depending on the option, further parameters have to be passed.
\describe{
\item{"equally":}{equally placed points in \code{min_t,max_t}. When
provided, user's \code{min_t} and \code{max_t} will be used (when
\code{Coinstrained == FALSE}).
}
\item{"NonOptAr":}{non optimal arithmetic placement.
}
\item{"uniformOpt":}{uniform optimal placement.
}
\item{"ArithOpt":}{arithmetic optimal placement.
}
\item{"Var Opt":}{optimal variance placement as explained above.
}
\item{"free":}{user needs to pass own set of points in \code{t_free}.
}
}

Parameter \code{WeightingMatrix}: One can choose among 3 different options:
\describe{
\item{"OptAsym":}{the optimal asymptotic choice.
}
\item{"DataVar":}{the covariance matrix of the data provided.
}
\item{"Id":}{the identity matrix.
}
}
}
\item{For Cgmm:}{Continuum Generalized Methods of Moments by Carrasco &
Kotchoni (2017). The parameters \code{algo, alphaReg, subdivisions,
  IntegrationMethod, randomIntegrationLaw, s_min, and s_max} must be
specified.
}
\item{For GMC:}{Generalized Method of Cumulants (GMC) by Massing, T.
(2022). The parameters \code{algo, alphaReg, regularization,
   WeightingMatrix, and ncond} must be specified.
}
}

\strong{IterationControl} If \code{algo = "IT..."} or \code{algo =
"Cue..."} the user can control each iteration by setting up the list
IterationControl which contains the following elements:
\describe{
\item{NbIter}{maximum number of iteration. The loop stops when NBIter is
reached; default = 10.}
\item{PrintIterlogical}{if set to TRUE, the value of the current parameter
estimation is printed to the screen at each iteration; default = TRUE.}
\item{RelativeErrMax}{the loop stops if the relative error between two
consecutive estimation steps is smaller than RelativeErrMax;
default = 1e-3.}
}

\strong{methodR} Random numbers must be generated for each MC study. For
each distribution, different methods are available for this (partly also
depending on alpha). For more information, the documentation of the
respective \code{r...()} distribution can be called up. By default, the fastest
method is selected. Since the deviation error can amplify to the edges of
alpha depending on the method, it is recommended to check the generated
random numbers once for each distribution using the density function before
starting the simulation.

\strong{Parallelization}  Parallelization of the function is possible with
using \code{\link[=parallelizeMCsimulation]{parallelizeMCsimulation()}}. If someone wants to parallelize the
function manually, the parameter \code{MCparam} must be set to \code{1} and
the parameter \code{SeedOption} must be changed for each iteration.

Since this package is structurally based on the \strong{"StableEstim"
package by Tarak Kharrat and Georgi N. Boshnakov}, more detailed
documentation can be found in their documentation.
}
\examples{
\donttest{
TemperedEstim_Simulation(ParameterMatrix = rbind(c(1.5,1,1,1,1,0),
                                                 c(0.5,1,1,1,1,0)),
                         SampleSizes = c(4), MCparam = 4,
                         TemperedType = "CTS", Estimfct = "ML",
                         saveOutput = FALSE)

TemperedEstim_Simulation(ParameterMatrix = rbind(c(1.5,1,1,1,1,0)),
                         SampleSizes = c(4), MCparam = 4,
                         TemperedType = "CTS", Estimfct = "GMM",
                         saveOutput = FALSE, algo = "2SGMM",
                         regularization = "cut-off",
                         WeightingMatrix = "OptAsym", t_scheme = "free",
                         alphaReg = 0.01,
                         t_free = seq(0.1,2,length.out=12))

TemperedEstim_Simulation(ParameterMatrix = rbind(c(1.45,0.55,1,1,1,0)),
                         SampleSizes = c(4), MCparam = 4,
                         TemperedType = "CTS", Estimfct = "Cgmm",
                         saveOutput = FALSE, algo = "2SCgmm",
                         alphaReg = 0.01, subdivisions = 50,
                         IntegrationMethod = "Uniform",
                         randomIntegrationLaw = "unif",
                         s_min = 0, s_max= 1)

TemperedEstim_Simulation(ParameterMatrix = rbind(c(1.45,0.55,1,1,1,0)),
                         SampleSizes = c(4), MCparam = 4,
                         TemperedType = "CTS", Estimfct = "GMC",
                         saveOutput = FALSE, algo = "2SGMC",
                         alphaReg = 0.01, WeightingMatrix = "OptAsym",
                         regularization = "cut-off", ncond = 8)
}

}
\references{
Massing, T. (2023), 'Parametric Estimation of Tempered Stable Laws'

Kim, Y. s.; Rachev, S. T.; Bianchi, M. L. & Fabozzi, F. J. (2008), 'Financial
market models with lévy processes and time-varying volatility'
\doi{10.1016/j.jbankfin.2007.11.004}

Hansen, L. P. (1982), 'Large sample properties of generalized method of
moments estimators' \doi{10.2307/1912775}

Hansen, L. P.; Heaton, J. & Yaron, A. (1996), 'Finite-Sample Properties of
Some Alternative GMM Estimators' \doi{10.1080/07350015.1996.10524656}

Feuerverger, A. & McDunnough, P. (1981), 'On the efficiency of empirical
characteristic function procedures'
\doi{10.1111/j.2517-6161.1981.tb01143.x}

Carrasco, M. & Kotchoni, R. (2017), 'Efficient estimation using the
characteristic function' \doi{10.1017/S0266466616000025};

Kuechler, U. & Tappe, S. (2013), 'Tempered stable distribution and processes'
\doi{10.1016/j.spa.2013.06.012}
}
\seealso{
\url{https://github.com/GeoBosh/StableEstim/blob/master/R/Simulation.R}
}
