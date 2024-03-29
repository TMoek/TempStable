% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MonteCarloFunction.R
\name{parallelizeMCsimulation}
\alias{parallelizeMCsimulation}
\title{Function to parallelize the Monte Carlo Simulation}
\usage{
parallelizeMCsimulation(
  ParameterMatrix,
  MCparam = 10000,
  SampleSizes = c(200),
  saveOutput = FALSE,
  cores = 2,
  SeedOptions = NULL,
  iterationDisplayToFileSystem = FALSE,
  ...
)
}
\arguments{
\item{ParameterMatrix}{The matrix is to be composed of vectors, row by row.
Each vector must fit the pattern of theta of the \code{TemperedType}.
Compared to the function \code{\link[=TemperedEstim_Simulation]{TemperedEstim_Simulation()}}, the matrix here may
contain only one parameter vector.}

\item{MCparam}{Number of Monte Carlo simulation for each couple of parameter,
default=100; integer}

\item{SampleSizes}{Sample sizes to be used to simulate the data. By default,
we use 200 (small sample size). Vector of integer. Compared to the function
\code{\link[=TemperedEstim_Simulation]{TemperedEstim_Simulation()}}, the vector here may contain only one integer.}

\item{saveOutput}{Logical flag: In the function \code{\link[=TemperedEstim_Simulation]{TemperedEstim_Simulation()}}
the argument can be true. Then an external csv file is created. Here the
argument must be false. The output of the values works in this function
exclusively via the return of the function.}

\item{cores}{size of cluster for parallelization. Positive Integer.}

\item{SeedOptions}{is an argument what can be used in
\code{\link[=TemperedEstim_Simulation]{TemperedEstim_Simulation()}} but must be NULL here.}

\item{iterationDisplayToFileSystem}{creates a text file in your file system
that displays the current iteration of the simulation.}

\item{...}{The function works only if all necessary arguments from the
function \code{\link[=TemperedEstim_Simulation]{TemperedEstim_Simulation()}} are passed. See description and
details.}
}
\value{
The return object is a list of 2. Results of the simulation are
listed in \code{$outputMat}.
}
\description{
Since the Monte Carlo Simulation is very computationally intensive, it may
be worthwhile to split it across all available processor cores. To do this,
simply pass all the parameters from the \code{\link[=TemperedEstim_Simulation]{TemperedEstim_Simulation()}}
function to this function in the same way.
}
\details{
In this function exactly the arguments must be passed, which are also needed
for the function \code{\link[=TemperedEstim_Simulation]{TemperedEstim_Simulation()}}. However, a few functions of
\code{\link[=TemperedEstim_Simulation]{TemperedEstim_Simulation()}} are not possible here. The restrictions are
described in more detail for the individual arguments.

In addition to the arguments of function \code{\link[=TemperedEstim_Simulation]{TemperedEstim_Simulation()}}, the
argument "cores" can be assigned an integer value. This value determines how
many different processes are to be parallelized. If value is \code{NULL}, R
tries to read out how many cores the processor has and passes this
value to "cores".

During the simulation, the progress of the simulation can be viewed in a
file in the workspace named "IterationControlForParallelization.txt".
}
