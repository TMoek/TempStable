## Resubmission
The following points were criticised in the last submission, which I have dealt 
with as follows:

* If there are references describing the methods in your package, please add 
  these in the description field of your DESCRIPTION file...

  Added all references mentioned in any exported functions to DESCRIPTION-file.
  
* You write information messages to the console that cannot be easily 
  suppressed... e.g.: R/MasterEstim.R
  
  In the R/MonteCarloFunction/TemperedEstim_Simulation() for the parameter algo 
  == “Cue…” && “IT…”, a sub-iteration control and results of the simulation were
  output to the console. This is not happening anymore (I changed PrintIter to
  “FALSE” in GMM_TSS/checkIterationControl()).
  
  In TemperedEstim() for TemperedType == “Normal” 2 lines were written to the 
  console. This happened in MethodOfCumulants_NTS/MoC_NTS() due to the external 
  function rootSolve::multiroot(). Now the output is suppressed by 
  capture.output(…).
  
  R/MonteCarloFunction/TemperedEstim_Simulation() also writes progress 
  information to console (remaining time and remaining iterations). This is 
  helpful for the user, as simulations cost long time. Therefore, it should stay 
  in the package.

* Please do not modifiy the .GlobalEnv... e.g.: R/MonteCarloFunction.R

  The global environment was modified in parallelizeMCsimulation() with the 
  function parallel::clusterExport(). Now the exported functions are directly 
  called in foreach::foreach(.export = (…)). ClusterExport() will not be used 
  anymore.

* Please ensure that your functions do not write by default in the user’s 
  filesystem:
  
  TemperedEstim_Simulation() was writing to file system for each iteration of 
  the Monte Carlo Simulation. As those simulation can take a lot of time, this 
  file was used to restart the simulation at a specific point when Rstudio 
  crashed. Now this function does not create a checkpoint file on the filesystem
  anymore and the values for the next iteration will be passed by processing 
  environment. Therefore, I changed the following things:
    - I added some start values for different loops in 
    TemperedEstim_Simulation(). This values are added in readCheckPoint() as 
    parameters.
    - The function writeCheckPoint() in TemperedEstim_Simulation() will no 
    longer be called.
    - readCheckPoint() is changed completely. Now this function only returns 
    values as a list which it got as input. 
    - deleteCheckPoint() does nothing anymore. 
    
  ParallelizeMCsimulation() wrote iteration control text file to filesystem by 
  default.
    - I have created a parameter that is "False" by default and that determines 
    whether a text file is created on the file system.
    
* Please ensure that you do not use more than 2 cores in your examples, 
  vignettes, etc.
  
  parallelizeMCsimulation() is the only function which is able to use more than 
  one core. For this function no examples exist. I changed the default values 
  to 2 cores and changed the description.


## Test environment
* local Windows 8.1, R 4.2.2
* using devtools::check_win_release(), R 4.2.2 and devtools::check_win_devel()
* using devtools::check_mac_release(), R 4.2.1
* using devtools::check_rhub(), which uses 3 different engines: 1 windows and 2
  linux.

## R CMD check results
There were no ERRORs or WARNINGs.

There were 2 NOTEs:

* checking CRAN incoming feasibility: Possibly misspelled words in DESCRIPTION:
  Carrasco (24:3)
  Feuerverger (25:3)
  Kawai (28:3)
  Kuechler (30:3)
  Masuda (28:11)
  Rachev (31:3)
  Tappe (30:14)
  al (26:13, 29:10, 31:13)
  et (26:10, 29:7, 31:10)
  subordinator (15:49)

  This note appears only for devtools::check_rhub().
  Terms are names and other terms that are used like this in the research field.
  
  
* checking for detritus in the temp directory ... NOTE: Found the following
  files/directories: 'lastMiKTeXException'
  
  This note appears only for devtools::check_rhub() (only on Windows, not Linux).
  This note seems common on rhub windows (based on Github search which
  turns up many cran-comments mentioning this).
  Once, I had this error with devtools::check_mac_release(), too. Then I deleted
  the function “\mathbb” in “\deqn” in my documentation for charTSS() and
  charCTS() and the note was gone. This notes seems to be related to using 
  mathematical functions in the documentation.
  

* First submission:
  While Check works with Windows (local; devtools::check_win_release(); 
  rhub Windows Server 2022) and on Mac (devtools::check_mac_release()), it works
  only partly on Linux (rhub Ubuntu Linux 20.04.1 LTS, R-release, GCC; rhub 
  Fedora Linux, R-devel, clang, gfortran). On Fedora Linux Check is running into
  an error ("Running ‘testthat.R’Build timed out (after 20 minutes). Marking the
  build as failed.") due to a building time out. On Ubuntu Linux running tests 
  last shorter ("Running ‘testthat.R’ [5m/18m] [5m/18m] OK") as on all other 
  platforms. I have tried to shorten the test cases. However, shortening them 
  further would now mean deleting test cases.
  
  First re-submission:
  check_rhub worked for all 3 environments.


* This is a new release.


## Downstream dependencies

First version, no downstream dependencies
