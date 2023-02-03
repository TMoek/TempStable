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
  subordinator (15:49)

  This note appears only for devtools::check_rhub().
  The term is used like this in the research field.
  
  
* checking for detritus in the temp directory ... NOTE: Found the following
  files/directories: 'lastMiKTeXException'
  
  This note appears only for devtools::check_rhub() (only on Windows, not Linux).
  This note seems common on rhub windows (based on Github search which
  turns up many cran-comments mentioning this).
  Once, I had this error with devtools::check_mac_release(), too. Then I deleted
  the function “\mathbb” in “\deqn” in my documentation for charTSS() and
  charCTS() and the note was gone. This notes seems to be related to using 
  mathematical functions in the documentation.
  

* While Check works with Windows (local; devtools::check_win_release(); 
  rhub Windows Server 2022) and on Mac (devtools::check_mac_release()), it works
  only partly on Linux (rhub Ubuntu Linux 20.04.1 LTS, R-release, GCC; rhub 
  Fedora Linux, R-devel, clang, gfortran). On Fedora Linux Check is running into
  an error ("Running ‘testthat.R’Build timed out (after 20 minutes). Marking the
  build as failed.") due to a building time out. On Ubuntu Linux running tests 
  last shorter ("Running ‘testthat.R’ [5m/18m] [5m/18m] OK") as on all other 
  platforms. I have tried to shorten the test cases. However, shortening them 
  further would now mean deleting test cases.


* This is a new release.

File-System:
Since tempered stable distributions are based on normal stable distributions, 
we use many code elements from the package 'stableEstim'. The calculation in our 
Monte Carlo simulation (MonteCarloFunction.TemperedEstim_Simulation()) extends 
the function 'StableEstim::Estim_Simulation()' in a statistical way, but is 
structurally similar. In the function, the user's file system is written to 
during each iteration of the simulation. 


## Downstream dependencies

First version, no downstream dependencies
