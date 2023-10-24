#### Update TempStable 0.2.1

## Test environment

* local Windows 11, R 4.3.0
* using devtools::check(remote = TRUE, manual = TRUE)
* using devtools::check_win_devel(), R 4.4
* using revdepcheck::revdep_check(num_workers = 4)

## R CMD check results
There were no ERRORs or WARNINGs.

There was 1 new NOTE:

*   Maintainer: 'Till Massing <till.massing@uni-due.de>'
  
  New maintainer:
    Till Massing <till.massing@uni-due.de>
  Old maintainer(s):
    Cedric Maximilian Juessen <cedric.juessen@vwl.uni-due.de>
  
  Package has a VignetteBuilder field but no prebuilt vignette index.
  
This is due to the change of maintainer and was expected.

There was 1 NOTE, that has been there also with the previous versions:

* checking for detritus in the temp directory ... NOTE
  Found the following files/directories:
  'lastMiKTeXException'

## revdepcheck results

We checked 0 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages


#### Update TempStable 0.2.0

## Test environment
* local Windows 8.1, R 4.3.0 || Windows 8.1 R 4.2.2
* using devtools::check_win_release()
* using devtools::check_win_oldrelease()
* using devtools::check_win_devel(), R 4.4
* using devtools::check_mac_release()
* using devtools::check_rhub(), which uses 3 different engines: 1 windows and 2
  linux.

## R CMD check results
There were no ERRORs or WARNINGs.

There were 2 new NOTEs:

* checking CRAN incoming feasibility ... [8s/28s] NOTE
  Maintainer: ‘Cedric Maximilian Juessen <cedric.juessen@vwl.uni-due.de>’
  Possibly misspelled words in DESCRIPTION:
  Bianchi (27:3, 28:3)
  Hofert (33:3)
  
  Those words are names.
  
* checking for non-standard things in the check directory ... NOTE
  Found the following files/directories: ''NULL''
  
  This note only appears with Windows. I did not find anything about it on the 
  internet and could not fix it.
  
  
There were 3 NOTEs, that have been there also with the previous versions:

* checking HTML version of manual ... NOTE
  Skipping checking HTML validation: no command 'tidy' found
  
  This Note occurs only for Fedora Linux, R-devel, clang, gfortran. 
  The internet says that this note is due to a package being missing on the 
  virtual machine.

* checking for detritus in the temp directory ... NOTE
  Found the following files/directories:
  'lastMiKTeXException'
  
  This note appears only for devtools::check_rhub() (only on Windows, not Linux).
  This note seems common on rhub windows (based on Github search which
  turns up many cran-comments mentioning this).
  Once, I had this error with devtools::check_mac_release(), too. Then I deleted
  the function “\mathbb” in “\deqn” in my documentation for charTSS() and
  charCTS() and the note was gone. This notes seems to be related to using 
  mathematical functions in the documentation.
  
## Downstream dependencies

There are no downstream dependencies yet (07.07.2023). 
Used tools::dependsOnPkgs(c("TempStable"), "all").

#### Update TempStable 0.1.1

## The basis of this submission is a request for amendment by CRAN with a due date of 06.04.23
The following points were criticised in TempStable 0.1.0, which I have dealt 
with as follows:
*	Errors for 3 of 13 test environments

  All errors originated in the test cases, where rounding inaccuracies occurred
  due to different functions being called. All test cases that were affected 
  have been adjusted.
  
* It is inadvisable to use a dependence on R with patchlevel (the third digit) 
  other than zero
  
  As far as it was possible, it was changed
  
* "methods" was part of the imports in DESCRIPTION file

  Moved to "Depends"


#### New submission 0.1.0

## Test environment
* local Windows 8.1, R 4.2.2
* using devtools::check_win_release(), R 4.2.2 and devtools::check_win_devel()
* using devtools::check_mac_release(), R 4.2.1
* using devtools::check_rhub(), which uses 3 different engines: 1 windows and 2
  linux.

## R CMD check results
There were no ERRORs or WARNINGs.

There were 4 NOTEs:

* checking CRAN incoming feasibility ... NOTE
  Maintainer: ‘Cedric Maximilian Juessen <cedric.juessen@vwl.uni-due.de>’
  Days since last update: 6
  
  The basis of this submission is a request for amendment by CRAN with a due 
  date of 06.04.23. All changes refer exclusively to this request.
  
* Found the following (possibly) invalid DOIs:
  DOI: 10.1002/9781118268070
    From: DESCRIPTION
    Status: Forbidden
    Message: 403
  DOI: 10.1111/j.2517-6161.1981.tb01143.x
    From: DESCRIPTION
    Status: Forbidden
    Message: 403
  DOI: 10.2307/1912775
    From: DESCRIPTION
    Status: Forbidden
    Message: 403

  This Note occurs only for Ubuntu Linux 20.04.1 LTS, R-release, GCC. 
  All DOIs work fine.  I have double checked them. This problem could be due to 
  the virtual machine not having access to the documents.

* checking HTML version of manual ... NOTE
  Skipping checking HTML validation: no command 'tidy' found
  
  This Note occurs only for Fedora Linux, R-devel, clang, gfortran. 
  The internet says that this note is due to a package being missing on the 
  virtual machine.

* checking for detritus in the temp directory ... NOTE: Found the following
  files/directories: 'lastMiKTeXException'
  
  This note appears only for devtools::check_rhub() (only on Windows, not Linux).
  This note seems common on rhub windows (based on Github search which
  turns up many cran-comments mentioning this).
  Once, I had this error with devtools::check_mac_release(), too. Then I deleted
  the function “\mathbb” in “\deqn” in my documentation for charTSS() and
  charCTS() and the note was gone. This notes seems to be related to using 
  mathematical functions in the documentation.
  

## Downstream dependencies

There are no downstream dependencies yet.
