## Test environment
* local Windows 8.1, R 4.2.2
* using devtools::check_win_release(), R 4.2.2 and devtools::check_win_devel()
* using devtools::check_mac_release(), R 4.2.1

## R CMD check results

0 errors | 0 warnings | 0 note

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
