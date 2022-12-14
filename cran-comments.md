## Test environment
* local Windows 8.1, R 4.2.2
* using devtools::check_win_release(), R 4.2.2
* using devtools::check_mac_release(), R 4.2.1

## R CMD check results

0 errors | 0 warnings | 1 note

* 1 Note: "Authors@R field gives persons with deprecated elements:
           Warning in person1(given = given[[i]], family = family[[i]], 
           middle = middle[[i]],  :
           It is recommended to use 'given' instead of 'middle'."
           
           This note only occurs with check_win_release(). According to 
           internet research, this is a problem that only occurs with Windows 
           servers and can be ignored. 
           

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
