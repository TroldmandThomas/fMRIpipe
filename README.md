# fMRIpipe
Python pipeline for automating rs-fMRI graph theory estimates.  
The pipeline assumes preprocessing and correlation matrix construction was performed in the MATLAB module 'Conn'.  
First step is to load the matrices from MATLAB into a numpy ndarray variable in Python,
and then use the [BrainConnectivityToolbox](https://github.com/aestrivex/bctpy) written in Python to estimate graph theory measures. 
Then we use packages from [scipy.stats](https://docs.scipy.org/doc/scipy/reference/stats.html) to perform t-tests and mannwitney u-tests.


## VirtualEnvironment
It is recommended to use the virtualenvironment provided by Python3+.  
To do so, simply type 
> python3.6 -m venv myenv  

where 'myenv' can be any name that you choose, to be the name of the virtualenvironment.

