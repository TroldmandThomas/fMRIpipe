# fMRIpipe
Python pipeline for automating rs-fMRI graph theory estimates.  
The pipeline assumes preprocessing and correlation matrix construction was performed in the MATLAB module 'Conn'.  
First step is to load the matrices from MATLAB into a numpy ndarray variable in Python,
and then use the [BrainConnectivityToolbox](https://github.com/aestrivex/bctpy) written in Python to estimate graph theory measures. 
Then we use packages from [scipy.stats](https://docs.scipy.org/doc/scipy/reference/stats.html) to perform t-tests, ks-tests and mannwhitney u-tests.

This software pipeline assumes that Python3+ is installed, in particular we used Python3.6 for all the work regarding the project. 


## VirtualEnvironment
It is recommended to use the virtualenvironment (venv) provided by Python3+.  
There is an excellent primer on venv's for Python at [realpython](https://realpython.com/python-virtual-environments-a-primer/),
if you want to know more about them. 
To setup a venv, simply type:
> python3.6 -m venv myenv  

where 'myenv' can be any name that you choose, to be the name of the venv.

Once the venv is setup, it may be activated like so:

>source fmrienv/bin/activate

A little parenthesis with (myenv) will show up in the terminal once the venv has been activated.  
With the venv activated, continue with the rest of the steps outlined in this document. 

To exit the venv, simply type:

>deactivate

## Installation

Once the GitHub folder containing the project has been downloaded, navigate to the folder and
run the _pipeinstall.sh_ shellscript. This will install all the 
On a Mac, typing the following would be sufficient:

>sh pipeinstall.sh

The installation will probably take a minute or two. 

## Data

The pipeline was built for MATLAB files following the 'Conn' module file structure. As such, it has been built for files
containing correlation matrices such as _resultsROI_Condition001.mat_ . Additionally, it is required to have the corresponding
group identification labels, which should contain at least subject status and scan season. This file could for example be named _groupID.csv_ .

Experimental work has been done using pure NumPy array files, but there are some issues regarding these. We only tested with a very low sample size, which might be the reason the statistical tests in this regard are returning errors. 


## Usage

A control script for the whole pipeline can be found in _entry.py_ . It has five modes:

1. 'full' (runs the whole pipeline)
2. 'estimate (runs only the graph theory estimates)
3. 'ttest' (runs only the t-tests and u-tests)
4. 'graphs' (runs both t-tests, u-tests and draws graphs(plots!) based on these tests
5. 'numpy' (experimental numpy mode, runs graph theory estimates only on provided numpy arrays)




















