# fMRIpipe
Python pipeline for automating rs-fMRI graph theory estimates.  
The pipeline assumes preprocessing and correlation matrix construction were performed in the MATLAB module 'Conn'.  
First step is to load the matrices from MATLAB into a numpy ndarray variable in Python,
and then use the [BrainConnectivityToolbox](https://github.com/aestrivex/bctpy) written in Python to estimate graph theory measures.
The MATLAB version of BCT has excellent documentation about the behavior, inputs and outputs at [BCT-MATLAB](https://sites.google.com/site/bctnet/measures/list). BCT-PY is a direct port of BCT-MATLAB, so this documentation is often very helpful.   
Then we use packages from [scipy.stats](https://docs.scipy.org/doc/scipy/reference/stats.html) to perform t-tests, ks-tests and mannwhitney u-tests. Based upon these tests and data, plots can be produced through [matplotlib](https://matplotlib.org).

## Prerequisites
This software pipeline assumes that Python3+ is installed, in particular we used Python3.6 for all the work regarding the project.   
Testing performed with other versions of Python has NOT been carried out.

## VirtualEnvironment
It is recommended to use the virtualenvironment (venv) provided by Python3+.  
Venv's are a way to manage Python packages, they are held in a container, thus the user only has access to packages within the venv. This is useful to simulate a clean installation of the packages, and guards against conflicting packages being installed. 
There is an excellent primer on venv's for Python at [realpython](https://realpython.com/python-virtual-environments-a-primer/),
if you want to know more about them. 
To setup a venv, simply type:
> python3.6 -m venv myenv  

where 'myenv' can be any name that you choose, to be the name of the venv.

Once the venv is setup, it may be activated like so:

>source myenv/bin/activate

A little parenthesis with (myenv) will show up in the terminal once the venv has been activated.  
With the venv activated, continue with the rest of the steps outlined in this document. 

To exit the venv, simply type:

>deactivate

## Installation

Once the GitHub folder containing the project has been downloaded, navigate to the folder and
run the **pipeinstall.sh** shellscript. This will install all the packages required automatically.
On a Mac, typing the following would be sufficient:

>sh pipeinstall.sh

The installation will probably take a minute or two. 

## Packages used

A list of the packages installed from the **install.sh**. This is equivalent to a **requirements.txt**, but note that the BCT-PY one can obtain through **pip** is outdated, it is recommend to download from [the authors GitHub page](https://github.com/aestrivex/bctpy) instead. Even though that they list the same version (0.5.0), they are _not_ the same package. 

* bctpy==0.5.0
* cycler==0.10.0
* kiwisolver==1.0.1
* matplotlib==2.2.2
* nibabel==2.2.1
* numpy==1.14.3
* pandas==0.23.0
* pyparsing==2.2.0
* python-dateutil==2.7.3
* pytz==2018.4
* rpy2==2.8.6
* scipy==1.1.0
* six==1.11.0


## Data

The pipeline was built for MATLAB files following the **Conn** module file structure. As such, it has been built for files
containing correlation matrices such as **resultsROI_Condition001.mat** . Additionally, it is required to have the corresponding
group identification labels, which should contain at least subject status and scan season. This file could for example be named **groupID.csv** .


## Usage

A control script for the whole pipeline can be found in **entry.py** . It has four modes:

1. 'estimate' (runs only the graph theory estimates)
2. 'ttest' (runs only the t-tests and u-tests)
3. 'plots' (runs both t-tests, u-tests and draws plots based on these tests
4. 'full' (runs the whole pipeline)

### estimate

This mode will only carry out the estimation of graph theory measures on the given dataset.
Assuming the subject file to be estimated is named **resultsROI_Condition001.mat**, 
the group file is labeled **groupID_Thomas.csv**,
the various thresholds to be estimated upon are from 40 to 42, with a 2 percent increase for each iteration,
and the directory to output the resulting graph theory estimate CSV files should be written to,
then the following command can be used:

>python3.6 entry.py estimate -mat resultsROI_Condition001.mat -id groupID_Thomas.csv -thr 40:42:2 -out ~/Desktop/PipeTest

Only estimate files are produced from this step, which are placed under the **auto_results** directory, with the naming convention **estimate.xx.csv**, where 'xx' denote the threshold percentage. 
This could be useful if one wishes to add or edit estimate CSV files, that later has to be tested once the user is ready for it. 

### ttest

The statistical results from the t-tests are also saved in CSV files under the **tests** folder. 
The following command can be used to run the statistical tests:

>python3.6 entry.py ttest -ws 'W' -dir ~/Desktop/PipeTest/auto_results -out ~/Desktop/PipeTest/

where **-ws** denotes the season ('W' for winter, 'S' for summer), **-dir** denotes the path to the files that should be testet (i.e. the estimate files obtained from running in _estimate_ mode) and **-out** is the path to where the resulting CSV files with the _p_-values should be written to. A folder named **tests** is created at the given path by -out. Within **tests**, two CSV files are created: **W_normality.csv** (which contains results of KS-tests) and **W_ttests.csv** (which contains the results of the two sample t-tests).  
This mode also prints the various results to the screen when run. 

### plots

Producing plots from the data estimated by the _estimate_ mode is also supported. Try out:

>python3.6 entry.py plots -dir ~/Desktop/EstimateTest/auto_results -out ~/Desktop/PipeTest/

where **-dir** denotes the path to the files that should be testet and **-out** is the path to where the resulting plots should be written to.
This will create a folder named **plots**, and write the plots to this folder with the naming convention **W_assortativity_wei-r.png**.
A plot will be drawn for each graph theory measure, and for each season. As such, a file named **S_assortativity_wei-r.png** will also be produced during this execution. 

### full

All the above modes can be combined into a single command:

>python3.6 entry.py full -mat resultsROI_Condition001.mat -id groupID_Thomas.csv -out ~/Desktop/FullTest/ -thr 24:26:2

This will run the graph theory estimation, apply statistical testing, and finally draw some plots based upon said statistical tests. Notice that this command takes exactly the same inputs as the _estimate_ mode. It will produce three folders, **auto_results**, **tests** and **plots** in the given path by **-out**. In these three folders, files will be placed as described previously. 

### optional clause: -cut

The graph theory estimate modes for _estimate_ and _full_ also have an additional, optional clause: **-cut**. This will take a specified subset of the matrix, and only use this in the graph theory estimations. It is useful if multiple correlation matrices are stored in the same file. For example, if a user only wanted to use the first _32x32_ indices of a given matrix, one could run the pipeline with:

>python3.6 entry.py full -mat resultsROI_Condition001.mat -id groupID.csv -thr 40:42:2 -cut 1:32x1:32

Note that the option assumes one-based indexing is used, this is to adhere to the MATLAB array indexing convention.

## Notes

As a final note, all the files can in the pipeline can also be used individually, like regular python scripts. This is the only way to run the **glm.py**, which contains our generalized linear models (this code depends on **R** being installed, and is imported into Python by the **rpy2** module). 


































