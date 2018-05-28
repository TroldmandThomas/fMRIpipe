# fMRIpipe
Python pipeline for automating rs-fMRI graph theory estimates.  
The pipeline assumes preprocessing and correlation matrix construction were performed in the MATLAB module 'Conn'.  
First step is to load the matrices from MATLAB into a numpy ndarray variable in Python,
and then use the [BrainConnectivityToolbox](https://github.com/aestrivex/bctpy) written in Python to estimate graph theory measures.
The MATLAB version of BCT has excellent documentation about the behavior, inputs and outputs at [BCT-MATLAB](https://sites.google.com/site/bctnet/measures/list). BCT-PY is a direct port of BCT-MATLAB, so this documentation is often very helpful. 
Then we use packages from [scipy.stats](https://docs.scipy.org/doc/scipy/reference/stats.html) to perform t-tests, ks-tests and mannwhitney u-tests.

This software pipeline assumes that Python3+ is installed, in particular we used Python3.6 for all the work regarding the project. 


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

The pipeline was built for MATLAB files following the 'Conn' module file structure. As such, it has been built for files
containing correlation matrices such as **resultsROI_Condition001.mat** . Additionally, it is required to have the corresponding
group identification labels, which should contain at least subject status and scan season. This file could for example be named **groupID.csv** .

Experimental work has been done using pure NumPy array files, but there are some issues regarding these. We only tested with a very low sample size, which might be the reason the statistical tests in this regard are returning errors. 


## Usage

A control script for the whole pipeline can be found in **entry.py** . It has five modes:

1. 'full' (runs the whole pipeline)
2. 'estimate' (runs only the graph theory estimates)
3. 'ttest' (runs only the t-tests and u-tests)
4. 'graphs' (runs both t-tests, u-tests and draws graphs(plots!) based on these tests
5. 'numpy' (experimental numpy mode, runs graph theory estimates only on provided numpy arrays)

### full

Assuming the subject file to be estimated is named **resultsROI_Condition001.mat**, 
the group file is labeled **groupID.csv** and
the various thresholds to be estimated upon are from 40 to 42, with a 2 percent increase for each iteration,
the following command can be used:

>python3.6 entry.py full -mat resultsROI_Condition001.mat -id groupID.csv -thr 40:42:2

This will run the graph theory estimation, apply statistical testing, and finally draw some graphs/plots based upon said statistical tests. To view the produced plots, navigate to the **graphs** directory and open up one of **.png** images and have a look.

### estimate

For running only the graph theory estimates, try out the following command:

>python3.6 entry.py estimate -mat resultsROI_Condition001.mat -id groupID.csv -thr 60:62:2

Only estimate files are produced from this step, which are placed under the **auto_results** directory, with the naming convention **estimate.60.csv**. This could be useful if one wishs to add or edit estimate CSV files, that later has to be tested once the user is ready for it. 

### graphs

Producing graphs/plots can also be done seperately. Try out:

>python3.6 entry.py graphs

Notice that this does not take any additional arguments. The path to the estimate files are hardcoded in, and simply performs statistical testing and plotting based on the **estimate.xx.csv** contained in the **auto_results** folder. This mode will also automatically apply statistical testing for both the summer and winter season. 

### ttest

The statistical results from the t-tests and u-tests are also saved in CSV files under the **tests** folder. However, the code is a bit messy and so is the resulting CSV files. This will be addressed in a future build. But for now, a simpler overview can be obtained by just running the **entry.py** script with the mode **ttest**. Try out:

>python3.6 entry.py ttest

This will print all the statistical result to the terminal. Not a great solution, but might be a better overview than the CSV files. 

### optional clause: -cut

The graph theory estimate modes also have an additional, optional clause: -cut. This will take a specified subset of the matrix, and only use this in the graph theory estimations. It is useful if multiple correlation matrices are stored in the same file. For example, if a user only wanted to use the first 32x32 indices of a given matrix, one could run the pipeline with:

>python3.6 entry.py full -mat resultsROI_Condition001.mat -id groupID.csv -thr 40:42:2 -cut 1:32x1:32

Note that the option assumes one based indexing is used, this is to adhere to the MATLAB array indexing convention.

### optional clause: matrices and thresholds provided as lists

It is also possible to provide multiple matrices in one go, and to provide a list of thresholds rather than a range (start":"end":"stride) notation. This could be typed as:

>python3.6 entry.py numpy -mat [resultsROI_Condition001_Subject001.mat,resultsROI_Condition001_Subject001.mat] -id NumPy_test/numpygroup.csv -thr [0.10,0.12,0.14,0.16,0.18,0.20]

The lists share a similiar notation with pure Python lists. 

## Notes

As a final note, all the files can in the pipeline can also be used individually, like regular python scripts. This is the only way to run the **glm.py**, which contains our genralized linear models (this code depends on **R** being installed, and is imported into Python by the **rpy2** module). 


































