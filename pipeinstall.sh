#make a virtual environment for python3.6
#python3.6 -m venv uenv

#activate the virtual environment
#source uenv/bin/activate

#upgrade pip to latest version
pip install --upgrade pip

#install packages
pip install numpy
pip install nibabel
pip install scipy
pip install matplotlib
pip install pandas

#newer version of rpy2 seems to be buggy(on Mac at least)
#only used for glm.py
pip install rpy2==2.8.6
#use these commands if getting wierd warnings from running glm.py
# export LANG=en_US.UTF-8
# export LC_ALL=en_US.UTF-8
# locale

#get the bctpy from GitHub 
#(note that pip will get an outdated version, even though both are version 0.5)
git clone https://github.com/aestrivex/bctpy.git

#install bctpy
cd bctpy
python3.6 setup.py install
cd ..


#deactivate the virtual environment
#deactivate