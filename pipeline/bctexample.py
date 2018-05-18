import scipy.io
import bct  #available at https://github.com/aestrivex/bctpy
            #to install:
            #python setup.py --install
import numpy as np


#all functions from BCT-py is a more or less a direct port of BCT-MATLAB,
#as such all functions should exist with the same names in BCT-MATLAB
#MATLAB BCT website
#https://sites.google.com/site/bctnet/measures/list

######################################################
#Stuff I do to get the matrix into Python
######################################################

#load matrix from MATLAB
fm = scipy.io.loadmat('AALresultsROI_Condition001.mat')
#reshape such that we can pull out matrices (another MATLAB->Python fix)
fmt = np.transpose(fm['Z'])
#Patrick saved the atlas116 in the same matrices as the networks.nii
#so the matrix for one subject has both parcellation schemes,
#and has dimensions 32+116 = 84 x 148 x 148
#it looks something likes this:
# 1 .. 32    ...  33 .. 149
# 1  networks  ...   networks x atlas116
# ..
# 33
# .. atlas116 x networks ...  atlas 116
# 149 

#so we pull out the atlas116 scheme, to test for against Borchardt
fmt = fmt[:][:,32:,32:]

#the fifth guy is pretty representative of the mean 
#the first guy is also a reasonable approximate
guy = fmt[4]

#this is just Python to make sure it is in row-major order,
#matrix is symmetric so does not matter much anyway
numpy_cm = np.array(guy, order='A', dtype=float)


####################################################
#Stuff I do to clean up matrix:
#1) Drop Grey Matter column
#2) Replace NaN with 0's in diagonal (self-self)
#3) Fisher 'Z' to Pearson 'r'
#4) Remove all negative weights
####################################################

#drop the last column (the Grey matter column is not used)
sym_cm = np.delete(numpy_cm, -1, 0)
#the diagonal is filled with NaN's 
no_NaN = np.nan_to_num(sym_cm)
#convert from Fisher coefficients to Pearson coefficients (just like Borchardt)
pearson_cm = np.tanh(no_NaN)
#threshold the matrix to remove negative weights
thresh_cm = bct.threshold_absolute(pearson_cm, 0.0)


##################################################################
#Stuff I do to exract metrics once the cleaned up matrix is ready
#1) threshold lowest n% percentage of edges away. 
#   I used my own threshold_connected, which ensures fully connected
#   graph at end of function. This is not part of BCT, but I 
#   included my own code at bottom.
#2) Invert the correlation matrix, such that high correlation 
#   coefficients are interpreted as low path lengths. 
#3) Compute shortest paths for all pairs of nodes
#4) Compute characteristic path length from shortest path matrix
##################################################################

#threshold bottom n% weakest links away,
#does NOT ensure connectedness. 
#I made my own homebrew for this, BCT-MATLAB does not have it, 
#but I included it in end of this script
cm = bct.threshold_proportional(thresh_cm, 0.5)

#the matrix is inverted because high correlation corresponds to low path lengths
cm_inv = bct.invert(cm)

#the shortest path between all pairs of nodes are needed as input to charpath.
#An entry (u,v) represents the length of shortest path from node u to node v.
#distance_wei returns a tuple, the first element is the shortest path matrix
distance_wei = bct.distance_wei(cm_inv)[0]

#charpath computes the characteristic path length from the distance matrix,
#basically just by taking the average of all these shortest path lengths.
#bct.charpath returns a 5-tuple, we want the charpath-lambda, which is
#the first element
charpath = bct.charpath(distance_wei, False, False)
charpath_lambda = charpath[0]
charpath_eff = charpath[1]

#clustering coefficent seems much lower in Borchardt (~2-4% vs my ~10-20%)
#clustering coefficient gets computed locally, we take mean to get global
clustering_coef_wu = bct.clustering_coef_wu(cm)
avg_clustering_coef_wu = np.mean(clustering_coef_wu)

#result on my machine is:
#6.275372
print('Characteristic path length:')
print(round(charpath_lambda,6))

#result on my machine:
#0.191485
print('Global efficiency:')
print(round(charpath_eff,6))


#result on my machine:
#0.11283
print('Clustering coefficent:')
print(round(avg_clustering_coef_wu,6))




def threshold_connected(W, p, copy=True):
    '''
    Parameters
    ----------
    W : np.ndarray
        weighted connectivity matrix
    p : float
        proportional weight threshold (0<p<1)
    copy : bool
        if True, returns a copy of the matrix. Otherwise, modifies the matrix
        in place. Default value=True.

    Returns
    -------
    W : np.ndarray
        thresholded connectivity matrix

    Notes
    -----
    Attempts to remove bottom p % of weakest links.
    If removal of a link results in a disconnected graph (number of components > 1)
    then the edge is restored, and the algorithm continues removal of next links. 
    This function calls bct.number_of_components(), which is also included in MATLAB.
    '''

    def teachers_round(x):
         if ((x > 0) and (x % 1 >= 0.5)) or ((x < 0) and (x % 1 > 0.5)):
             return int(np.ceil(x))
         else:
             return int(np.floor(x))


    if p > 1 or p < 0:
        raise BCTParamError('Threshold must be in range [0,1]')
    if copy:
        W = W.copy()
    n = len(W)                      # number of nodes

    W = np.around(W, decimals=6)
    
    ind = np.where(W)                   # find all links



    z = np.where(W==0)
    zlen = len(z[0])

    I = np.argsort(W[ind])#[::-1]        # sort indices by magnitude

    
    ind_len = len(ind[0])

    en = int(teachers_round(ind_len * (1.0-p)))      # number of links to be discarded

    for i in range(0,en,2):                         #matrix is symmetric, so skip every other link
        temp = W[ind[0][I][i]][ind[1][I][i]]        #save the link in case for restoration

        W[ind[0][I][i]][ind[1][I][i]] = 0
        W[ind[1][I][i]][ind[0][I][i]] = 0


        if bct.number_of_components(W) > 1:
            W[(ind[0][I][i], ind[1][I][i])] = temp
            W[(ind[1][I][i], ind[0][I][i])] = temp

    return W











# import scipy.io 
# import numpy as np
# import bct

# #load in the data structure provided by Conn
# matlab_data = scipy.io.loadmat('resultsROI_Condition001.mat')
# #reshape array such that it is in row major
# mdt = np.transpose(matlab_data['Z'])
# #ensure memory layout is in row major
# numpy_cm = np.array(mdt, order='A', dtype=float)
# #drop the last column (the Grey matter column is not used)
# sym_cm = np.delete(numpy_cm, -1, 0)
# #the diagonal is filled with NaN's, replace with zeroes
# no_NaN = np.nan_to_num(sym_cm)
# #convert from Fisher coefficients to Pearson coefficients
# #(just like Borchardt)
# pearson_cm = np.tanh(no_NaN)
# #threshold the matrix to remove negative weights
# thresh_cm = bct.threshold_absolute(pearson_cm, 0.0)
