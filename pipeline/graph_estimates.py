import bct #the meat of the project
import numpy as np
from collections import OrderedDict


'''
    This function "thresholds" the connectivity matrix by preserving a
    proportion p (0<p<1) of the strongest weights. 
    Entries in the connectivity matrix which are already set to 0 
    will not be considered a link. Thus only the p strongest actual 
    links will be preserved. 

    If copy is not set, this function will *modify W in place.*

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
    If removal of a link would result in a disconnected graph, 
    i.e. the graph has more than 1 component, then the link will
    be reinserted into the graph, despite having a low weight. 
    The algorithm will continue afterwards.

'''

def threshold_connected(W, p, copy=True):

    
    #rounding function which is also used in threshold_proportional,
    #probably fine to just use floor or ceiling instead
    def teachers_round(x):
         if ((x > 0) and (x % 1 >= 0.5)) or ((x < 0) and (x % 1 > 0.5)):
             return int(np.ceil(x))
         else:
             return int(np.floor(x))


    #p must be in the interval ]0;1[
    if p > 1 or p < 0:
        raise BCTParamError('Threshold must be in range [0,1]')
    if copy:
        W = W.copy()
    n = len(W)                      # number of nodes

    W = np.around(W, decimals=6)    # number of components tries to check if matrix is symmetric,
                                    # but seems to be disrupted by possible rounding errors when 
                                    # too many decimals. Hence we limit the floating point numbers
    
    ind = np.where(W)                   # find all links (entries which are not 0 (assuming all 
                                        # negative links has been dsicarded already))
    ind_len = len(ind[0])               # number of links found

    # z = np.where(W==0)                # doesn't seem to be needed, remove after testing
    # zlen = len(z[0])                  # we don't really need to know about how many zeroes our 
                                        # connectivity matrix has. 

    I = np.argsort(W[ind])#[::-1]        # sort indices by magnitude (should be ascending order)

    

    en = int(teachers_round(ind_len * (1.0-p)))   # number of links to be discarded

    #iterate until we have discarded p% of weakest links
    for i in range(0,en,2):
        #store the entry which we attempt to remove, in case of disconnectedness
        temp = W[ind[0][I][i]][ind[1][I][i]]

        #remove the link, along with it's mirror element. 
        #numpy syntax used for indexing into the connectivity matrix.
        W[ind[0][I][i]][ind[1][I][i]] = 0
        W[ind[1][I][i]][ind[0][I][i]] = 0

        #take advantage of bct already having a function to check graph connectedness by
        #checking number of components. BCT-py has this in clustering.py.
        if bct.number_of_components(W) > 1:
            #if the graph was disconnected, restore the elements we removed,
            #then continue iterations. 
            W[(ind[0][I][i], ind[1][I][i])] = temp
            W[(ind[1][I][i], ind[0][I][i])] = temp

    return W


'''
Parameters:
-----------

cm : NxN np.ndarray
     undirected weighted/binary connection matrix


Returns:
--------

d : OrderedDict
    An ordered dictionary of all the graph estimates listed in Rubinov/Sporns 2010. 
    The indices are meant to be names for the resulting .nii files (vector estimates)
    and for the column names in the .csv file (for singular values).


Notes
-----

This is the function which utilizes bctpy to extract the graph theory measures we
are interested in examining. 

'''

def graph_estimates(cm, th):

    #dictionary for storing our results
    d = OrderedDict()

    #thresholding moved here for other matrices than MatLab matrices
    #removes negative weights
    cm = bct.threshold_absolute(cm, 0.0)

    cm = threshold_connected(cm, th)

    
    #for binarizing the connectivity matrices, 
    #we work with weighted so this is turned off
    #bin_cm = bct.binarize(cm)
    
    #invert the connectivity for computing shortest paths
    cm_inv = bct.invert(cm)

    #modularity_und is found in modularity.py
    modularity_und = bct.modularity_und(cm)

    #the community_affiliation vector that gets input to some of the functions
    community_affiliation = modularity_und[0]
    
    #distance_wei and charpath is found in distance.py
    distance_wei = bct.distance_wei(cm_inv)
    charpath = bct.charpath(distance_wei[0], False, False)

    #clustering_coef_wu is found in clustering.py
    clustering_coef_wu = bct.clustering_coef_wu(cm)
    avg_clustering_coef_wu = np.mean(clustering_coef_wu)


    #assortativity_wei is found in core.py
    d['assortativity_wei-r'] = bct.assortativity_wei(cm, flag=0)

    #just taking the average of clustering_coef_wu
    d['avg_clustering_coef_wu:C'] = avg_clustering_coef_wu

    d['charpath-lambda'] = charpath[0]
    #d['charpath-efficiency'] = charpath[1]   
    #d['charpath-ecc'] = charpath[2]           
    #d['charpath-radius'] = charpath[3]
    #d['charpath-diameter'] = charpath[4]

    d['clustering_coef_wu-C'] = clustering_coef_wu


    d['efficiency_wei-Eglob'] = bct.efficiency_wei(cm)
    #d['efficiency_wei-Eloc'] = bct.efficiency_wei(cm, True)

    #d['modularity_und-ci'] = modularity_und[0]
    d['modularity_und-Q'] = modularity_und[1]

    d['small_worldness:S'] = compute_small_worldness(cm,
                                                     avg_clustering_coef_wu,
                                                     charpath[0])

   
   #transitivity_wu can be found in clustering.py
    d['transitivity_wu-T'] = bct.transitivity_wu(cm)


    #EXAMPLES for local measures and binary measures. Comment in to use. 

    #VECTOR MEASURES
    #d['betweenness_wei-BC'] = bct.betweenness_wei(cm_inv)
    # d['module_degree_zscore-Z'] = bct.module_degree_zscore(cm, community_affiliation)
    #d['degrees_und-deg'] = bct.degrees_und(cm)
    #d['charpath-ecc'] = charpath[2]


    #BINARIES
    # d['clustering_coef_bu-C'] = bct.clustering_coef_bu(bin_cm)
    # d['efficiency_bin-Eglob'] = bct.efficiency_bin(bin_cm)
    # d['efficiency_bin-Eloc'] = bct.efficiency_bin(bin_cm, True)
    #  d['modularity_und_bin-ci'] = modularity_und_bin[0]
    #  d['modularity_und_bin-Q'] = modularity_und_bin[1]
    # d['transitivity_bu-T'] = bct.transitivity_bu(bin_cm)
    #  d['betweenness_bin-BC'] = bct.betweenness_bin(bin_cm)
    #  modularity_und_bin = bct.modularity_und(bin_cm)
    #d['participation_coef'] = bct.participation_coef(cm, community_affiliation)


    ######## charpath giving problems with ecc, radius and diameter
    # np.seterr(invalid='ignore')


    return d




'''
Parameters
----------

cm : NxN np.nadarray
     undirected binary/weighted connection matrix

cc : float
     Clustering coefficient of the connectitvity matrix
     Takes this as input as it is computed in graph_estimates()

cpl : float
      Characteristic path length of the ocnnectivity matrix.
      Also takes this as input rather than computing it again.


Returns:
--------

S_W : float
      the Small-Worldness value of the network tested against a random
      network also of size NxN, as explained in Rubinov/Sporns 2010.


Notes
-----

Very low scores of S seems to be reported at the moment, while RS10 reports networks which
does exhibit small-worldness has S >> 1 (far greater than). This matrix seems to have very low
small-worldness property, which might be true actually. 


'''


def compute_small_worldness(cm, cc, cpl):

    #randmio_und_connected can be found in reference.py
    #second argument is number of iterations
    #construct a random network for comparison with our real network
    rand_network = bct.randmio_und_connected(cm,5)[0]

    #clustering_coef_wu is found in clustering.py
    #make sure that C_rand is non-zero, so we avoid division with zero
    #could probably be made more correct by taking the average of
    #some number of random networks. 
    #we did not do this to keep run time at a minimum
    C_rand = np.mean(bct.clustering_coef_wu(rand_network))
    while C_rand == 0.0:
        rand_network = bct.randmio_und_connected(cm,5)[0]
        C_rand = np.mean(bct.clustering_coef_wu(rand_network))

    #invert can be found in other.py
    rand_inv = bct.invert(rand_network)
    
    #distance_wei and charpath can be found in distance.py
    distance_rand = bct.distance_wei(rand_inv)
    charpath_rand = bct.charpath(distance_rand[0])
 
    #compute the small worldness index according to Rubinov
    C = cc
    L = cpl
    L_rand = charpath_rand[0]

    Ctemp = C/C_rand
    Ltemp = L/L_rand

    S_W = Ctemp / Ltemp

    return S_W






