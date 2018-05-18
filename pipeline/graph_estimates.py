import bct
import numpy as np
from collections import OrderedDict


def threshold_connected(W, p, copy=True):
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


'''

def graph_estimates(cm, th):

    #Missing: number of triangles from RUBINOV
    #Missing: closeness centrality maybe covered by some distance case?
    #TODO: Anatomical and functional motifs (ALL MOTIFS: motif z-score, motif fingerprint)
    # motifs generally not used on undirected networks(Rubinov/Sporns)
    #TODO: Degree distribution
    #TODO: Average neighbour degree
    #TODO: other concepts: Small worldness, degree distribution preserving randomization
    #(done): module_degree_zscore, in particular, what to do with community affiliation vector
    #what should really be done with the binary matrices??

    #Some of the functions does not have different binary/weighted versions,
    #should we just input a binarized matrix to these?

    d = OrderedDict()

    #does not work nicely with anything
    #cm = bct.threshold_absolute(cm, 0.0)

    #dont really now how good this one works
    

 
   # cm = bct.threshold_absolute(cm, 0.0)
    #cm = bct.binarize(cm)

    #works nicely for top 50% strongest weights.
    #cm = bct.threshold_proportional(cm, 0.5)
    cm = threshold_connected(cm, th)
   # cm = bct.threshold_proportional(cm, 0.4)
    #0.45 is low as u can go

   # cm = bct.weight_conversion(cm, 'normalize')

   # cm = bct.autofix(cm)

    
    #bin_cm = np.where(cm > 0, 1, cm)
    bin_cm = bct.binarize(cm)
    
    #cm_inv = np.absolute(np.linalg.inv(cm))
    cm_inv = bct.invert(cm)

    modularity_und = bct.modularity_und(cm)
  #  modularity_und_bin = bct.modularity_und(bin_cm)

    #the community_affiliation vector that gets input to some of the functions
    community_affiliation = modularity_und[0]
    


    distance_wei = bct.distance_wei(cm_inv)
    charpath = bct.charpath(distance_wei[0], False, False)
    clustering_coef_wu = bct.clustering_coef_wu(cm)
    avg_clustering_coef_wu = np.mean(clustering_coef_wu)


    # Bug: RuntimeWarning: invalid value encountered in double_scalars
    #r = (term1 - term2) / (term3 - term2)
    #d['assortativity_bin:r'] = bct.assortativity_bin(bin_cm)
    d['assortativity_wei-r'] = bct.assortativity_wei(cm, flag=0)
    d['avg_clustering_coef_wu:C'] = avg_clustering_coef_wu

  #  d['betweenness_wei-BC'] = bct.betweenness_wei(cm_inv)
  #  d['betweenness_bin-BC'] = bct.betweenness_bin(bin_cm)

    #efficiency is different from efficiency_wei
    #efficiency_wei seems to use more sophisticated inversion
    #than just np.linalg.inv() FIXED
    d['charpath-lambda'] = charpath[0]
    #d['charpath-efficiency'] = charpath[1]   
    #d['charpath-ecc'] = charpath[2]           
    #d['charpath-radius'] = charpath[3]
    #d['charpath-diameter'] = charpath[4]

   # d['clustering_coef_bu-C'] = bct.clustering_coef_bu(bin_cm)
    d['clustering_coef_wu-C'] = clustering_coef_wu

    #d['degrees_und-deg'] = bct.degrees_und(cm)

    d['efficiency_wei-Eglob'] = bct.efficiency_wei(cm)
    d['efficiency_wei-Eloc'] = bct.efficiency_wei(cm, True)

   # d['efficiency_bin-Eglob'] = bct.efficiency_bin(bin_cm)
   # d['efficiency_bin-Eloc'] = bct.efficiency_bin(bin_cm, True)

    d['modularity_und-ci'] = modularity_und[0]
    d['modularity_und-Q'] = modularity_und[1]
    #just passing a binarized version of the matrix
  #  d['modularity_und_bin-ci'] = modularity_und_bin[0]
  #  d['modularity_und_bin-Q'] = modularity_und_bin[1]

   # d['module_degree_zscore-Z'] = bct.module_degree_zscore(cm, community_affiliation)

    #d['participation_coef'] = bct.participation_coef(cm, community_affiliation)

    d['small_worldness:S'] = compute_small_worldness(cm,
                                                     avg_clustering_coef_wu,
                                                     charpath[0])



   # d['transitivity_bu-T'] = bct.transitivity_bu(bin_cm)
    d['transitivity_wu-T'] = bct.transitivity_wu(cm)

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

    # rs = cm.shape
    # rand = np.random.rand(rs[0],rs[1])
    # rand_network = np.tril(rand) + np.tril(rand, -1).T
    # np.fill_diagonal(rand_network, 0)
    #rand_network = bct.null_model_und_sign(cm)[0]
    rand_network = bct.randmio_und_connected(cm,5)[0]

    C_rand = np.mean(bct.clustering_coef_wu(rand_network))
    while C_rand == 0.0:
        rand_network = bct.randmio_und_connected(cm,5)[0]
        C_rand = np.mean(bct.clustering_coef_wu(rand_network))

    # rand_network = bct.threshold_absolute(cm, 0.0)
    # rand_network = bct.threshold_proportional(cm, 0.2)

    #rand_network = bct.randmio_dir(cm, 5)[0]
    rand_inv = bct.invert(rand_network)
    #MODIFIED null_model_und_sign line 1009 with int(...)
    distance_rand = bct.distance_wei(rand_inv)
    charpath_rand = bct.charpath(distance_rand[0])
 


    C = cc
    #MIGHT BE 0.0 FOR LOW THRESHOLDING; PROBABLY NO CLUSTERING TO DO 
    L = cpl
    L_rand = charpath_rand[0]

 
   
    # print("This is L_rand:" + str(L_rand))
    # print("This is C_rand:" + str(C_rand))
    # print("This is C:" + str(C))

    Ctemp = C/C_rand
    Ltemp = L/L_rand

    S_W = Ctemp / Ltemp

    return S_W







# def threshold_homebrew(W, p, copy=True):
#     '''

#     ATTEMPTS TO REMOVE LINKS, PRESERVES ANY LINKS
#     WHICH WILL RESULT IN DISCONNECTED GRAPH

#     This function "thresholds" the connectivity matrix by preserving a
#     proportion p (0<p<1) of the strongest weights. All other weights, and
#     all weights on the main diagonal (self-self connections) are set to 0.

#     If copy is not set, this function will *modify W in place.*

#     Parameters
#     ----------
#     W : np.ndarray
#         weighted connectivity matrix
#     p : float
#         proportional weight threshold (0<p<1)
#     copy : bool
#         if True, returns a copy of the matrix. Otherwise, modifies the matrix
#         in place. Default value=True.

#     Returns
#     -------
#     W : np.ndarray
#         thresholded connectivity matrix

#     Notes
#     -----
#     The proportion of elements set to 0 is a fraction of all elements
#     in the matrix, whether or not they are already 0. That is, this function
#     has the following behavior:

#     That is, the 50% thresholding of x_25 does nothing because >=50% of the
#     elements in x_25 are aleady <=0. This behavior is the same as in BCT. Be
#     careful with matrices that are both signed and sparse.
#     '''
#     #W = np.float64(W)


#     if p > 1 or p < 0:
#         raise BCTParamError('Threshold must be in range [0,1]')
#     if copy:
#         W = W.copy()
#     n = len(W)                      # number of nodes
#    # np.fill_diagonal(W, 0)          # clear diagonal

#     # if np.allclose(W, W.T):             # if symmetric matrix
#     #     W[np.tril_indices(n)] = 0       # ensure symmetry is preserved
#     #     ud = 2                      # halve number of removed links
#     # else:
#     #     ud = 1

#     #W = np.tril(W)
    
#     #W[np.tril_indices(n)] = 0
#     W = np.around(W, decimals=6)
    
#     ind = np.where(W)                   # find all links

#     # z = np.where(W==0)
#     # zlen = len(z[0])

#     I = np.argsort(W[ind])#[::-1]        # sort indices by magnitude

#     # def teachers_round(x):
#     #         if ((x > 0) and (x % 1 >= 0.5)) or ((x < 0) and (x % 1 > 0.5)):
#     #             return int(np.ceil(x))
#     #         else:
#     #             return int(np.floor(x))

    

#     en = int(round((n * n) * (1-p)) - zlen)      # number of links to be discarded

#     #en = len(ind)

#    # W[W==0] = 1

# #    test = W.copy()
#     # print("ind is : " + str(ind))
#     # print("I is : " + str(I))
#     # print("en is : " + str(en))
#     # print("ind is : " + str(type(ind)))
#     # print("I is : " + str(len(I)))
#     # print("en is : " + str(en))
#     for i in range(en):
#         temp = W[ind[0][I][i]][ind[1][I][i]]
#         # print(temp)
#         # print(W[ind[0][I][i]][ind[1][I][i]])
#         # print(W[ind[1][I][i]][ind[0][I][i]])
#         # print(W[(ind[0][I][0], ind[1][I][0])])
#         # print(W[(ind[0][I][1], ind[1][I][1])])
#         #temp2 = W[(ind[0][I][i], ind[1][I][i])]
#        #  print(temp)
#        #  print(np.amin(W))
#        # # print(temp2)
#         W[ind[0][I][i]][ind[1][I][i]] = 0
#         W[ind[1][I][i]][ind[0][I][i]] = 0
#         # if ud == 2:                     # if symmetric matrix
#         #    W[:, :] = W + W.T                       # reconstruct symmetry
#        #  print(np.amin(W))
#         if bct.number_of_components(W) > 1:
#             W[(ind[0][I][i], ind[1][I][i])] = temp
#             W[(ind[1][I][i], ind[0][I][i])] = temp

#     #W[:, :] = W + W.T

#     #W[(ind[0][I][en:], ind[1][I][en:])] = 0  # apply threshold
#     # for x in np.nditer(a, order='K', op_flags=['readwrite']):
#     #     if 

#     #W[np.ix_(ind[0][I][en:], ind[1][I][en:])]=0


#     # if ud == 2:                     # if symmetric matrix
#     #     W[:, :] = W + W.T                       # reconstruct symmetry

#     return W







