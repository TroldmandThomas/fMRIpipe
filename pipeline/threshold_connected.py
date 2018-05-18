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
        #for example, ind[0][I][i] should read:
        # the i'th element in I represents the i'th smallest link in ind[0], 
        # where ind[0] is a row of W. Similarly for ind[1].
        # so W[ind[0][I][i]][ind[1][I][i]] is the i'th smallest link in W. 
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