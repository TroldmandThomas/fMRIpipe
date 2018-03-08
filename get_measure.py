import bct

'''
modularity_und(A, gamma=1.0, kci=None)
'''


'''
Parameters:
-----------

measure : string
          Name of the measure to extract.
          Remember to state weighted/binary,
          directed/undirected,
          Ex. "modularity_und"
          EXTEND COMMENT
cm : NxN np.ndarray
     undirected weighted/binary connection matrix


Notes
-----
try to make the dictionary for lambda funtioncs global instead of local 
to vec_measure     
'''


def vec_measure(measure, cm):

    vec_res = {}
                                  # [0] is because of modularity_und returning more than one result value
    m_dict = {'modularity_und' : lambda x: bct.modularity_und(x),
              'clustering_coef_wu' : lambda x: bct.clustering_coef_wu(x),
              'efficiency_wei' : lambda x: bct.efficiency_wei(x)} #second argument has 'defualt: local=False', hence why float is returned

    #vec_res[measure] = m_dict[measure](cm)

    return m_dict[measure](cm)

























