import scipy.io #scipy.io.loadmat
import numpy as np 
import sys #commnad line arguments
#import bct #thresholding negative weights (MOVED TO graph_estimates)



'''
Parameters:
-----------

file_list : string list
            This is the tail of sys.argv, i.e. a list containing 
            all of the files a user may want to process through the 
            pipeline.


Returns:
--------

prepared_matrices : list of NxN numpy array
                    A list of connectivity matrices will be returned,
                    each extracted from the given files as specified by the user.

Notes:
------

Currently, the function only supports .mat files, which can either be of the form
'resultsROI_Subject001_Condition001.mat' (for a single 2D connectivity matrix),
or of the form 'resultsROI_Condition001.mat' (for a collection of 2D connectivity
matrices). 
The functions is meant to be extended in the future for use by connectivity matrices
obtained through other means than MATLAB Conn. 


'''

def conn_interface(file_list, size='full'):

    prepared_matrices = []

    for f in file_list:
        #check if the given file is a .mat file
        token = f.split('.')[-1]
        if token != 'mat':
            print(str(f) + ' was not a .mat file, closing..')
            exit()

        try:
            #load the matrix given by the .mat fle
            fm = scipy.io.loadmat(f)
            #transpose the matrix to give row order for 3D matrices,
            #2D matrices will also be tranposed, but treated the
            #same as 3D matrices
            fmt = np.transpose(fm['Z'])
            fm_len = len(fm['Z'].shape)

            #check if user gave any other dimensions than just the full matrix,
            #modify the matlab matrix appropriately
            #WILL ASSUME 1-INDEXING IS USED
            if size != 'full':

                dim_tok = size.split('x')
                ns = int(dim_tok[0].split(':')[0]) - 1
                ne = int(dim_tok[0].split(':')[1])       #include the greymatter column, will be dropped later
                                                         #for compatibility with full mode
                ms = int(dim_tok[1].split(':')[0]) - 1
                me = int(dim_tok[1].split(':')[1]) - 1

                #extract the necessary part of the matrices
                fmt = fmt[:][:,ns:ne,ms:me]

            #check whether a single matrix or multiple matrices
            #was given as user input
            if fm_len == 3:
                print('Found multiple matrices in given MATLAB file')
                for i in range(fmt.shape[0]):
                    #prepare the individual matrices
                    temp = prepare_conn_matrix(fmt[i])
                    prepared_matrices.append(temp)
            
            
            elif fm_len == 2:
                print('Found a single matrix in given MATLAB file')
                prepared_matrices.append(prepare_conn_matrix(fmt))

            #case for user input is a .mat file, 
            #but it has either 1D or >3D. Just close program for now
            else:
                print('The file is some unknown collection of matrices')
                exit()
        except:
            print("Unexpected error occured, closing.")
            exit()


    return prepared_matrices


'''
Parameters
----------

mat_matrix : dict
             The dictionary that contains all the data from a .mat file,
             only the Fisher coefficient 'Z' matrix will be used.
             It is ok to pass the whole dict, since only a pointer to object
             is passed (call-by-object Python style)

Returns
-------

pearson_cm : A NxN numpy array,
             which is the connectivity matrix
             expressed by the Pearson correlation coefficient.
             Will also be column major, as the Conn matrix is. 
             NaN's will be overwritten with zeroes, 
             NaN's messes up numpy computations in bctpy.


Notes
-----

Extracts the connectivity matrix from the dictionary,
'Z' is for Fisher correlation coefficient, 
which Conn always outputs.
The last column is grey matter, which needs to
be dropped to make the matrix NxN, i.e. symmetric.


'''
def prepare_conn_matrix(conn_cm):
    
    #convert the connectivity matrix to a 2D numpy array
    #(remember numpy matrices(data type) is another thing, more like Linear Algebra)
    #https://docs.scipy.org/doc/numpy-dev/user/numpy-for-matlab-users.html
    #we store the connecitivity matrix in 'FORTRAN' order,
    #which corresponds to MATLAB's column major order.
    #should not make a difference however, since the connectivity
    #matrix is symmetric(for undirected matrices at least)

    #just to make sure, we let numpy open the matrix in the most 
    #memory efficient way as possible.
    #This is done in O(1) by maniuplating numpy metadata, so should'nt impact run time.
    #https://stackoverflow.com/questions/39264196/read-mat-file-in-python-but-the-shape-of-the-data-changed
    numpy_cm = np.array(conn_cm, order='A', dtype=float)

    #drop the last column to make it symmetric and get rid
    #of the gray matter column
 

    sym_cm = np.delete(numpy_cm, -1, 0)

    #Conn places NaN in the reflexive connectivity of nodes
    #which screws up Numpy's computation. 
    #we replace the NaN's with zeroes,

    no_NaN = np.nan_to_num(sym_cm)
 
    #perform the inverse Fisher transformation on 'Z' to obtain
    #the Pearson correlation coefficients 'r' in the matrix
    #https://en.wikipedia.org/wiki/Fisher_transformation

    pearson_cm = np.tanh(no_NaN)

    #THRESHOLDING OF NEGATIVE WEIGHTS MOVED TO graph_estimate.py,
    #FOR COMPATIBILITY WITH OTHER MATRICES THAN FROM MATLAB
    #threshhold the matrix to remove negative weights
    #thresh_cm = bct.threshold_absolute(pearson_cm, 0.0)

    #"The network matrices should not contain self-self connections. 
    #In other words, all values on the main diagonal 
    #of these matrices should be set to 0." 
    #From: https://sites.google.com/site/bctnet/Home/help

    return pearson_cm



if __name__ == "__main__":
    head, *tail = sys.argv
    pm = conn_interface(tail)
    print(len(pm))
    













