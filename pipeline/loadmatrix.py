import scipy.io 
import numpy as np 
import sys



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

def conn_interface(file_list):

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

            #check whether a single matrix or multiple matrices
            #was given as user input
            if fm_len == 3:
                print('Found Alot of matrices!!')
                for i in range(fmt.shape[0]):
                    temp = prepare_conn_matrix(fmt[i])
                    prepared_matrices.append(temp)
            
            
            elif fm_len == 2:
                print('Found a single matrix!')
                prepared_matrices.append(prepare_conn_matrix(fmt))

            #case for user input is a .mat file, 
            #but it has either 1D or >3D. Just close program for now
            else:
                print('The file is some unknow collection of matrices')
                exit()
        except:
            print("Either the file was not found, or wrong usage of program.")
            print("Usage: python3.6 filename.mat ROI_template ")
            print("Alternative: python3.6 [filename1.mat, filename2.mat] ROI_template")
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
'''
def prepare_conn_matrix(conn_cm):
    
    #extracts the connectivity matrix from the dictionary,
    #'Z' is for Fisher correlation coefficient, 
    #which Conn always(claim?) outputs
    #The last column is grey matter, which needs to
    #be dropped to make the matrix NxN

    #conn_cm = mat_matrix['Z']

    #convert the connectivity matrix to a 2D numpy array
    #(remember numpy matrices(data type) is another thing, more like Linear Algebra)
    #https://docs.scipy.org/doc/numpy-dev/user/numpy-for-matlab-users.html
    #we store the connecitivity matrix in 'FORTRAN' order,
    #which corresponds to MATLAB's column major order.
    #should not make a difference however, since the connectivity
    #matrix is symmetric(for undirected matrices at least)

    numpy_cm = np.array(conn_cm, order='F', dtype=float)

    #drop the last column to make it symmetric and get rid
    #of the gray matter column
    #TEST and see if correct column is dropped

    sym_cm = np.delete(numpy_cm, 32, 0)

    #Conn places NaN in the reflexive connectivity of nodes
    #which screws up Numpy's computation. 
    #we replace the NaN's with zeroes for now,
    #but should maybe be replaced with ones instead.
    #np.nan_to_num only replaces with zeroes however

    no_NaN = np.nan_to_num(sym_cm)

    #sym_cm[np.isnan(sym_cm)] = 6.667
    #no_NaN = sym_cm   #np.absolute(sym_cm)
 
    #perform the inverse Fisher transformation on 'Z' to obtain
    #the Pearson correlation coefficients 'r' in the matrix
    #https://en.wikipedia.org/wiki/Fisher_transformation
    #FIND some OTHER source than wikipedia

    #pearson_cm = no_NaN
    pearson_cm = np.tanh(no_NaN)

    #fill diagonal with 1's, since each node is correlated 100%
    #with itself.
    #np.fill_diagonal(pearson_cm, 1)

    #"The network matrices should not contain self-self connections. 
    #In other words, all values on the main diagonal 
    #of these matrices should be set to 0." 
    #From: https://sites.google.com/site/bctnet/Home/help

    #CLEAN UP COMMENTS
    #compute absolute value elementwise
    #Probbaly dont use absolute values
    #pearson_cm = np.absolute(pearson_cm)
   # print(np.amin(pearson_cm))


    return pearson_cm

   # return conn_cm




if __name__ == "__main__":
    head, *tail = sys.argv
    pm = conn_interface(tail)
    print(len(pm))
    # print(pm[0] == pm[84])
    # print(np.amin(pm[0]))
    # print(np.amin(pm[84]))
    # print(np.amax(pm[0]))
    # print(np.amax(pm[84]))
    













