import scipy.io
import numpy as np
import bct
#import ast     #ast.literal_eval
import pprint  #pretty printer
import nibabel as nib
import get_measure as gm



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
def prepare_conn_matrix(mat_matrix):
    
    #extracts the connectivity matrix from the dictionary,
    #'Z' is for Fisher correlation coefficient, 
    #which Conn always(claim?) outputs
    #The last column is grey matter, which needs to
    #be dropped to make the matrix NxN
    conn_cm = mat_matrix['Z']

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
    sym_cm = np.delete(numpy_cm, 32, 1)

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
    pearson_cm = np.absolute(pearson_cm)


    return pearson_cm



'''
Parameters
----------

tm : np.ndarray
     the matrix that needs to be thresholded
th_val : float
         the thresholding value

Returns : np.ndarray
          the matrix with every value
          under the threshold filtered out,
          replaced by 0's    

'''
def threshold_matrix(tm, th_val):
    return np.where(th_val > tm, 0, tm)



'''
Parameters
----------

vec : Nx1 np.ndarray
      any vector result measure obtained from bctpy
size : tuple(int, int, int)
       size/shape for the 3D image to be embedded in NiftI file.
       Default is (91,109,91).
ROI_template : string
               The filename which shall be used as the template for
               the ROI's. Namely, what masks shall be used for composing
               all of the ROI's together in the same matrix with
               values filled in from the vector result
affine : 4x4 np.ndarray
         The affine transformation to be applied when 
         calling nib.nifti1.Nifti1Image().
         Default is the identity matrix.

Returns
-------

img : nibabel.nifti1.Nifti1Image
      The nifti1 file that will be returned.
      The function will not save the nifti1 image, 
      only return it. 

'''

def make_nifti_image(vector, 
                     size=(91,109,91), 
                     affine=None,
                     ROI_template='networks.nii'):
     
    #create the matrix to hold the values
    tmp = np.zeros(size)

    #load the nib file into Python
    nib_template = nib.load(ROI_template)
    #extract the mask matrices
    nib_data = nib_template.get_data()
    #remember that networks.nii was created by MATLAB,
    #hence it is column major(FORTRAN style),
    #we need row major(C style) for NumPy
    masks = np.transpose(nib_data)

    vlen = len(vector)

    #each ROI will encode the obtained value from the vector results
    for i in range(vlen):
        tmp = np.where(masks[i] == 1, vector[i], tmp)

    #create an identity transformation as default
    if affine==None:
        affine = np.eye(4,4)
    
    #convert numpy array to Nifti image,
    #might wanna revisit this/reread documentation
    img = nib.nifti1.Nifti1Image(tmp, affine)

    return img

    


'''
Parameters
----------

img : nibabel.nifti1.Nifti1Image
      The image to be saved to disk
name : string
       The name the image file should have when being saved

Returns
-------

(void) : This function will only save to disk, 
         and print to screen.


Notes: Try and see if a folder can be save to instead just of current directory.
       Should be needed for multiple vector results.

'''


def save_image(img, name):

    try:
        nib.save(img, name)
        print("Succesfully saved to disk: **" + str(name) + "** ")
    except:
        print("Could not save to disk: **" + str(name) + "** ")

    return




'''
Parameters
----------

res_meas : dictionary
           The dictionary containing all the bct measures extracted
           from the connectivity matrix. 

Returns
-------

(void) : save_results will not return any value, instead it will save
         any vector measures as a .nii file, any other value(measures consisting
         of just a single value) will be saved to a .txt file.  

Notes
-----
    save_results uses a nested helper function, save_results_helper, which has the purpose
    to type check the arguments of the dictionary indices. As of now, there are three possible
    types returned by bct: 1) tuples, where each element in the tuple is recursively called, 
    such that they may be saved correctly 2) np.nadarray, which are the vector results obtained
    from bct 3) singular values, like floats and ints, which are all accumulated in a list, 
    before being saved to a .txt file. 
    The function could probably be made more elegant in it's overall structure. 
    The builtin Python function, isinstance(), was chosen, even though it is not very Pythonic,
    becuase it made for simpler code than 'multimethods', and ducktyping was not really helpful
    becuase we would have to extend the types with methods to be able to utilize it properly. 

'''




def save_results(res_meas):

    res_list = []

    ##############################################
    def save_results_helper(meas, key, index):

        res_list_help = []
    
        if (isinstance(meas, tuple)): 
            for i in range(len(meas)):
                save_results_helper(meas[i], key, i)

        elif (isinstance(meas, np.ndarray)):
            img = make_nifti_image(meas)
            save_image(img, key)
        else:
            res_list.append((key, meas, index))
    ###############################################

    for key in res_meas:
        save_results_helper(res_meas[key], key, 0)

    if (res_list == []):
        print('No singular values were found...')
        return

    else: 
        print('Writing to measures.txt...')
        print(res_list)
        with open("measures.txt", 'w') as f:
            for tupl in res_list:
                f.write(str(tupl[0]) + ':' + str(tupl[1]) + ':' + str(tupl[2]) + '\n')
        return







if __name__ == "__main__":
    #hmm, seems like hack (for pesky NaN values)
    # np.seterr(invalid='ignore')
    conn_filename = 'resultsROI_Subject001_Condition001.mat'

    #load the .mat file into a Python dictionary
    mat_matrix = scipy.io.loadmat(conn_filename)
  # print(mat_matrix)
    data = prepare_conn_matrix(mat_matrix)

   # data = prepare_conn_matrix('resultsROI_Subject001_Condition001.mat')
   # mat_matrix = scipy.io.loadmat('resultsROI_Subject001_Condition001.mat')
    # names = mat_matrix['names']
    md = {}
    # md['Eloc_wei'] = bct.efficiency_wei(data, True)
    # roi = make_nifti_image(md['Eloc_wei'])
    md['modularity_und'] = bct.modularity_und(data)[0]
    md['clustering_coef_wu'] = bct.clustering_coef_wu(data)
    md['efficiency_wei'] = bct.efficiency_wei(data)

    #print(type(mat_matrix))
    #print(roi)
    #print(type(roi))
   # nib.viewers.OrthoSlicer3D(roi.dataobj).show()
   # save_image(roi, 'heste.nii')
    pp = pprint.PrettyPrinter(depth=6)
    test = {}
    test_list = ['modularity_und', 'clustering_coef_wu', 'efficiency_wei']

    for key in test_list:
        test[key] = gm.vec_measure(key, data)
    print('test')
    pp.pprint(test)
    print('')
    print('md')
    pp.pprint(md)

    pp.pprint(type(test))
    print('type of modularity_und:')
    pp.pprint(type(test['modularity_und']))
    pp.pprint(type(test['modularity_und'][0]))
    pp.pprint(type(test['modularity_und'][1]))
    pp.pprint(type(test['modularity_und'][0][0]))
    print('type of efficiency_wei:')
    pp.pprint(type(test['efficiency_wei']))
    
    print('isinstance test:')
    print(isinstance(test['modularity_und'], tuple))
    print(isinstance(test['modularity_und'][0], np.ndarray))


    save_results(test)



    # img = nib.load('modularity_und.nii')
    # print(img.header)

    # nib.viewers.OrthoSlicer3D(img.get_data()).show()


        








