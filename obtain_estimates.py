import scipy.io
import numpy as np
import bct
import pathlib #only Python 3.5+ 
import pprint  #pretty printer
import nibabel as nib
from collections import OrderedDict
import graph_estimates as ge
import pandas as pd
import progressbar as pb
import ast #passing lists from command line args

#ONLY FOR TEST, DISABLE WHEN DONE!
import sys
import loadmatrix as lm



##########################################################################
#GLOBAL VARIABLES
#at the moment, just some path names to save results, 
#and for reading the mask for ROI's

result_dir = 'auto_results'
vector_dir = '/vector_measures/'
csv_dest = '/estimate.csv'
ROI_template = 'networks.nii'
mask_name = ROI_template.split('/')[-1].split('.')[0]
dest = result_dir + vector_dir


##########################################################################



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
                     mask_template=ROI_template):
     
    #create the matrix to hold the values
    tmp = np.zeros(size)

    #load the nib file into Python
    nib_template = nib.load(mask_template)
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
mn : string
       The name of the measure to be saved, see below for format
sn : string
     The subject name, for saving in the following format:
     measuretype_mask_subject.nii


Returns
-------

(void) : This function will only save to disk, 
         and print to screen.


Notes: Try and see if a folder can be save to instead just of current directory.
       Should be needed for multiple vector results.

'''


def save_image(img, mn, sn):
    
    formatted_name = mn + '_' + mask_name + '_' + sn
    f = dest + '/' + formatted_name
    
    try:
        nib.save(img, f)
        print("Succesfully saved to disk: **" + str(formatted_name) + "** ")
    except:
        print("Could not save to disk: **" + str(formatted_name) + "** ")

    return



'''
Parameters
----------

est_dic : OrderedDict
          the ordered dictionary which contains the graph estimates
          obtained from bct

sn : string
     Subject name, the name of the subject being processed.
     This is is just for passing along to save_image

Returns
-------

filt_dic : OrderedDict
           the ordred input dictionary, which no longer
           contains any vector estimates, since they are filtered out
           and stored as a .nii file
'''

def filter_singular_values(est_dic, sn):

    key_list = []

    def filter_helper(est, key):

        if (isinstance(est, np.ndarray)):
            #REINCLUDE COMMENTED CODE; JUST SAVES TIME FOR CSV ESTIMATES
            # img = make_nifti_image(est)
            # save_image(img, key, sn)
            key_list.append(key)

        #restrict some of all these decimals
        elif (isinstance(est, np.float64)):
            est_dic[key] = np.around(est, decimals=5)
            

        # elif (isinstance(meas, tuple)): 
        #     for i in range(len(meas)):
        #         save_results_helper(meas[i], key, i)
    for key in est_dic:
        filter_helper(est_dic[key], key)

    for key in key_list:
        del est_dic[key]
    
    return est_dic





'''
Parameters
----------

file_name_list : list
                 A list of all the filenames which should be estimated graph
                 measures upon. 
                 Right now, sample input could be: 
                 [resultsROI_Subject001_Condition001.mat,
                  resultsROI_Subject001_Condition001.mat]

groupIDcsv : csv file
             The accompying csv file to generate the ID tags for each 
             subject in the scan file. 

Returns
-------

(void) : No results are returned, as this is really the "main()"
         function of the program, tying it all together
         #TODO: refine comments


'''

def obtain_estimates(cm_list, groupIDcsv, th):

    dic_list = []
    iddf = pd.read_csv(groupIDcsv)
    #'~/Desktop/Samlemappen/Bachelorprojekt/NewBCT/auto_results/estimate.csv'
    
    #'auto_results/vector_measures'
    pathlib.Path(dest).mkdir(parents=True, exist_ok=True) 
    #pathlib.Path('/auto_results/singular_csv').mkdir(parents=True, exist_ok=True) 
 
    #filter out the singular values to be saved to .csv file,
    #the vector measures will be saved to disk
    #For each connection matrix file, 
    #run the graph estimates and store them properly
    th_p = int(th * 100)
    i = 0 
    l = len(cm_list)
    pb.printProgressBar(0, l, prefix = 'Progress:', suffix = 'Complete', length = 50)
    for cm in cm_list:

            dic = ge.graph_estimates(cm,th)

            subject_name = str(i)
            filt_dic = filter_singular_values(dic, subject_name)
            #filt_dic['ID'] = i
            filt_dic['Threshold'] = th_p
            filt_dic['Group'] = iddf['group'][i]
            filt_dic['Season'] = iddf['season'][i]
            dic_list.append(filt_dic)
            pb.printProgressBar(i + 1, l, prefix = 'Progress:', suffix = 'Complete', length = 50)
            i = i + 1 # TESTING ONLY, DISABLE

    #store the singular values in Pandas dataframe,
    #for convienient conversion to .csv file
    df = pd.DataFrame(dic_list)

    #insert the threshold into the csv name
    csv_split = csv_dest.split('.')
    csv_name = csv_split[0] + '.' + str(th_p) + '.' +  csv_split[1]

    df.to_csv(result_dir + csv_name)

    pp = pprint.PrettyPrinter(depth=6)
    pp.pprint(dic_list)

    
    return





if __name__ == "__main__":

    #converting file names to a Python list of file names
    cms = list(map(str, sys.argv[1].strip('[]').split(',')))
    thresh = float(sys.argv[2])
    print(thresh)
    print(type(thresh))

    pm = lm.conn_interface(cms)
    group_csv = 'groupID_Thomas.csv'
   # print(pm)
   # print(len(pm))
    obtain_estimates(pm, group_csv,thresh)

    #print(cms)

    pp = pprint.PrettyPrinter(depth=6)

    
    


        








