import numpy as np 
import pathlib #only Python 3.5+ 
import pprint  #pretty printer
import nibabel as nib  #for saving .nii nifti images
from collections import OrderedDict
import pandas as pd #for csv file creation
import utils.progressbar as pb #progressbar courtesy of stackoverflow
import sys #for getting commandline arguments
import pipeline.loadmatrix as lm #getting the connectivity matrices from Conn
import pipeline.graph_estimates as ge


##########################################################################
#GLOBAL VARIABLES
#at the moment, just some path names to save results, 
#and for reading the mask for ROI's

ROI_template = 'networks.nii'
mask_name = ROI_template.split('/')[-1].split('.')[0]
# result_dir = 'auto_results'
# vector_dir = '/vector_measures/'
# csv_dest = '/estimate.csv'
# dest = result_dir + vector_dir


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
dest : string
       the path to the directory where the resulting vector files will be put


Returns
-------

(void) : This function will only save to disk, 
         and print to screen.


Notes: Try and see if a folder can be save to instead just of current directory.
       Should be needed for multiple vector results.

'''


def save_image(img, mn, sn, dest):

    #make a directory to contain our vector estimate files,
    #does nothing if the directory already exists.
    #Python 3+ dependent
    pathlib.Path(dest).mkdir(parents=True, exist_ok=True) 
    
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

    #function for checking if we have acquired a local measure
    # or a global measure. Local measures will be saved as 
    #a nifti image file. 
    def filter_helper(est, key):

        #case for local measure
        if (isinstance(est, np.ndarray)):
            #REINCLUDE COMMENTED CODE; JUST SAVES TIME FOR CSV ESTIMATES
            # img = make_nifti_image(est)
            # save_image(img, key, sn)
            key_list.append(key)

        #case for global measure
        elif (isinstance(est, np.float64)):
            #restrict some of all these decimals
            est_dic[key] = np.around(est, decimals=5)
            
    #run the nested function to filter out local measures
    for key in est_dic:
        filter_helper(est_dic[key], key)

    #delete the indices for local measures in our given dictionary,
    #such that we can save the global measures (which are only a single number)
    # to a CSV file via a pandas dataframe 
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
dest : string
       the path to the directory where the resulting estimate files will be put

Returns
-------

(void) : No results are returned, as this is really the "main()"
         function of the program, tying it all together
         #TODO: refine comments


'''

def obtain_estimates(cm_list, groupIDcsv, th, path):

    dic_list = []
    
    #the CSV file used to identify and label the subjects in our matrix file
    iddf = pd.read_csv(groupIDcsv)
    
    #filter out the singular values to be saved to .csv file,
    #the vector measures will be saved to disk
    #For each connection matrix file, 
    #run the graph estimates and store them properly
    th_p = int(th * 100)
    #counter used for the progressbar
    i = 0 
    l = len(cm_list)
    #initialize progressbar so we have a feeling of the time consumed by the script
    pb.printProgressBar(0, l, prefix = 'Progress:', suffix = 'Complete', length = 50)
    for cm in cm_list:

            #perform the actual graph theory estimations
            dic = ge.graph_estimates(cm,th)

            subject_name = str(i)

            #filter out the local measures that won't fit in a CSV file
            filt_dic = filter_singular_values(dic, subject_name)
            
            #note the threshold percentage that the estimates were performed under
            filt_dic['Threshold'] = th_p
            #note the group and season for the subject
            filt_dic['Group'] = iddf['group'][i]
            filt_dic['Season'] = iddf['season'][i]

            dic_list.append(filt_dic)
            pb.printProgressBar(i + 1, l, prefix = 'Progress:', suffix = 'Complete', length = 50)
            #increment counter for progressbar
            i = i + 1 

    #store the singular values in Pandas dataframe,
    #for convienient conversion to .csv file
    df = pd.DataFrame(dic_list)

    #insert the threshold into the csv name
    csv_dest = '/estimate.csv'
    csv_split = csv_dest.split('.')
    csv_name = csv_split[0] + '.' + str(th_p) + '.' +  csv_split[1]

    result_dir = '/auto_results'
    est_dir = path + result_dir
    pathlib.Path(est_dir).mkdir(parents=True, exist_ok=True)

    #save the estimates to our CSV file
    df.to_csv(est_dir + csv_name)
    
    return




#usage for when the script is called by itself
if __name__ == "__main__":

    #checking if the user gave the correct inputs
    #case for when size of given matrix is given as input
    if len(sys.argv) == 5:
        size = sys.argv[4]
    #case for when the number of arguments given is wrong (i.e. not 4 or 5)
    elif len(sys.argv) != 4:
        print('Usage: (python3.6) obtain_estimates.py (resultsROI.mat) (label_csv) (threshold_list) [size]')
        exit()
    #case for when number of arguments is exactly 4, then just use the whole matrices
    else:
        size = 'full'

    #converting file names to a Python list of file names
    cms = list(map(str, sys.argv[1].strip('[]').split(',')))
    #extracting the CSV file for subject labelling
    label_csv = sys.argv[2]
    #converting the threshold list from args to floats
    
    #case for when interval is given
    try:
        thresh_tok = sys.argv[3].split(':')
        start = int(thresh_tok[0])
        end = int(thresh_tok[1])
        stride = int(thresh_tok[2])
        thresh_list = []
        for i in range(start,(end+stride),stride):
            thresh_list.append(float(i)/100)

    #case for when list is given
    except:
        thrs = list(map(str, sys.argv[3].strip('[]').split(',')))
        thresh_list = list(map(float, thrs))

    #convert the matrices from MATLAB to Python
    pm = lm.conn_interface(cms, size)

    #perform the actual graph theory estimates on the given matrices
    for thresh in thresh_list:
        print('Now processing threshold: ' +str(round(100 * thresh,2)) + '%')
        obtain_estimates(pm, label_csv,thresh)
    #done
    print("Graph theory estimates completed on all thresholds")


    
    


        








