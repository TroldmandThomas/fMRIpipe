import numpy as np
from scipy import stats #ttest_ind, ks_test
import pandas as pd 
import glob #for regex like functionality of finding files
import pprint
import sys
import os
from collections import OrderedDict
from operator import itemgetter
import pathlib #only Python 3.5+

pp = pprint.PrettyPrinter(depth=6)


'''
Parameters
----------

d : OrderedDict,
    A dictionary containing the dataframes along with the threshold percentage
    they were estimated with.

alpha : float,
        the level of significance to be tested against

g : string
    the group to test, e.g. 'Case' for SAD or 'Healthy Control' for HC's

s : string
    the season to test, e.g. 'S' for summer and 'W' for winter

nor_t : string
        whether to use Kolmogorov-Smirnov test ('ks') or Shapiro-Wilks test ('shapiro')
        for testing the samples for normality distribution

Returns
-------

rad : OrderedDict,
      the 'Rejected and Accepted Dictionary', contains smaller dicts with 
      every sample that failed test of normality (in rad['rejected']) and
      every sample which passed the test of normality (in rad['accepted']),
      along with the threshold percentage the data/estimates came from. 

Notes
-----

It is hardcoded that this file runs tests from auto_results/estimate.??.csv', 
i.e. the CSV files that were produced from obtain_estimates.py. 

'''


def get_norm_dist(d, alpha, g, s,nor_t='ks'):

    #group the data with regards to HC/SAD and Winter/Summer
    guys = d['groups'].get_group((g, s))
    
    #drop the unused columns, so we can just iterate over the data structure
    nd = guys.drop(columns=['Unnamed: 0', 'Threshold', 'Group', 'Season'])
    
    #dictionary of accepted and rejected hypothesises
    rad = OrderedDict()
    rd = OrderedDict()
    ad = OrderedDict()

    for item in nd:

            #normalize our data first, to ttest against the implementation of N ~ (0,1)
            normalized = (nd[item] - nd[item].mean()) / nd[item].std()
            
            #pick either ks or shapiro test, can be expanded by other tests
            if nor_t == 'ks':
                norm_test = stats.kstest(normalized,'norm')
            elif nor_t == 'shapiro':
                norm_test = stats.shapiro(normalized)

            #null hypothesis is data is normal distributed. 
            #when p > alpha, null hypothesis cannot be rejected
            #when p < alpha, null hypothesis is rejected
            if norm_test[1] > alpha:
                ad[item] = norm_test
            else:
                print("Rejected normal distribution of : " \
                    + str(item) + " in threshold of : " + str(d['thresh_percent']) + '%' + \
                    ' \n with ' + str(norm_test) )
                rd[item] = norm_test

    #store the rejecte and accepted dictionaries into a single dictionary (the rad)
    rad['accepted'] = ad
    rad['rejected'] = rd 
    rad['thresh_percent'] = d['thresh_percent']

    return rad


'''
Parameters
----------

d : OrderedDict,
    A dictionary containing the dataframes along with the threshold percentage
    they were estimated with.

hc_rad : OrderedDict,
         the 'Rejected and Accepted Dict' for the Healthy Control group

sad_rad : OrderedDict,
          the 'Rejected and Accepted Dict' for the Case group

alpha : float,
        the level of significance to test against

s : string
    the season to test ('S' for summer, 'W' for winter)


Returns
-------

dr : OrderedDict,
     the dictionary of rejections, mainly containing those samples which
     were not normally distributed under both the HC sampls AND the SAD samples.

ttest_csv : list
            data list for the results of the t-testing, will be saved to a csv files
            with raw, unadjusted p-values.

Notes
-----

The return value is not the essential result here. Rather, the function will print
out the hypotheses that were rejected, or report if there were no significant 
p-values in the given threshold. This might be counterintuitive, and could perhaps be made 
as a return value in the future, rather than just printing them out. 

'''



def compute_ttest(d, hc_rad, sad_rad, alpha, s):

    #group our subjects according to HC, SAD and the season (summer/winter)
    SAD_guys = d['groups'].get_group(('Case', s)).drop(columns=['Unnamed: 0', 'Threshold', 'Group', 'Season'])
    HC_guys = d['groups'].get_group(('Healthy Control', s)).drop(columns=['Unnamed: 0', 'Threshold', 'Group', 'Season'])

    #following code is to check that the attribute is normally distributed in BOTH groups
    both_norm = []
    rejected_norm = []

    for item in hc_rad['accepted']:
        if item in sad_rad['accepted']:
            both_norm.append(item)         #both samples were normally distributed
        else:
            #some sample was not normally distributed in the sad_rad 
            print('The measure : ' + str(item) + ' was not normally distributed in both lists at ' 
                  + str(hc_rad['thresh_percent'] + '% threshold'))
            rejected_norm.append(item)

    for item in sad_rad['accepted']:
        if item not in hc_rad['accepted']:
            #some sample was not normally distributed in the hc_rad 
             print('The measure : ' + str(item) + ' was not normally distributed in both lists at ' 
                  + str(hc_rad['thresh_percent'] + '% threshold'))
             rejected_norm.append(item)
    
    #dr is the 'dictionary of rejections (for nomality distributions in both samples)'
    dr = OrderedDict()
    dr['rejected_norm'] = rejected_norm
    dr['thresh_percent'] = d['thresh_percent']

    #USED FOR UNADJUSTED ALPHA VALUE T-TESTING
    rd = OrderedDict()
    ad = OrderedDict()

    #null hypothesis is data is similar. 
    #when p > alpha, null hypothesis cannot be rejected
    #when p < alpha, null hypothesis is rejected

    ttest_dict = OrderedDict()
    ttest_res_temp = []
    pval_list = []

    #those samples which were both normally distributed, 
    #we run a t-test on them. 
    for item in both_norm:

        #DO NOT NORMALIZE DATA BEFORE T-TESTING, MEANS WILL BE IDENTICAL
        # #normalize our data before t-testing
        # normalized_SAD = (SAD_guys[item] - SAD_guys[item].mean()) / SAD_guys[item].std()
        # normalized_HC = (HC_guys[item] - HC_guys[item].mean()) / HC_guys[item].std()

        #use Welch's t-test, i.e. unknown variance and unequal sample sizes
        ttest_res = stats.ttest_ind(SAD_guys[item], HC_guys[item], equal_var=False)

        #store the resulting tuple (t-statistic, p-value)
        ttest_dict[item] = ttest_res

        #for the csv file:
        #ttest_res_temp.append(ttest_res[1])

        #make a list of the metric and the results
        pval_list.append((item, ttest_res[0],ttest_res[1]))

    #data for our ttest csv file
    ttest_lister = pval_list
    ttest_lister.sort(key=itemgetter(0))
    for i in range(len(ttest_lister)):
        ttest_res_temp.append(ttest_lister[i][2])
    ttest_csv = [d['thresh_percent'],s] + ttest_res_temp 

    #BENJAMIN-HOCHBERG ALGORITHM
    #used for FDR correction
    pval_list.sort(key=itemgetter(2))
    
    m = len(pval_list)
    crit_list = []
    for i in range(m):
        crit_val = ((i+1)/m) * alpha
        crit_list.append(crit_val)
    
    comp = -1
    for j in range(m):
        if crit_list[j] > pval_list[j][2]:
            comp = j



    #the above computes the q-value from our p-values, 
    #now we report what hypotheses we can reject/fail-to-reject

    #code for fail-to-reject
    if comp == -1:
        print('No significant results in threshold of : ' +str(d['thresh_percent']) + '%')
        dr['accepted_ttest'] = pval_list
        dr['rejected_ttest'] = []
    #code for appending all those metrics which reject the null hypothesis
    else:
        comp += 0
        ad_list = pval_list[comp:]
        rd_list = pval_list[:comp]
        dr['accepted_ttest'] = ad_list
        dr['rejected_ttest'] = rd_list 
        for k in range(len(rd_list)):
            print('T-testing rejected null hypothesis of: ' \
                +str(rd_list[k][0]) + ' with a p-value of :' +str(round(rd_list[k][2],6)) \
                  + ' in the threshold with ' + str(d['thresh_percent']) + '%')

    return (dr,ttest_csv) 



'''
Parameters
----------

d : OrderedDict,
    A dictionary containing the dataframes along with the threshold percentage
    they were estimated with.

hc_rad : OrderedDict,
         the 'Rejected and Accepted Dict' for the Healthy Control group

sad_rad : OrderedDict,
          the 'Rejected and Accepted Dict' for the Case group

alpha : float,
        the level of significance to test against

s : string
    the season to test ('S' for summer, 'W' for winter)


Returns
-------

res_list : list,
           the list of the results from the Mann-Whitney U-test

'''


def compute_mannwhitney(d, hc_rad, sad_rad, alpha, s):

    #group our subjects according to HC, SAD and the season (summer/winter)
    SAD_guys = d['groups'].get_group(('Case', s)).drop(columns=['Unnamed: 0', 'Threshold', 'Group', 'Season'])
    HC_guys = d['groups'].get_group(('Healthy Control', s)).drop(columns=['Unnamed: 0', 'Threshold', 'Group', 'Season'])
    
    res_list = []


    for item in hc_rad['rejected']:
        res = stats.mannwhitneyu(SAD_guys[item], HC_guys[item], alternative='two-sided')
        res_list.append((res,item))


    for item in sad_rad['rejected']:
        res = stats.mannwhitneyu(SAD_guys[item], HC_guys[item], alternative='two-sided')
        res_list.append((res,item))


    # #BENJAMIN-HOCHBERG ALGORITHM
    # #used for FDR correction
    res_list.sort(key=itemgetter(0))
    res_list = res_list[::-1]
    m = len(res_list)
    crit_list = []
    for i in range(m):
        crit_val = ((i+1)/m) * alpha
        crit_list.append(crit_val)
    comp = -1
    for j in range(m):
        if crit_list[j] >= res_list[j][0][1]:
            comp = j
    # #the above computes the q-value from our p-values, 
    # #now we report what hypotheses we can reject/fail-to-reject

    # #code for fail-to-reject
    if comp == -1:
        print('No significant differences in non-normal distributed samples through ' \
             'Mann-Whitney U-test in threshold of : ' + str(hc_rad['thresh_percent']) + '%')

    #code for appending all those metrics which reject the null hypothesis
    else:
        comp += 1
        ad_list = res_list[comp:]
        rd_list = res_list[:comp]
        for k in range(len(rd_list)):
            print('Mann Whitney U-test found significant differences in: ' \
                +str(rd_list[k][1]) + ' with a p-value of :' +str(round(rd_list[k][0][1],6)) \
                  + ' in the threshold with ' + str(d['thresh_percent']) + '%')


    return res_list
    


'''
Parameters
----------

WS : string,
     'Winter or Summer', 'S' for summer, 'W' for winter,
     which season in question are to be statistically analyzed

alpha_norm : float,
             which level of significance should be used for testing 
             that the samples are normally distributed

alpha_ttest : float,
              which level of significance should be used for testing
              that the samples are significantly different from each other

nt : string
     'Normality Test' to be used. 
     Currently 'ks' for Kolmogorov-Smirnov and 'shapiro' for 
     Shapiro-Wilks test is supported (same as the ones for get_norm()) 

Returns
-------

ct_list : OrderedDict,
          contains the results from the computed t-testings

dfl : list
      a list of dataframes that were used in the testing

rad_dict : OrderedDict,
           the 'Rejected and Accepted Dict' which contains an overview
           of which samples which failed tests for normality, and which 
           that did not

thl_list : list
           a list of which threshold percentages were used for the analysis


Notes
-----           

Currently the file path has to be hardcoded in here, this should be 
fixed soon. 

'''

def gtt_main(WS='S',alpha_norm=0.05,alpha_ttest=0.05,nt='ks', path=None, dest=None):

    #     path = os.path.dirname(os.path.dirname( __file__ ))
    if path == None:
        print(' ')
        print('**Please provide a path to the estimate files**')
        exit()

    # find the estimate files in the given directory
    fl = glob.glob(str(path) + '/estimate.??.csv')
    #'auto_results/estimates/AALestimate.??.csv'
    
    #create a list of dictionaries, 
    #then load the sorted groups along with the threshold into it
    dfl = []     #pandas dataframe list 
    thl = []     #threshold percentage list

    for f in fl:
        df = pd.read_csv(f)
        thp = f.split('.')[::-1][1]               #get the threshold percentage
        groups = df.groupby(['Group', 'Season'])  

        d = OrderedDict()
        d['thresh_percent'] = thp
        d['groups'] = groups

        thl.append(thp)
        dfl.append(d)

    
    #compute the t-test for the various thresholds
    ct_list = []
    rad_dict = OrderedDict()
    t_csv = []
    frames = []
    

    for item in dfl:
        #perform the tests for normality distribution
        hc_rad = get_norm_dist(item, alpha_norm, 'Healthy Control', WS,nor_t=nt)
        sad_rad = get_norm_dist(item, alpha_norm, 'Case', WS,nor_t=nt)

        #perform the actual t-testing
        ttest_result,ttest_csv = compute_ttest(item, hc_rad, sad_rad, alpha_ttest, WS)
        

        #save the results in dictionaries for return value
        ct_list.append(ttest_result)
        rad_dict[item['thresh_percent']] = OrderedDict()
        rad_dict[item['thresh_percent']]['HC'] = hc_rad
        rad_dict[item['thresh_percent']]['SAD'] = sad_rad

        #build up the dataframes to save as a CSV file
        #CSV for normality tests
        frames.append(pd.DataFrame(hc_rad))
        frames.append(pd.DataFrame(sad_rad))
        #CSV for ttests
        t_csv.append(ttest_csv)

        #perform Wilcoxon rank sum if there are any rejections of normal distributions
        if len(hc_rad['rejected'])!= 0 or len(sad_rad['rejected'])!=0:
            compute_mannwhitney(item, hc_rad, sad_rad, alpha_ttest, WS)
    

    if dest != None:    
        #save the results to CSV files
        #first make a directory to save the files to
        dest = dest + '/tests'
        pathlib.Path(dest).mkdir(parents=True, exist_ok=True)
    
        norms = pd.concat(frames)

        ttests = pd.DataFrame(t_csv, columns=['Threshold', 'Season', \
            'Assortativity', 'ClusteringCoefficient', 'CharPath', 'GlobalEfficiency', \
            'Modularity', 'SmallWorld', 'Transitivity'])

        #paths to the CSV files
        #save the normality tests to a CSV file
        norms.to_csv(dest + '/' + WS + '_normality.csv')
        #save the ttests results to a CSV file
        ttests.to_csv(dest + '/'+ WS + '_ttests.csv')


    return (ct_list, dfl, rad_dict, thl)

#if used as an individual file
if __name__ == "__main__":
    
    #get the season and the alpha level from command line args
    try: 
        season = sys.argv[1]
        if len(sys.argv) == 2:
            alpha = sys.argv[2]
        else:
            alpha = 0.05
    except:
        print('Usage: python3.6 (season=\'W\'/\'S\') (alpha)')
        exit()

    #run the script
    ret = gtt_main(WS=season,alpha_ttest=alpha, nt='ks')

    













