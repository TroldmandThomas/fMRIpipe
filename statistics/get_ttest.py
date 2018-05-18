import numpy as np
from scipy import stats
import pandas as pd
import glob
import pprint
from collections import OrderedDict
from operator import itemgetter
pp = pprint.PrettyPrinter(depth=6)



def get_norm_dist(d, alpha, g, s,nor_t='ks'):

    #group the data with regards to HC/SAD and Winter/Summer
    guys = d['groups'].get_group((g, s))
    
    #drop the unuse columns, so we can just iterate over the data structure
    nd = guys.drop(columns=['Unnamed: 0', 'Threshold', 'Group', 'Season'])
    
    #dictionary of accepted and rejected hypothesises
    rad = OrderedDict()
    rd = OrderedDict()
    ad = OrderedDict()

    for item in nd:

            normalized = (nd[item] - nd[item].mean()) / nd[item].std()
            
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

    rad['accepted'] = ad
    rad['rejected'] = rd 
    rad['thresh_percent'] = d['thresh_percent']

    return rad



def compute_ttest(d, hc_rad, sad_rad, alpha, s):

    #group our subjects according to HC, SAD and the season (summer/winter)
    SAD_guys = d['groups'].get_group(('Case', s)).drop(columns=['Unnamed: 0', 'Threshold', 'Group', 'Season'])
    HC_guys = d['groups'].get_group(('Healthy Control', s)).drop(columns=['Unnamed: 0', 'Threshold', 'Group', 'Season'])

    #following code is to check that the attribute is normally distributed in BOTH groups
    both_norm = []
    rejected_norm = []

    for item in hc_rad['accepted']:
        if item in sad_rad['accepted']:
            both_norm.append(item)
        else:
            print('The measure : ' + str(item) + ' was not normally distributed in both lists at ' 
                  + str(hc_rad['thresh_percent'] + '% threshold'))
            rejected_norm.append(item)

    for item in sad_rad['accepted']:
        if item not in hc_rad['accepted']:
             print('The measure : ' + str(item) + ' was not normally distributed in both lists at ' 
                  + str(hc_rad['thresh_percent'] + '% threshold'))
           #  print(str(item) + ' had a p-value of ' + str(sad_rad['accepted'][item][1]))
           #  print(str(item) + ' had a p-value of ' + str(hc_rad['rejected'][item][1]))
             rejected_norm.append(item)
         
    dr = OrderedDict()
    dr['rejected_norm'] = rejected_norm
    dr['thresh_percent'] = d['thresh_percent']

    #USED FOR UNADJUSTED ALPHA VALUE T-TESTING
    rd = OrderedDict()
    ad = OrderedDict()

    #should we Normalize data here also ?? 
    #null hypothesis is data is similar. 
    #when p > alpha, null hypothesis cannot be rejected
    #when p < alpha, null hypothesis is rejected

    #BENJAMIN-HOCHBERG ALGORITHM
    #MIGHT BE BREAKING SINGLE RESPONSIBILITY PRINCIPLE
    #BUT RESPONSIBILITY IS TO DECIDE WHETHER OR NOT T-TESTS CAN BE REJECTED
    #USED FOR MULTIPLE COMPARISON TESTS
    ttest_dict = OrderedDict()
    pval_list = []

    for item in both_norm:

        normalized_SAD = (SAD_guys[item] - SAD_guys[item].mean()) / SAD_guys[item].std()
        normalized_HC = (HC_guys[item] - HC_guys[item].mean()) / HC_guys[item].std()
        ttest_res = stats.ttest_ind(normalized_SAD, normalized_HC, equal_var=False)

        #ttest_res = stats.ttest_ind(SAD_guys[item], HC_guys[item], equal_var=False)
        ttest_dict[item] = ttest_res
        pval_list.append((item, ttest_res[0],ttest_res[1]))
        

        #CODE FOR T-TESTING WITH UNADJUSTED ALPHA VALUE
        #POSSIBLE TO CORRECT WITH CONSERVATIVE BONFERRONI CORRECTION BEFORE CALLING FUNCTION
        # if ttest_res[1] > alpha:
        #     ad[item] = ttest_res
        # else:
        #     print("T-testing of : " +str(item) + \
        #         " with p-value of: " +str(round(ttest_res[1],6)) + \
        #         " was rejected in threshold of " +str(d['thresh_percent']) + "%")
        #     rd[item] = ttest_res

    # print('OG pval')
    # print(pval_list)
    # print('Sorted pval')
    #og_pval_list = pval_list
    pval_list.sort(key=itemgetter(2))
    m = len(pval_list)
    crit_list = []
    for i in range(m):
        crit_val = (i+1)/m * alpha
        crit_list.append(crit_val)
    
    comp = -1
    for j in range(m):
        if crit_list[j] > pval_list[j][2]:
            comp = j



    if comp == -1:
        print('No significant results in threshold of : ' +str(d['thresh_percent']) + '%')
        dr['accepted_ttest'] = pval_list
        dr['rejected_ttest'] = []
    else:
        comp += 1
        ad_list = pval_list[comp:]
        rd_list = pval_list[:comp]
        dr['accepted_ttest'] = ad_list
        dr['rejected_ttest'] = rd_list 
        #print(len(rd_list))
        for k in range(len(rd_list)):
            print('T-testing rejected null hypothesis of: ' \
                +str(rd_list[k][0]) + ' with a p-value of :' +str(round(rd_list[k][2],6)) \
                  + ' in the threshold with ' + str(d['thresh_percent']) + '%')
     
    

    # dr['accepted_ttest'] = ad
    # dr['rejected_ttest'] = rd 

    return dr 




def gtt_main(WS='S',alpha_norm=0.05,alpha_ttest=0.05,nt='ks'):

    #find the estimate csv files
    fl = glob.glob('auto_results/estimates/NETestimate.??.csv')
    
    #create a list of dictionaries, 
    #then load the sorted groups along with the threshold into it
    dfl = []
    thl = []

    for f in fl:
        df = pd.read_csv(f)
        thp = f.split('.')[::-1][1]
        groups = df.groupby(['Group', 'Season'])

        d = OrderedDict()
        d['thresh_percent'] = thp
        d['groups'] = groups

        thl.append(thp)
        dfl.append(d)


 
    # #BONFERRONI CORRECTION INTERMEZZO
    # count = dfl[0]['groups'].get_group(('Case', WS))\
    #         .drop(columns=['Unnamed: 0', 'Threshold', 'Group', 'Season']).count(axis=1)
    # bonferroni = alpha_ttest / count[0]
    # alpha_ttest = bonferroni
    # print("Computing with a Bonferroni corrected alpha-value of: " +str(bonferroni))
    #print(bonferroni)

    
    #compute the t-test for the various thresholds
    ct_list = []
    rad_dict = OrderedDict()
    

    for item in dfl:
        hc_rad = get_norm_dist(item, alpha_norm, 'Healthy Control', WS,nor_t=nt)
        sad_rad = get_norm_dist(item, alpha_norm, 'Case', WS,nor_t=nt)
        ct_list.append(compute_ttest(item, hc_rad, sad_rad, alpha_ttest, WS))
        rad_dict[item['thresh_percent']] = OrderedDict()
        rad_dict[item['thresh_percent']]['HC'] = hc_rad
        rad_dict[item['thresh_percent']]['SAD'] = sad_rad


    return (ct_list, dfl, rad_dict, thl)


    #pp.pprint(ct_list)

#gtt_main(WS='W',alpha_ttest=0.05)













# def draw_graphs(data, ttest, metric, s='S',ymin=0, ymax=1, xmin=10, xmax= 42):
    
    
#     for i in range(len(data)):
#         thresh = data[i]['thresh_percent']
#         SADS = data[i]['groups'].get_group(('Case', s))
#         HCS = data[i]['groups'].get_group(('Healthy Control', s))
#         plt.plot(thresh, SADS[metric].mean(), 'ro')
#         plt.plot(thresh, HCS[metric].mean(), 'bo')
    
#     for j in range(len(ttest)):
#         try:
#             temp = t[j]['rejected_ttest'][metric]
#             plt.plot(t[j]['thresh_percent'], ymax-0.5, 'g*', markersize=12)
#         except:
#             continue
    
#     plt.axis([xmin, xmax, ymin, ymax])
#     plt.xlabel('Sparsity threshold %')
#     plt.ylabel('Global mean ')
#     plt.title(metric)




# def draw_graphs(data, ttest, metric, s='S',ymin=0, ymax=1, xmin=10, xmax= 42):
    
#     for i in range(len(data)):
        
#         thresh = data[i]['thresh_percent']
#         SADS = data[i]['groups'].get_group(('Case', s))
#         HCS = data[i]['groups'].get_group(('Healthy Control', s))
        
#  #       plt.plot(thresh, SADS[metric].mean(), 'ro')
#  #       plt.plot(thresh, HCS[metric].mean(), 'bo')
 
#         y1error = SADS[metric].std()
#         y2error = HCS[metric].std()
#         plt.errorbar(np.array(thresh), SADS[metric].mean(), yerr=y1error, color='red', marker='D')
#         plt.errorbar(np.array(thresh), HCS[metric].mean(), yerr=y2error, color='blue', marker='D')
        
    
#     for j in range(len(ttest)):
#         try:
#             temp = t[j]['rejected_ttest'][metric]
#             plt.plot(t[j]['thresh_percent'], ymax-0.5, 'k*', markersize=12)
#         except:
#             continue
    
#     plt.axis([xmin, xmax, ymin, ymax])
#     plt.xlabel('Sparsity threshold %')
#     plt.ylabel('Global mean ')
#     plt.title(metric)






# draw_graphs(kd,kt,'charpath-lambda', kr, kp,ymin=5,ymax=12.5)

# draw_graphs(kwd,kwt,'charpath-lambda', kwr, kwp,s='W',ymin=5,ymax=12.5)

# draw_graphs(kd,kt,'efficiency_wei-Eglob', kr, kp, ymin=0,ymax=0.3)

# draw_graphs(kwd,kwt,'efficiency_wei-Eglob', kwr, kwp,s='W', ymin=0,ymax=0.3)

# draw_graphs(kd,kt,'assortativity_wei-r',kr,kp,ymin=0,ymax=1.0)

# draw_graphs(kwd,kwt,'assortativity_wei-r',kwr,kwp,s='W',ymin=0,ymax=1.0)

# draw_graphs(kd,kt,'avg_clustering_coef_wu:C',kr,kp,ymin=0,ymax=0.4)

# draw_graphs(kwd,kwt,'avg_clustering_coef_wu:C',kwr,kwp,s='W',ymin=0,ymax=0.4)

# draw_graphs(kd,kt,'modularity_und-Q',kr,kp,ymin=0,ymax=1.0)

# draw_graphs(kwd,kwt,'modularity_und-Q',kwr,kwp,s='W',ymin=0,ymax=1.0)

# draw_graphs(kd,kt,'small_worldness:S',kr,kp,ymin=0,ymax=8.0)

# draw_graphs(kwd,kwt,'small_worldness:S',kwr,kwp,s='W',ymin=0,ymax=8.0)

# draw_graphs(kd,kt,'transitivity_wu-T',kr,kp,ymin=0,ymax=0.5)

# draw_graphs(kwd,kwt,'transitivity_wu-T',kwr,kwp,s='W',ymin=0,ymax=0.5)













