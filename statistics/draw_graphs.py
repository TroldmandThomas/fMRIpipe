import numpy as np
from scipy import stats
import pandas as pd
import pathlib 
import matplotlib.pyplot as plt
import glob
from collections import OrderedDict
import statistics.get_ttest as gtt

#where to put the graphs
dest = 'graphs'

#makes a directory if it does not already exist
#does nothing if the directory exists
#Python3+ dependent
pathlib.Path(dest).mkdir(parents=True, exist_ok=True)

'''
Parameters
----------

data : list,
       a list of dataframes used in the ttesting

ttest : OrderedDict,
        the results of the ttesting, used to display those
        where significant differences were found

metric : string,
         the name of the graph theory metric that are to be plotted

rad : OrderedDict,
      the 'Rejected and Accepted Dict' which contains an overview
      of which samples which failed tests for normality, and which 
      that did not

thrs : list
        a list of which threshold percentages were used for the analysis

s : string,
    the season to plot, can be 'S' or 'W' (summer or winter)


Returns
-------

(void) : does not return anything

Notes
-----

Code tries to figure out the best axis for y and x axis
by taking the min and max of the metrics for the y-axis,
and the min and max for the threshold percentages.


'''

def draw_graphs(data, ttest, metric, rad, thrs, s='S'):

    #initialize the plot axis variables
    ymax = -100
    ymin = 100
    xmax = -1
    xmin = 101

    #clear the plot for good measure
    plt.clf()
    #make the font size on the graph labels a bit bigger
    plt.rcParams.update({'font.size': 20})
    
    #loop over the data
    for i in range(len(data)):
        
        thresh = data[i]['thresh_percent']
        SADS = data[i]['groups'].get_group(('Case', s))
        HCS = data[i]['groups'].get_group(('Healthy Control', s))
 
        #errorbars for our two sample groups
        y1error = SADS[metric].std()
        y2error = HCS[metric].std()
        plt.errorbar(float(thresh), SADS[metric].mean(), yerr=y1error, color='red', marker='D')
        plt.errorbar(float(thresh), HCS[metric].mean(), yerr=y2error, color='blue', marker='D')

        #figure out the best axis values for both our x and y axis
        hc_max = HCS[metric].max()
        case_max = SADS[metric].max()
        hc_min = HCS[metric].min()
        case_min = SADS[metric].min()

        temp_max = max(hc_max,case_max)
        temp_min = min(hc_min,case_min)
        temp_thresh = int(thresh)

        #find the minimum and maximum for our x and y axis
        if ymax < temp_max:
            ymax = temp_max
        if ymin > temp_min:
            ymin = temp_min
        if temp_thresh > xmax:
            xmax = temp_thresh
        if temp_thresh < xmin:
            xmin = temp_thresh

    #adjust the x-axis slightly
    xmin = xmin - 2
    xmax = xmax + 2
        
    
    #code for putting markers on the graphs to show which samples that failed the ttest
    for j in range(len(ttest)):
        
        temp2 = ttest[j]['rejected_norm']
        if metric in temp2:
             plt.plot(ttest[j]['thresh_percent'], ymax-(0.05*ymax), color='yellow', marker='^', markersize=12)
        try:
            temp = ttest[j]['rejected_ttest'][metric]
            plt.plot(ttest[j]['thresh_percent'], ymax-(0.05*ymax), 'k*', markersize=12)
        except:
            pass
        
  #COMMENT IN FOR SHOWING WHICH OF THE NORMAL DISTRIBUTION T-TEST THAT FAILS     
    # for k in range(len(thrs)):
       
    #    try:
    #        temp3 = rad[thrs[k]]['HC']['rejected'][metric]
    #        #plt.plot(float(thrs[k]), ymin+(0.95*ymin), marker='p', color='black', markersize=12)
    #        plt.plot(float(thrs[k]), ymax-(0.05*ymax), marker='p', color='green', markersize=12)
          
    #    except:
    #        pass
       
    #    try:
    #        temp4 = rad[thrs[k]]['SAD']['rejected'][metric]
    #       # plt.plot(float(thrs[k]), ymin+(0.05*ymin), marker='p', color='black', markersize=12)
    #        plt.plot(float(thrs[k]), ymax-(0.05*ymax), marker='p', color='green', markersize=12)
    #    except:
    #        pass

    if s == 'S':
        season = 'Summer'
    else:
        season = 'Winter'
    
    #dictionary for what to label on the plot title
    d = OrderedDict()

    d['assortativity_wei-r'] = 'Assortativity'
    d['avg_clustering_coef_wu:C'] = 'Clustering coefficient'
    d['charpath-lambda'] = 'Characteristic pathlength'
    d['efficiency_wei-Eglob'] = 'Global Efficiency'
    d['modularity_und-Q'] = 'Modularity'
    d['small_worldness:S'] = 'Small Worldness'
    d['transitivity_wu-T'] = 'Transitivity'

    #label our plot axis
    plt.axis([xmin, xmax, ymin, ymax])
    plt.xlabel('Sparsity threshold %')
    plt.ylabel('Global mean ')
    plt.title(d[metric])

    met = metric.split(':')[0]

    #build the name which will be used as the stored file name
    filename = str(dest) + '/' + str(s) + '_' + str(met) + '.png'

    #save the file to disk
    plt.savefig(filename, format='png', bbox_inches='tight')


def execute():
    #extract the summer season data,
    #run the t-tests on the data before and use the return
    #value of get_ttest.py to draw graphs upon
    kwdata = gtt.gtt_main(WS='W',nt='ks')
    kwt = kwdata[0]
    kwd = kwdata[1]
    kwr = kwdata[2]
    kwp = kwdata[3]


    kdata = gtt.gtt_main(WS='S',nt='ks')
    kt = kdata[0]
    kd = kdata[1]
    kr = kdata[2]
    kp = kdata[3]


    #draw the actual graphs
    draw_graphs(kd,kt,'charpath-lambda', kr, kp)

    draw_graphs(kwd,kwt,'charpath-lambda', kwr, kwp,s='W')

    draw_graphs(kd,kt,'efficiency_wei-Eglob', kr, kp)

    draw_graphs(kwd,kwt,'efficiency_wei-Eglob', kwr, kwp,s='W')

    draw_graphs(kd,kt,'assortativity_wei-r',kr,kp)

    draw_graphs(kwd,kwt,'assortativity_wei-r',kwr,kwp,s='W')

    draw_graphs(kd,kt,'avg_clustering_coef_wu:C',kr,kp)

    draw_graphs(kwd,kwt,'avg_clustering_coef_wu:C',kwr,kwp,s='W')

    draw_graphs(kd,kt,'modularity_und-Q',kr,kp)

    draw_graphs(kwd,kwt,'modularity_und-Q',kwr,kwp,s='W')

    draw_graphs(kd,kt,'small_worldness:S',kr,kp)

    draw_graphs(kwd,kwt,'small_worldness:S',kwr,kwp,s='W')

    draw_graphs(kd,kt,'transitivity_wu-T',kr,kp)

    draw_graphs(kwd,kwt,'transitivity_wu-T',kwr,kwp,s='W')



#for use independent of other files
if __name__ == "__main__":
    execute()




