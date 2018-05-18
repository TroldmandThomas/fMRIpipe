import numpy as np
from scipy import stats
import pandas as pd
import pathlib 
import matplotlib.pyplot as plt
import glob
from collections import OrderedDict
import get_ttest as gtt


dest = 'graphs'

pathlib.Path(dest).mkdir(parents=True, exist_ok=True)




def draw_graphs(data, ttest, metric, rad, thrs, s='S',ymin=0, ymax=1, xmin=8, xmax= 42):

    plt.clf()
    plt.rcParams.update({'font.size': 20})
    
    for i in range(len(data)):
        
        thresh = data[i]['thresh_percent']
        SADS = data[i]['groups'].get_group(('Case', s))
        HCS = data[i]['groups'].get_group(('Healthy Control', s))
 
        y1error = SADS[metric].std()
        y2error = HCS[metric].std()
        plt.errorbar(float(thresh), SADS[metric].mean(), yerr=y1error, color='red', marker='D')
        plt.errorbar(float(thresh), HCS[metric].mean(), yerr=y2error, color='blue', marker='D')
        
    
    for j in range(len(ttest)):
        
        temp2 = ttest[j]['rejected_norm']
        if metric in temp2:
             plt.plot(ttest[j]['thresh_percent'], ymax-(0.05*ymax), color='yellow', marker='^', markersize=12)
            # pass
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

    d = OrderedDict()

    d['assortativity_wei-r'] = 'Assortativity'
    d['avg_clustering_coef_wu:C'] = 'Clustering coefficient'
    d['charpath-lambda'] = 'Characteristic pathlength'
    d['efficiency_wei-Eglob'] = 'Global Efficiency'
    d['modularity_und-Q'] = 'Modularity'
    d['small_worldness:S'] = 'Small Worldness'
    d['transitivity_wu-T'] = 'Transitivity'





    plt.axis([xmin, xmax, ymin, ymax])
    plt.xlabel('Sparsity threshold %')
    plt.ylabel('Global mean ')
    plt.title(d[metric])

    met = metric.split(':')[0]

    filename = str(dest) + '/' + str(s) + '_' + str(met) + '.png'

    

    plt.savefig(filename, format='png', bbox_inches='tight')







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




draw_graphs(kd,kt,'charpath-lambda', kr, kp,ymin=5,ymax=14.5)

draw_graphs(kwd,kwt,'charpath-lambda', kwr, kwp,s='W',ymin=5,ymax=12.5)

draw_graphs(kd,kt,'efficiency_wei-Eglob', kr, kp, ymin=0,ymax=0.3)

draw_graphs(kwd,kwt,'efficiency_wei-Eglob', kwr, kwp,s='W', ymin=0,ymax=0.3)

draw_graphs(kd,kt,'assortativity_wei-r',kr,kp,ymin=0,ymax=0.7)

draw_graphs(kwd,kwt,'assortativity_wei-r',kwr,kwp,s='W',ymin=0,ymax=0.7)

draw_graphs(kd,kt,'avg_clustering_coef_wu:C',kr,kp,ymin=0,ymax=0.4)

draw_graphs(kwd,kwt,'avg_clustering_coef_wu:C',kwr,kwp,s='W',ymin=0,ymax=0.4)

draw_graphs(kd,kt,'modularity_und-Q',kr,kp,ymin=0,ymax=1.0)

draw_graphs(kwd,kwt,'modularity_und-Q',kwr,kwp,s='W',ymin=0,ymax=1.0)

draw_graphs(kd,kt,'small_worldness:S',kr,kp,ymin=0,ymax=8.0)

draw_graphs(kwd,kwt,'small_worldness:S',kwr,kwp,s='W',ymin=0,ymax=8.0)

draw_graphs(kd,kt,'transitivity_wu-T',kr,kp,ymin=0,ymax=0.4)

draw_graphs(kwd,kwt,'transitivity_wu-T',kwr,kwp,s='W',ymin=0,ymax=0.4)


