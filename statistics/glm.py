
#https://stackoverflow.com/questions/43043519/is-it-possible-to-do-glmm-in-python
#https://stats.idre.ucla.edu/other/mult-pkg/introduction-to-generalized-linear-mixed-models/

#use these in terminal on Mac(maybe UNIX also) if getting alot of wierd warnings from rpy2
# export LANG=en_US.UTF-8
# export LC_ALL=en_US.UTF-8
# locale

#packages used: pandas, glob and rpy2

import matplotlib.pyplot as plt
import pandas as pd
import glob
from rpy2 import robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
from rpy2.robjects import Formula
from rpy2.robjects.vectors import FloatVector
import sys


def glm(path, s='S'):
    #activate R functions in Python using rpy2
    pandas2ri.activate()
    base = importr('base')
    stats = importr('stats')
    lme4 = importr('lme4')
    lmtest = importr('lmtest')
    R = ro.r

    if s == 'S':
        season = 1
    else:
        season = 0

    df = pd.read_csv(path)


    #remap the the group labels into numerical values
    df['Group'] = df['Group'].map({'Case': 1, 'Healthy Control': 0})
    df['Season'] = df['Season'].map({'S': 1, 'W': 0})

    #the csv files contained the return value names of bctpy, statsmodels.MixedLM does not seem to like -:_ characters,
    #thus we rename all the dictionary keys so that statsmodels.MixedLM is happy
    df.rename(index=str, columns={'charpath-lambda': 'charpath', 'avg_clustering_coef_wu:C' : 'cc','efficiency_wei-Eglob' : 'effGlob', \
                                   'assortativity_wei-r' : 'assor', 'modularity_und-Q' : 'modul', 'transitivity_wu-T' : 'trans', \
                                   'small_worldness:S' : 'smwor',   'Unnamed: 0' : 'ID'  }, inplace=True)



    #get only data from one season (should this be a random effect?)
    groups = df.groupby(['Season'])
    #change to 0 for summer, 1 for winter
    df = groups.get_group(season)


    #Global efficiency makes model unpredictable, so we remove it
    full = 'Group ~ assor + charpath + cc + modul + trans + smwor'
    a = stats.glm(full, family = 'binomial', data = df)
    print(R.summary(a))


    #testing for trans
    print('TRANSITIVITY')
    null1 = 'Group ~ assor + charpath + cc + modul + smwor'
    b = stats.glm(null1, family = 'binomial', data = df)
    lr_t = lmtest.lrtest(a, b)
    print(lr_t)


    #testing for smwor
    print('SMALL WORLDNESS')
    null2 = 'Group ~ assor + charpath + cc + modul + trans'
    c = stats.glm(null2, family = 'binomial', data = df)
    lr_s = lmtest.lrtest(a, c)
    print(lr_s)

    #testing for modul
    print('MODULARITY')
    null3 = 'Group ~ assor + charpath + cc + trans + smwor'
    d = stats.glm(null3, family = 'binomial', data = df)
    lr_m = lmtest.lrtest(a, d)
    print(lr_m)

    #testing for clustering coefficient
    print('CLUSTERING COEFFICIENT')
    null4 = 'Group ~ assor + charpath + modul + trans + smwor'
    e = stats.glm(null4, family = 'binomial', data = df)
    lr_cc = lmtest.lrtest(a, e)
    print(lr_cc)

    #testing for charpath
    print('CHARACTERISTIC PATH LENGTH')
    null5 = 'Group ~ assor + cc + modul + trans + smwor'
    f = stats.glm(null5, family = 'binomial', data = df)
    lr_cpl = lmtest.lrtest(a, f)
    print(lr_cpl)

    #testing for assor
    print('ASSORTATIVITY')
    null6 = 'Group ~ charpath + cc + modul + trans + smwor'
    g = stats.glm(null6, family = 'binomial', data = df)
    lr_a = lmtest.lrtest(a, g)
    print(lr_a)

    sum_win = 'DUMMY'

    if season == 0:
        sum_win = 'Winter'
    else:
        sum_win = 'summer'

    print(' ')
    print('GLM performed with season: ' + str(sum_win))

    print(' ')
    n = 12
    print('Adjusted p-values after Bonferroni correction for multiple comparison testing with '\
           + str(n) + ' tests')
    lst = [lr_a[4][1],lr_cpl[4][1],lr_cc[4][1],lr_m[4][1],lr_s[4][1],lr_t[4][1]]
    flist = [float(i) for i in lst]
    flist = [round(i,4) for i in flist]

    print(flist)
    print(stats.p_adjust(FloatVector(flist), method='bonferroni', n=n))








