
#https://stackoverflow.com/questions/43043519/is-it-possible-to-do-glmm-in-python
#https://stats.idre.ucla.edu/other/mult-pkg/introduction-to-generalized-linear-mixed-models/

#use these in terminal on Mac(maybe UNIX also) if getting alot of wierd warnings from rpy2
# export LANG=en_US.UTF-8
# export LC_ALL=en_US.UTF-8
# locale

#packages used: pandas, glob and rpy2

import pandas as pd
import glob
from rpy2 import robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
from rpy2.robjects import Formula
import sys

#check if season and parcellation scheme is given
if len(sys.argv) != 3:
    print('Usage: python [season:0/1] [parcellationscheme]')
    exit()

#get season from input
season = int(sys.argv[1])
#get parcellation scheme from input
scheme = sys.argv[2]

#activate R functions in Python using rpy2
pandas2ri.activate()
base = importr('base')
stats = importr('stats')
lme4 = importr('lme4')
R = ro.r

#load in all the estimate csv files as pandas dataframes
filename = 'auto_results/' + scheme + '/estimate.??.csv'

fl = glob.glob(filename)
frames = []

for item in fl: 
    tdf = pd.read_csv(item) 
    frames.append(tdf)

#concatenate them into one large dataframe
df = pd.concat(frames)

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
alt = 'Group ~ charpath + cc + assor + modul + trans + smwor + (1 | Threshold)'
a = lme4.glmer(alt, family = 'binomial', data = df, REML='FALSE')
print(R.summary(a))

#testing for trans
print('TRANSITIVITY')
null1 = 'Group ~ charpath + cc + assor + modul + smwor + (1 | Threshold)'
b = lme4.glmer(null1, family = 'binomial', data = df, REML='FALSE')
anova_c = stats.anova(a, b)
print(anova_c)

#testing for smwor
print('SMALL WORLDNESS')
null2 = 'Group ~ charpath + cc + assor + modul + trans + (1 | Threshold)'
c = lme4.glmer(null2, family = 'binomial', data = df, REML='FALSE')
anova_b = stats.anova(a, c)
print(anova_b)

#testing for modul
print('MODULARITY')
null3 = 'Group ~ charpath + cc + assor + trans + smwor + (1 | Threshold)'
d = lme4.glmer(null3, family = 'binomial', data = df, REML='FALSE')
anova_d = stats.anova(a, d)
print(anova_d)

#testing for clustering coefficient
print('CLUSTERING COEFFICIENT')
null4 = 'Group ~ assor + charpath + modul + trans + smwor + (1 | Threshold)'
e = lme4.glmer(null4, family = 'binomial', data = df, REML='FALSE')
anova_e = stats.anova(a, e)
print(anova_e)

#testing for charpath
print('CHARACTERISTIC PATH LENGTH')
null5 = 'Group ~ assor + cc + modul + trans + smwor + (1 | Threshold)'
f = lme4.glmer(null5, family = 'binomial', data = df, REML='FALSE')
anova_f = stats.anova(a, f)
print(anova_f)

#testing for assor
print('ASSORTATIVITY')
null6 = 'Group ~ charpath + cc + modul + trans + smwor + (1 | Threshold)'
g = lme4.glmer(null6, family = 'binomial', data = df, REML='FALSE')
anova_g = stats.anova(g, a)
print(anova_g)

sum_win = 'DUMMY'

if season == 0:
    sum_win = 'Winter'
else:
    sum_win = 'summer'

print(' ')
print('GLMM performed with parcellation: ' +str(scheme) + ' and with season: ' + str(sum_win))





# result = R.round(base.summary(s, digits=4).rx2('coefficients'),7)
# print(result)














