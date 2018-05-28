import numpy as np
import pandas as pd

a = np.array([[0,0.12,0.13,0.14,0.15,0.16],
	          [0.12,0,0.23,0.24,0.25,0.26],
	          [0.13,0.23,0,0.34,0.35,0.36],
	          [0.14,0.24,0.34,0,0.45,0.46],
	          [0.15,0.25,0.35,0.45,0,0.56],
	          [0.16,0.26,0.36,0.46,0.56,0]])

b = np.array([[0,0.19,0.63,0.74,0.35,0.56],
	          [0.19,0,0.43,0.84,0.25,0.26],
	          [0.63,0.43,0,0.64,0.35,0.36],
	          [0.74,0.84,0.64,0,0.45,0.86],
	          [0.35,0.25,0.35,0.45,0,0.56],
	          [0.56,0.26,0.36,0.86,0.56,0]])



dfa = pd.DataFrame(a)
dfb = pd.DataFrame(b)
dfa.to_csv('anumpy.csv')
dfb.to_csv('bnumpy.csv')

# ok = pd.read_csv('anumpy.csv', sep=',',header=None, skiprows=1)
# pk = np.delete(ok.values, 0, 1)

# print(pk.shape)

