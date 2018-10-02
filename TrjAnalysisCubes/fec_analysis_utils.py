from random import choice
import math
from scipy import stats as sp
import numpy as np


def Bootstrap(seq):  # return a bootstrap sample with replacement of the supplied sequence
    bsList = []
    for i in range(len(seq)):
        bsList.append(choice(seq))
    return bsList


def GetCI(numlist):  # get 95% interval from a list
    scipy25 = sp.mstats.mquantiles(numlist, 0.025, 0, 0)
    scipy975 = sp.mstats.mquantiles(numlist, 0.975, 0, 0)
    return np.round(scipy25, 3)[0], np.round(scipy975, 3)[0]


def linreg(lis1, lis2):
    slope, intercept, rvalue, pvalue, stderr = sp.linregress(lis1, lis2)
    rsq = float(rvalue ** 2)
    return np.round(slope, 4), np.round(rsq, 4), np.round(stderr, 4)


# Regression coefficient confidence interval according to Fisher.
def GetrCI(rsqd, num):
    stderr = 1.0 / math.sqrt(num - 3)
    delta = 1.96 * stderr  # 1.96 = 95% CI for r. Set other multiplier for different CI.
    lower = math.tanh(math.atanh(rsqd) - delta)
    upper = math.tanh(math.atanh(rsqd) + delta)

    return np.round(lower, 3), np.round(upper, 3)
