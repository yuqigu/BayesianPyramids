#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 01:51:43 2020

@author: yuqi
"""

# --- splice gene data --- #

from corels import *
import numpy as np
import csv
# Load the dataset
X, Y1, _, _ = load_from_csv("/Users/yuqi/Dropbox/Research20/Code/multilayer/corels/Y1_splice_K_upper7_alpha2_eff.csv")
_, Y2, _, _ = load_from_csv("/Users/yuqi/Dropbox/Research20/Code/multilayer/corels/Y2_splice_K_upper7_alpha2_eff.csv")
_, Y3, _, _ = load_from_csv("/Users/yuqi/Dropbox/Research20/Code/multilayer/corels/Y3_splice_K_upper7_alpha2_eff.csv")


# Create the model, with 10000 as the maximum number of iterations 
c = CorelsClassifier(max_card=2, n_iter=10000)

# Fit, and score the model on the training set
a1 = c.fit(X, Y1).score(X, Y1)
a2 = c.fit(X, Y2).score(X, Y2)
a3 = c.fit(X, Y3).score(X, Y3)


# Print the model's accuracy on the training set
print(a1)

print(a2)

print(a3)

###### make Y a three-category variable ######
Y = np.zeros((np.size(Y1), 1))
Y[Y1 == 1] = 1
Y[Y2 == 1] = 2
Y[Y3 == 1] = 3
Y = Y.flatten()

a_all = c.fit(X, Y).score(X, Y)



###################### Promoters Dataset ########################
from corels import *
import csv
# Load the dataset
X, Y, _, _ = load_from_csv("/Users/yuqi/Dropbox/Research20/Code/multilayer/corels/Y_promoter_K_upper7_alpha2_eff.csv")


# Create the model, with 10000 as the maximum number of iterations 
c = CorelsClassifier(max_card=2, n_iter=10000)

# Fit, and score the model on the training set
a1 = c.fit(X, Y).score(X, Y)

# Print the model's accuracy on the training set
print(a1)





