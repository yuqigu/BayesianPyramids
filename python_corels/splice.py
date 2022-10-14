#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script performs downstream classification using the CORELS package 
for interpretable rule-list classification

@author: yuqi
"""

# --- splice gene data --- #

from corels import *
import numpy as np
import csv
import os
import random
import time

random.seed(1)
os.chdir("/Users/yuqigu/Dropbox/Since2020/Code/Bayesian_Pyramids/corels")

# Bayesian Pyramid features
X, Y1, _, _ = load_from_csv("Y1_splice_K_upper7_alpha2_eff.csv")
_, Y2, _, _ = load_from_csv("Y2_splice_K_upper7_alpha2_eff.csv")
_, Y3, _, _ = load_from_csv("Y3_splice_K_upper7_alpha2_eff.csv")


# DX2009 features
X_dx2009, _, _, _ = load_from_csv("Y1_splice_dx2009_k10.csv")


# IndepBinary features
X_indepbin, _, _, _ = load_from_csv("Y1_splice_csp_indep.csv")


# RAW sequence data
X_all, _, _, _ = load_from_csv("splice_bin_Y1.csv")


# -- train-test split -- #
n, p = X.shape
num_train = np.ceil(n*0.8).astype(int)
train_ind = np.array(random.sample(range(n), num_train));
test_ind = np.setdiff1d(range(n), train_ind);
#
Y1_train = Y1[train_ind]; Y1_test = Y1[test_ind]; 
Y2_train = Y2[train_ind]; Y2_test = Y2[test_ind]; 
Y3_train = Y3[train_ind]; Y3_test = Y3[test_ind]; 


# Create the model, with 10000 as the maximum number of iterations 
# c = CorelsClassifier(max_card=2, n_iter=10000)
c = CorelsClassifier(max_card=3, n_iter=10000)


# -- 1. RAW sequence: look at splice data itself for building rule-lists -- #
start_time = time.time()
X_all_train = X_all[train_ind]; X_all_test = X_all[test_ind]

# Fit, and score the model on the training set
# Print the model's accuracy on the training set
raw_train1 = c.fit(X_all_train, Y1_train).score(X_all_train, Y1_train)
#print(raw_train1)

raw_train2 = c.fit(X_all_train, Y2_train).score(X_all_train, Y2_train)
#print(raw_train2)

raw_train3 = c.fit(X_all_train, Y3_train).score(X_all_train, Y3_train)
#print(raw_train3)


# Test set accuracy
# c_orig_train = CorelsClassifier(max_card=2, n_iter=10000)
raw_test1 = c.fit(X_all_train, Y1_train).score(X_all_test, Y1_test)
#print(raw_test1); print('\n')
# print(np.mean(orig1.predict(X_all_test) == Y1_test))

raw_test2 = c.fit(X_all_train, Y2_train).score(X_all_test, Y2_test)
#print(raw_test2); print('\n')
# print(np.mean(orig2.predict(X_all_test) == Y2_test))

raw_test3 = c.fit(X_all_train, Y3_train).score(X_all_test, Y3_test)
#print(raw_test3); print('\n')
# print(np.mean(orig3.predict(X_all_test) == Y3_test))

print("--- %s seconds for RAW ---" % (time.time() - start_time))



# -- 2. DX2009 -- #
start_time = time.time()
X_train = X_dx2009[train_ind]; X_test = X_dx2009[test_ind]

# Fit, and score the model on the training set
dx1_train = c.fit(X_train, Y1_train).score(X_train, Y1_train)
# print(dx1_train); print('\n')

dx2_train = c.fit(X_train, Y2_train).score(X_train, Y2_train)
# print(dx2_train); print('\n')

dx3_train = c.fit(X_train, Y3_train).score(X_train, Y3_train)
# print(dx3_train); print('\n')


# test set accuracy
dx1_test = c.fit(X_train, Y1_train).score(X_test, Y1_test)
# print(dx1_test); print('\n')

dx2_test = c.fit(X_train, Y2_train).score(X_test, Y2_test)
# print(dx2_test); print('\n')

dx3_test = c.fit(X_train, Y3_train).score(X_test, Y3_test)
# print(dx3_test); print('\n')
print("--- %s seconds for DX2009 ---" % (time.time() - start_time))


# -- 3. IndepBinary -- #
start_time = time.time()
X_train = X_indepbin[train_ind]; X_test = X_indepbin[test_ind]

# Fit, and score the model on the training set
indep1_train = c.fit(X_train, Y1_train).score(X_train, Y1_train)
# print(indep1_train); print('\n')

indep2_train = c.fit(X_train, Y2_train).score(X_train, Y2_train)
# print(indep2_train); print('\n')

indep3_train = c.fit(X_train, Y3_train).score(X_train, Y3_train)
# print(indep3_train); print('\n')


# test set accuracy
indep1_test = c.fit(X_train, Y1_train).score(X_test, Y1_test)
# print(indep1_test); print('\n')

indep2_test = c.fit(X_train, Y2_train).score(X_test, Y2_test)
# print(indep2_test); print('\n')

indep3_test = c.fit(X_train, Y3_train).score(X_test, Y3_test)
# print(indep3_test); print('\n')
print("--- %s seconds for IndepBinary ---" % (time.time() - start_time))



# -- 4. Bayesian Pyramid -- #
start_time = time.time()
X_train = X[train_ind]; X_test = X[test_ind]

# Fit, and score the model on the training set
bp1_train = c.fit(X_train, Y1_train).score(X_train, Y1_train)
# print(bp1_train); print('\n')

bp2_train = c.fit(X_train, Y2_train).score(X_train, Y2_train)
# print(bp2_train); print('\n')

bp3_train = c.fit(X_train, Y3_train).score(X_train, Y3_train)
# print(bp3_train); print('\n')


# test set accuracy
bp1_test = c.fit(X_train, Y1_train).score(X_test, Y1_test)
# print(bp1_test); print('\n')

bp2_test = c.fit(X_train, Y2_train).score(X_test, Y2_test)
# print(bp2_test); print('\n')

bp3_test = c.fit(X_train, Y3_train).score(X_test, Y3_test)
# print(bp3_test); print('\n')
print("--- %s seconds for BayesPyramid ---" % (time.time() - start_time))



# ---- look at accuracy table ---- %
print('\n\nRaw sequence accuracy')
print(np.round(np.array([raw_train1, raw_train2, raw_train3, raw_test1, raw_test2, raw_test3]), 3))

print('\n\nDX2009 accuracy')
print(np.round(np.array([dx1_train, dx2_train, dx3_train, dx1_test, dx2_test, dx3_test]), 3))


print('\n\nIndepBinary accuracy')
print(np.round(np.array([indep1_train, indep2_train, indep3_train, indep1_test, indep2_test, indep3_test]), 3))

print('\n\nBayesPyramid accuracy')
print(np.round(np.array([bp1_train, bp2_train, bp3_train, bp1_test, bp2_test, bp3_test]), 3))


# ###### make Y a three-category variable ######
# Y = np.zeros((np.size(Y1), 1))
# Y[Y1 == 1] = 1
# Y[Y2 == 1] = 2
# Y[Y3 == 1] = 3
# Y = Y.flatten()

# a_all = c.fit(X, Y).score(X, Y)




# ###################### Promoters Dataset ########################
# from corels import *
# import csv
# # Load the dataset
# X, Y, _, _ = load_from_csv("Y_promoter_K_upper7_alpha2_eff.csv")


# # Create the model, with 10000 as the maximum number of iterations 
# c = CorelsClassifier(max_card=2, n_iter=10000)

# # Fit, and score the model on the training set
# a1 = c.fit(X, Y).score(X, Y)

# # Print the model's accuracy on the training set
# print(a1)





