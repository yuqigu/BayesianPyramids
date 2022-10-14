# -*- coding: utf-8 -*-
"""
Created on Wed Oct 28 22:29:22 2020

@author: yuqi
"""

from corels import CorelsClassifier


# --- toy example --- #
# ["loud", "samples"] is the most verbose setting possible
C = CorelsClassifier(max_card=2, c=0.0, verbosity=["loud", "samples"])


# 4 samples, 3 features
X = [[1, 0, 1], [0, 0, 0], [1, 1, 0], [0, 1, 0]]
y = [1, 0, 0, 1]
# Feature names
features = ["Mac User", "Likes Pie", "Age < 20"]



# Fit the model
C.fit(X, y, features=features, prediction_name="Has a dirty computer")



# Print the resulting rulelist
print(C.rl())

# Predict on the training set
print(C.predict(X))


# --- compass data --- #
from corels import *
# Load the dataset
X, y, _, _ = load_from_csv("compas.csv")

# Create the model, with 10000 as the maximum number of iterations 
c = CorelsClassifier(n_iter=10000)

# Fit, and score the model on the training set
a = c.fit(X, y).score(X, y)

# Print the model's accuracy on the training set
print(a)



