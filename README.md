# Identifying Interpretable Discrete Latent Structures from Discrete Data


This repository is associated with the paper [Gu, Y. and Dunson, D.B. (2021), Identifying Interpretable Discrete Latent Structures from Discrete Data](https://arxiv.org/abs/2101.10373).

### Abstract
High dimensional categorical data are routinely collected in biomedical and social sciences. It is of great importance to build interpretable models that perform dimension reduction and uncover meaningful latent structures from such discrete data. Identifiability is a fundamental requirement for valid modeling and inference in such scenarios, yet is challenging to address when there are complex latent structures. In this article, we propose a class of interpretable discrete latent structure models for discrete data and develop a general identifiability theory. Our theory is applicable to various types of latent structures, ranging from a single latent variable to deep layers of latent variables organized in a sparse graph (termed a Bayesian pyramid). The proposed identifiability conditions can ensure Bayesian posterior consistency under suitable priors. As an illustration, we consider the two-latent-layer model and propose a Bayesian shrinkage estimation approach. Simulation results for this model corroborate identifiability and estimability of the model parameters. Applications of the methodology to DNA nucleotide sequence data uncover discrete latent features that are both interpretable and highly predictive of sequence types. The proposed framework provides a recipe for interpretable unsupervised learning of discrete data, and can be a useful alternative to popular machine learning methods.

### Main Contents
Main code implementing the methodology is in `source_code/`. The nucleotide sequence data are in `real_data/`.
