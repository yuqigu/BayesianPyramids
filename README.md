# Identifying Interpretable Discrete Latent Structures from Discrete Data


This repository is associated with the paper [Gu, Y. and Dunson, D.B. (2021), Bayesian Pyramids: Identifiable Multilayer Discrete Latent Structure Models for Discrete Data](https://arxiv.org/abs/2101.10373).


### Abstract
High dimensional categorical data are routinely collected in biomedical and social sciences. It is of great importance to build interpretable parsimonious models that perform dimension reduction and uncover meaningful latent structures from such discrete data. Identifiability is a fundamental requirement for valid modeling and inference in such scenarios, yet is challenging to address when there are complex latent structures. In this article, we propose a class of identifiable multilayer (potentially deep) discrete latent structure models for discrete data, termed Bayesian pyramids. We establish the identifiability of Bayesian pyramids by developing novel transparent conditions on the pyramid-shaped deep latent directed graph. The proposed identifiability conditions can ensure Bayesian posterior consistency under suitable priors. As an illustration, we consider the two-latent-layer model and propose a Bayesian shrinkage estimation approach. Simulation results for this model corroborate the identifiability and estimability of model parameters. Applications of the methodology to DNA nucleotide sequence data uncover useful discrete latent features that are highly predictive of sequence types. The proposed framework provides a recipe for interpretable unsupervised learning of discrete data, and can be a useful alternative to popular machine learning methods.


### For simulations:
Note that all Matlab functions and source code files are under `matlab_source_code/`. Run `generate_2layer_data_strong.m` to generate simulated data. Then run the function `bp_csp_simu(K, alpha0)` to perform simulations with Gibbs sampling to estimate the parameters. See the meaning and discussion of K and alpha0 in the paper. Then run `simu_post_process.m` and `simu_figure.m` to evaluate and visualize the estimation results from the simulations.


### For real data analysis:
Take the splice junction dataset analyzed in the paper as an example. Run the function `splice_csp(K, alpha0)` to analyze the splice junction data using the proposed method. Then use `splice_processing.m` and `python_corels/splice.py` to perform downstream classification using the rule-list classifier. Then run `splice_csp_figure.m` to generate figures to evaluate and visualize the data analysis results.


### Generate some simple data and then run the proposed method on it:
Matlab script `demo_run.m` is a simple, readable, and short script which first generates a simulated dataset `Y_data.csv`, and then runs the proposed method on it using the `bayes_pyramid('Y_data.csv')` function and returns the `saved_file` containing all the estimation results. Finally, the last line of the script `load(saved_file)` loads all the estimation results into the current Matlab workspace.

