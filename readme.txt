This folder contains example code for analysing 2D NMR relaxation experiments in the time domain using the bicomplex signal formulation. The details of the method as been presented in the article "Time-domain signal modelling in multidimensional NMR experiments for estimation of relaxation parameters".

Description of the folder contents:
_data\  --  a folder containing exemplar experimental datasets in the Azara format for T1 (1150), T1rho (1151), and T2 (1155) experiments along with matlab files with their parameters (chemical shifts and peak widths).

bicomplex.m  -- a Matlab class for working with bicomplex data. Credit: Adriaen Verheyleweghen

examples_experimental.m  --  a script for processing experimental spectra. Please choose a specific example to run.

examples_fitting.m  --  a scripts with examples of non-linear model fitting by minimization of the variable projection functional.

getModelMatrix.m  --  a function to compute the model matrix given a set of parameters.

costFunc_VP.m  --  computes the value of variable projection functional based on the difference between the measured data and the model defined by a specific set of parameters.

ffthc.m and iffthc.m  --  functions to compute the forward and inverse Fourier transforms for bicomplex data.
