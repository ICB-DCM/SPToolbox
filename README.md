# About SPToolbox
SPToolbox allows for different Dirac mixture distribution methods to approximate Gaussian distributions. By propagate the approximated Dirac points through linear or nonlinear functions, mean and covariance matrix of the responses can be computed.
# Functions
## getSigmaPointApp.m
Function to compute mean and covariance matrix of the function.
Methods implemented: cmd method, Monte Carlo method, Halton, Sobol quasi Monte Carlo methods, Julier1, Julier2, Menegaz, Lerner, Charalampidis and Merwe sigma point methods.
## CompCMD_Location.m
Function to compute the locations of CMD points, needs to be run before using the 'cmd'method in function 'getSigmaPointApp.m'. Locations are saved in the folder 'CMDTrueMeanInfo' to be used later.
# Non-uniform weights
Non-uniform weighting is also implemented by optimizing weight and the locations simultaneously. To use non-uniform weights, different corresponding functions need to be used:
'getSigmaPointAppOPTw.m' instead of 'getSigmaPointApp.m',
'CompCMD_LocationTM_OPTw.m' instead of 'CompCMD_Location.m'.
Locations and non-uniform weights are saved in the folder 'CMDTM_OPTwInfo' to be used later.
