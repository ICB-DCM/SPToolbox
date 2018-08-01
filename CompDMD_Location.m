function [opt]=CompDMD_Location(N,L)
%% 
% CompDMD_Location.m provides a method to calculate Sigma Points for a
%   user-defined dimensionality and number of components for each dimension.
%
%   This function calculates the Sigma Points by minimizing the distance 
%   between dirac mixture approximation and normal distribution in
%   terms of Localized Cumulative Distribution(LCD). Multistart method with
%   20 strart points is used here.
%
% Parameters:
%   N: dimension of dirac mixture approximation
%   L: number of component for each dimension 
%
% Return values:
%   optx: optimizition result of the location of dirac distribution
%
% Additional toolbox needed:
%   Pesto: https://github.com/ICB-DCM/PESTO
%
% History:
% * 2018/01/04 Dantong Wang

%% generate parameter field for getMultiStarts.m
parameters.number = N*(L-1);
parameters.min = -3*ones(1,parameters.number);
parameters.max = 3*ones(1,parameters.number);

%% Log-likelihood function
objectiveFunctionx = @(x) distanceTrueMean(x,N,L);

%% Set options
optionsPesto = PestoOptions();
optionsPesto.n_starts = 20;
optionsPesto.obj_type = 'negative log-posterior';
optionsPesto.save=false;
%optionsPesto.foldername=sprintf('%s%i%s%i','MultiInfo\dim',N,'points',L);

%% optimize locations
parameters = getMultiStarts(parameters, objectiveFunctionx, optionsPesto);
optx = reshape(parameters.MS.par(:,1),[N,L-1]);
optx = [optx,-sum(optx,2)];

%% optimize weights for covariance
parametersw.number = L-1;
parametersw.min = zeros(1,parametersw.number);
parametersw.max = ones(1,parametersw.number);
parametersw.constraints.A = ones(1,parametersw.number);
parametersw.constraints.b = 1;

%% Log-likelihood function
objectiveFunctionw = @(w_c) obj_variance(w_c,optx,N,L);

%% Set options
optionsPesto = PestoOptions();
optionsPesto.localOptimizer = 'fmincon';
optionsPesto.n_starts = 20;
optionsPesto.obj_type = 'negative log-posterior';
optionsPesto.save=false;

%% results
parametersw = getMultiStarts(parametersw, objectiveFunctionw, optionsPesto);
optw = reshape(parametersw.MS.par(:,1),[1,L-1]);
optw = [optw,1-sum(optw,2)];
opt = [optx;optw];
[SPToolboxFolder,~,~]=fileparts(which('CompDMD_Location'));
filename=sprintf('%s%i%s%i%s','DMDw_cPosInfo\B_SP_dim',N,'points',L,'.csv');
dlmwrite(fullfile(SPToolboxFolder,filename),opt,'delimiter',',','precision',12);
end