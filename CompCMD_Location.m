function [splocations]=CompCMD_Location(N,L)
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
%   splocations: optimizition result of the location of dirac distribution
%
% Additional toolbox needed:
%   Pesto: https://github.com/ICB-DCM/PESTO
%
% History:
% * 2018/01/04 Dantong Wang

%% generate parameter field for getMultiStarts.m
%parameters.number = N*L;
parameters.number = N*(L-1);
parameters.min = -3*ones(1,parameters.number);
parameters.max = 3*ones(1,parameters.number);

%% Log-likelihood function
objectiveFunction = @(x) distanceTrueMean(x,N,L);

%% Set options
optionsPesto = PestoOptions();
optionsPesto.n_starts = 20;
optionsPesto.obj_type = 'negative log-posterior';
optionsPesto.save=false;
%optionsPesto.foldername=sprintf('%s%i%s%i','MultiInfo\dim',N,'points',L);

%% output
parameters = getMultiStarts(parameters, objectiveFunction, optionsPesto);
splocations = reshape(parameters.MS.par(:,1),[N,L-1]);
splocations = [splocations,-sum(splocations,2)];
[SPToolboxFolder,~,~]=fileparts(which('CompDMD_Location'));
filepath = fullfile(SPToolboxFolder,'DMDTrueMeanInfo');
filename=sprintf('%s%i%s%i%s','B_SP_dim',N,'points',L,'.csv');
dlmwrite(fullfile(filepath,filename),splocations,'delimiter',',','precision',12);
end