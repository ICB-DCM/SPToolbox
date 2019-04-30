function CompCMD_LocationTM_OPTw(N, L)
%% 
% CompCMD_Location.m provides a method to calculate Sigma Points for a
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
%   optx: optimization result of the location of dirac distribution
%
% Additional toolbox needed:
%   Pesto: https://github.com/ICB-DCM/PESTO
%
% History:
% * 2018/04/03 Dantong Wang

%% generate parameter field for getMultiStarts.m
parameters.number = (N+1)*(L-1);
parameters.min = [-3*ones(1,N*(L-1)),zeros(1,(L-1))];
parameters.max = [3*ones(1,N*(L-1)),ones(1,(L-1))];

%% Log-likelihood function
objectiveFunction = @(x) distanceTM_OPTw(x,N,L);

%% Set options
optionsPesto = PestoOptions();
optionsPesto.n_starts = 20;
optionsPesto.obj_type = 'negative log-posterior';
optionsPesto.save = false;
%optionsPesto.foldername=sprintf('%s%i%s%i','MultiInfo\dim',N,'points',L);

%% output
parameters = getMultiStarts(parameters, objectiveFunction, optionsPesto);
optxw = parameters.MS.par(:,1);
optw = optxw(N*(L-1)+1:(N+1)*(L-1))';
optx = reshape(optxw(1:N*(L-1)),[N,L-1]);
optx = [optx,-sum(optw.*optx,2)/(1-sum(optw,2))];
optw = [optw,1-sum(optw,2)];
opt = [optx;optw];
[SPToolboxFolder,~,~]=fileparts(which('CompCMD_Location'));
filepath = fullfile(SPToolboxFolder,'CMDTM_OPTwInfo');
filename=sprintf('%s%i%s%i%s','B_SP_dim',N,'points',L,'.csv');
dlmwrite(fullfile(filepath,filename),opt,'delimiter',',','precision',12);
end
