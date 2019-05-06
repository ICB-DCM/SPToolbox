function [my, Cy] = average_time(dApp_MCt)
%% 
% average_difference.m computes the average the absolute differences of 
% mean and covariance matrix across measured time points.
%
% Parameters:
%   dApp_MCt: absolute difference between approximation and true values
%
% Return values:
%   my: average value of mean
%   Cy: average value of covariance matrix
my = sum(dApp_MCt.my, 1)/size(dApp_MCt.my, 1);
Cy = sum(dApp_MCt.Cy, 1)/size(dApp_MCt.Cy, 1);
end