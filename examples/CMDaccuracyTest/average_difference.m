function [dApp_MCa] = average_difference(dApp_MC, dim)
%% 
% average_difference.m computes the average the absolute differences of 
% multi-dimension mean and covariance matrix.
%
% Parameters:
%   dApp_MC: absolute difference between approximation and true values
%   dim: dimension of mean and covariance matrix
%
% Return values:
%   dApp_MCa: average value
dApp_MCa.my = sum(dApp_MC.my, 2) / dim;
dApp_MCa.Cy = zeros(size(dApp_MC.Cy, 1), 1);
for iCy = 1:size(dApp_MC.Cy, 1)
    Cy = permute(dApp_MC.Cy(iCy, :, :), [3, 2, 1]);
    sum_diag = sum(diag(Cy));
    sum_ndiag = sum(sum(Cy - diag(diag(Cy))))/2;
    dApp_MCa.Cy(iCy) = sum_diag + sum_ndiag;
end
dApp_MCa.Cy = dApp_MCa.Cy / ((dim^2 + dim) / 2);
end