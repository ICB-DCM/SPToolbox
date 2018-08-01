function [Dt,Gt] = distanceTrueMean(x,N,L)
%% 
% distanceTrueMean.m constraints the mean of DMD simgma point locations to be 0, 
% while calculating the distance between dirac mixture
% approximation and normal distribution in terms of Localized Cumulative
% Distribution. 
% Parameters:
%   x: location of dirac distributions
%   N: dimension of dirac mixture approximation
%   L: number of component for each dimension 
%
% Return values:
%   Dt: distance between dirac mixture approximation and normal distribution
%   Gt: gradient of D
%
% History:
% * 2018/01/08 Dantong Wang
%% reshape the input x with (L-1) components
x = reshape(x,[N,L-1]);
X = [x,-sum(x,2)];  % the last component is calculated to ensure the mean to be 0
X = reshape(X,[1,N*L]);
%% calculate distance and gradiant
[D,G] = distanceDiracGaussian(X,N,L);
%% calculate distance and gradiant for (L-1) components
Dt = D;
for i=1:L-1
    dXdx(:,:,i) = [zeros(N,i-1),ones(N,1),zeros(N,L-1-i),-ones(N,1)];
end
Gt = permute(sum(bsxfun(@times, reshape(G,[N,L]),dXdx),2),[1,3,2]);
Gt = reshape(Gt,[1,N*(L-1)]);