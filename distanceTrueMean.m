function [distanceTrue,gradientTrue] = distanceTrueMean(locations,N,L)
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
%   distanceTrue: distance between dirac mixture approximation and normal distribution
%   gradientTrue: gradient of D
%
% History:
% * 2018/01/08 Dantong Wang
%% reshape the input x with (L-1) components
locations = reshape(locations,[N,L-1]);
allLocations = [locations,-sum(locations,2)];  % the last component is calculated to ensure the mean to be 0
allLocations = reshape(allLocations,[1,N*L]);
%% calculate distance and gradiant
[distanceTrue,gradient] = distanceDiracGaussian(allLocations,N,L);
%% calculate distance and gradiant for (L-1) components
for iLocation=1:L-1
    dXdx(:,:,iLocation) = [zeros(N,iLocation-1),ones(N,1),zeros(N,L-1-iLocation),-ones(N,1)];
end
gradientTrue = permute(sum(bsxfun(@times, reshape(gradient,[N,L]),dXdx),2),[1,3,2]);
gradientTrue = reshape(gradientTrue,[1,N*(L-1)]);