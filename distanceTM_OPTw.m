function [Dt,Gt] = distanceTM_OPTw(xw,N,L)
% distanceTM_OPTw.m forces the mean of Dirac mixtrue distribution approximation
% to be 0 and non-uniform weights are implemented.
%
% Parameters:
%   xw: location of dirac distributions and the non-uniform weights
%   N: dimension of dirac mixture approximation
%   L: number of component for each dimension 
%
% Return values:
%   Dt: distance between dirac mixture approximation and normal distribution
%   Gt: gradient of Dt
%
% History:
% * 2019/04/30 Dantong Wang
%%
x = reshape(xw(1:N*(L-1)),[N,L-1]);
w = reshape(xw(N*(L-1)+1:(N+1)*(L-1)),[1,(L-1)]);
W = [w,1-sum(w,2)];
X = [x,-sum(w.*x,2)/(1-sum(w,2))];
XW = [reshape(X,[1,N*L]),W];
[D,G] = distanceDiracGaussianOPTw(XW,N,L);
Dt = D;
dXdx = repmat(permute([eye(L-1),(-w/(1-sum(w,2)))'],[3,2,1]),N,1,1);
Gtx = permute(sum(bsxfun(@times, reshape(G(1:N*L),[N,L]),dXdx),2),[1,3,2]);
Gtx = reshape(Gtx,[1,N*(L-1)]);
dWdw = [eye(L-1),-ones(L-1,1)];
Gtw = sum(bsxfun(@times,G(N*L+1:(N+1)*L),dWdw),2)';
Gt = [Gtx,Gtw];