function G1Int=dmd_G1Integral(N,x,b)
%% 
% dmd_G1Integral.m calculates the value, which need to be integrated from 0 to
% bmax(upper bound) to get the value of G1, when dimention is odd.
%
% This routine will only be run when the dimension is not even.
%
% Parameters:
%   N: dimension
%   x: location of dirac distributions
%   b: form 0 to upper bound of integral
%
% Return values:
%   G1Int: value needed to be integrated
%
% History:
% * 2018/01/05 Dantong Wang

%% 
sigma=zeros(N,1)+1;
PI=prod(1./((sigma.^2+2*b.^2).^(1/2)));
SIGMA=exp((-1/2)*sum(x.^2./(sigma.^2+2*b.^2)));
G1Int=b.^(N+1)./(sigma.^2+2*b.^2)*PI*SIGMA;
end