function D2Int=dmd_D2Integral(N,L,x,b,w)
%% 
% D2Integral.m calculates the value, which need to be integrated from 0 to
% bmax(upper bound) to calculate D2, when dimention is odd.
%
% This routine will only be run when the dimension is not even.
%
% Parameters:
%   N: dimension
%   L: number of component for each dimension
%   x: location of dirac distributions
%   b: form 0 to upper bound of integral
%   w: weight of dirac distribution
%
% Return values:
%   D2Int: value needed to be integrated
%
% History:
% * 2018/01/04 Dantong Wang

%% calculate D2Int
omega=b.^(1-N);
sigma=zeros(N,1)+1;
PI=prod(1./((sigma.^2+2.*b.^2).^(1/2)));

SIGMA=0;
for i=1:L
    s1=0;
    for j=1:N
        s1=s1+x(j,i)^2./(sigma(j)^2+2.*b.^2);
    end
    SIGMA=SIGMA+exp((-1/2)*s1);
end
P2=(2*pi).^(N/2).*b.^(2*N).*PI.*w.*SIGMA;
D2Int=omega.*P2;
end