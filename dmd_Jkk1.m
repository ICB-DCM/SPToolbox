function Jkk1=dmd_Jkk1(x,k,b)
%% 
% dmd_Jkk1.m calculates the value of Jkk1, which is used in the calculation
% of G1
%
% this routine will only be run if dimension is even.
%
% Parameters:
%   x: location of dirac distributions
%   k: if dimension is even, k=D/2
%   b: upper bound of integral
%
% Return value:
%   Jkk1: value to be used in the calculation of G1
% History:
% * 2018/01/04 Dantong Wang

%% calculate Jkk1
SIGMA=0;
for j=0:k
    SIGMA=SIGMA+(-1)^j*nchoosek(k,j).*dmd_J0l(x,b,j+1);
end
Jkk1=1/(2^k)*SIGMA;

end