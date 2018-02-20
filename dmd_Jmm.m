function Jmm=dmd_Jmm(m,x,b)
%% 
% dmd_Jmm.m calculates the value of Jmm, which is used in the calculation
% of D2
%
% Parameters:
%   m: when dimension D is even, m=D/2
%   x: location of dirac distributions
%   b: upper bound for integral
%
% Return values:
%   Jmm: value used in the calculation of D2
%
% History:
% * 2018/01/04 Dantong Wang

%% calculate Jmm
SIGMA=0;
for j=0:m
    SIGMA=SIGMA+(-1)^j*nchoosek(m,j).*dmd_J0l(x,b,j);
end
Jmm=1/2^m.*SIGMA;

end

        