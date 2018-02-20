function J0l=dmd_J0l(x,b,l)
%% 
% dmd_J0l.m calculates the value of J0l, which is used in the calculation
% of Jmm
%
% Parameters:
%   x: location of dirac distributions
%   b: upper bound for integral
%   l: from 0 to m, where m is dimension D/2
%
% Return values:
%   J0l: value used in the calculation of Jmm
%
% History:
% * 2018/01/04 Dantong Wang

%% calculate J01
c=sum(x.^2);
divc=-c./(2+4*b^2);
ex_in=ei(divc);
if l==0
    J0l=(1+2*b^2)/4.*exp(divc)+c./8.*ex_in;
else if l==1
        J0l=-1/4.*ex_in;
    else
        SIGMA=0;
        for j=2:l
            SIGMA=SIGMA+(factorial(l-2)*2^(l-j-1))./(factorial(l-2).*c.^(l-j+1).*(1+2*b^2)^(j-2));
        end
        J0l=exp(divc).*SIGMA;
    end
end
end