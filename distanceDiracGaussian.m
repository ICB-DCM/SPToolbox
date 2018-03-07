function [D,G] = distanceDiracGaussian(x,N,L)
%% 
% distanceDiracGaussian.m calculates the distance between dirac mixture
% approximation and normal distribution in terms of Localized Cumulative
% Distribution. Distance consists of three terms: D1,D2 and D3.
% Detailed inference can be found in:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% "Efficient Deterministic Dirac Mixture Approximation of Gaussian
% Distributions", Igor Gilitschenski and Uwe D. Hanebeck
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Parameters:
%   x: location of dirac distributions
%   N: dimension of dirac mixture approximation
%   L: number of component for each dimension 
%
% Return values:
%   D: distance between dirac mixture approximation and normal distribution
%   G: gradient of D
%
% History:
% * 2018/01/04 Dantong Wang

%% weight and range for integral
x = reshape(x,[N,L]);
w = 1/L;
bmax = 3;

%% calculate D1
D1 = pi^(N/2)*dmd_IN(N,bmax);

%% calculate D2
k = floor(N/2);
l = N - k*2;
% if N is even
if l == 0
    D2 = (2*pi)^(N/2)*sum(w*(dmd_Jmm(k,x,bmax)-dmd_Jmm(k,x,0)),2);
% if N is odd
else
    D2 = integral(@(b)dmd_D2Integral(N,x,b,w),0,3);
end
% calclate Tij as a matrix
Tij = permute(sum(bsxfun(@minus,x,permute(x,[1,3,2])).^2,1),[2,3,1])+1e-5;

%% calculate D3
D31 = bmax^2/2*exp(-1/2*(Tij./(2*bmax^2)));
D32 = (Tij/8).*ei(-1/2*(Tij./(2*bmax^2)));
D3 = pi^(N/2)*sum(sum(w*w*(D31+D32),1),2);

%% Gradient
%% calculate G1
% if N is even
if l == 0
    G1 = 2*(2*pi)^(N/2)*(w*x.*(dmd_Jkk1(x,k,bmax)-dmd_Jkk1(x,k,0)));
% if N is odd
else
    G1 = 2*(2*pi)^(N/2)*w.*x.*integral(@(b)dmd_G1Integral(N,x,b),0,3,'ArrayValue',true);
end

%% calculate G2
%calculate exponential interation of Tij
en_inT = ei(-1/2*(Tij./(2*bmax^2)));
%calculate the sigma in G2
SIGMA = permute(w*sum(w*bsxfun(@minus,permute(x,[2,3,1]),permute(x,[3,2,1])).*en_inT,2),[3,1,2]);
G2 = pi^(N/2)/2.*SIGMA;

%% output
D = D1-2*D2+D3;
G = reshape(G1+G2,[1,N*L]);
end
%%
function D2Int = dmd_D2Integral(N,x,b,w)
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
omega = b.^(1-N);
sigma = zeros(N,1)+1;
PI = prod(1./((sigma.^2+2.*b.^2).^(1/2)));
SIGMA = permute(sum(exp((-1/2)*sum(bsxfun(@rdivide,x.^2,permute((sigma.^2+2.*b.^2),[1,3,2])),1)),2),[1,3,2]);
P2 = (2*pi).^(N/2).*b.^(2*N).*PI.*w.*SIGMA;
D2Int = omega.*P2;
end
%%
function G1Int = dmd_G1Integral(N,x,b)
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
sigma = zeros(N,1)+1;
PI = prod(1./((sigma.^2+2*b.^2).^(1/2)));
SIGMA = exp((-1/2)*sum(x.^2./(sigma.^2+2*b.^2),1));
G1Int = b.^(N+1)./(sigma.^2+2*b.^2)*PI*SIGMA;
end
%%
function IN = dmd_IN(N,b)
%% 
% dmd_IN.m calculates IN, which is used in the calculation of D1
%
% Parameters:
%   N: dimension of dirac mixture approximation
%   b: upper bound for integral
%
% Return values:
%   IN: value used in the calculation of D1
%
% History:
% * 2018/01/04 Dantong Wang
% * 2018/01/25 Dantong Wang

%% I1 and I2
k = floor(N/2);
l = N-k*2;
I1 = 1/2*(b*(1+b^2)^(1/2)-asinh(b));
I2 = 1/2*(b^2-log(1+b^2));

%% calculate IN
%if N <= 4
if N == 1
    IN = I1;
elseif N == 2
    IN = I2;
elseif N == 4
    IN = k*I2-(b^(2*k)*((1+b^2)^(1/2))^(2-2*k))/(2*k-2);
%if N is even
elseif l == 0
    i = (2:k-1)';
    SIGMA = sum((b.^(2*i).*((1+b^2)^(1/2)).^(2-2*i))./(2*i.*(2*i-2)),1);
    IN = 2*k*(I2-SIGMA)-(b^(2*k)*((1+b^2)^(1/2))^(2-2*k))/(2*k-2);
%if N is odd
else
    i = (1:k-1)';
    SIGMA=sum((b.^(2*i+1).*((1+b^2)^(1/2)).^(1-2*i))./((2*i+1).*(2*i-1)),1);
    IN = (2*k+1)*(I1-SIGMA)-(b^(2*k+1)*((1+b^2)^(1/2))^(1-2*k))/(2*k-1);
end
end
%%
function J0l = dmd_J0l(x,b,l)
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
c = sum(x.^2,1);
divc = -c./(2+4*b^2);
ex_in = ei(divc);
if l == 0
    J0l = (1+2*b^2)/4.*exp(divc)+c./8.*ex_in;
elseif l == 1
    J0l = -1/4.*ex_in;
else
    j = (2:l)';
    SIGMA = sum((factorial(l-2)*2.^(l-j-1))./(factorial(l-2).*c.^(l-j+1).*(1+2*b^2).^(j-2)),1);
    J0l = exp(divc).*SIGMA;
end
end
%%
function Jkk1 = dmd_Jkk1(x,k,b)
%% 
% dmd_Jkk1.m calculates the value of Jkk1, which is used in the calculation
% of G1
%
% this routine will only be run if dimension is even.
%
% Parameters:
%   x: location of dirac distributions
%   k: if dimension is even, k = D/2
%   b: upper bound of integral
%
% Return value:
%   Jkk1: value to be used in the calculation of G1
% History:
% * 2018/01/04 Dantong Wang

%% calculate Jkk1
SIGMA = 0;
for j = 0:k
    SIGMA = SIGMA+(-1)^j*nchoosek(k,j).*dmd_J0l(x,b,j+1);
end
Jkk1 = 1/(2^k)*SIGMA;

end
%%
function Jmm = dmd_Jmm(m,x,b)
%% 
% dmd_Jmm.m calculates the value of Jmm, which is used in the calculation
% of D2
%
% Parameters:
%   m: when dimension D is even, m = D/2
%   x: location of dirac distributions
%   b: upper bound for integral
%
% Return values:
%   Jmm: value used in the calculation of D2
%
% History:
% * 2018/01/04 Dantong Wang

%% calculate Jmm
SIGMA = 0;
for j = 0:m
    SIGMA = SIGMA+(-1)^j*nchoosek(m,j).*dmd_J0l(x,b,j);
end
Jmm = 1/2^m.*SIGMA;

end