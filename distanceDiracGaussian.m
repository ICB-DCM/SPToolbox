function [D,G]=distanceDiracGaussian(x,N,L)
%% 
% distanceDiracGaussian.m calculates the distance between dirac mixture
% approximation and normal distribution in terms of Localized Cumulative
% Distribution.
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
x=reshape(x,[N,L]);
w=1/L;
bmax=3;

%% calculate D1
IN=dmd_IN(N,bmax);
D1=pi^(N/2)*IN;

%% calculate D2
k=floor(N/2);
l=N-k*2;
% if N is even
if l==0
    D2i=w*(dmd_Jmm(k,x,bmax)-dmd_Jmm(k,x,0));
    D2=(2*pi)^(N/2)*sum(D2i);
% if N is odd
else
    D2=integral(@(b)dmd_D2Integral(N,L,x,b,w),0,3);
end
% calclate Tij as a matrix
for i=1:L
    for j=1:L
        Tij(i,j)=sum((x(:,i)-x(:,j)).^2)+1e-5;
    end
end

%% calculate D3
D31=bmax^2/2*exp(-1/2*(Tij./(2*bmax^2)));
D32=(Tij/8).*ei(-1/2*(Tij./(2*bmax^2)));
D3=pi^(N/2)*sum(sum(w*w*(D31+D32)));

%% Gradient
%% calculate G1
% if N is even
if l==0
    J=dmd_Jkk1(x,k,bmax)-dmd_Jkk1(x,k,0);
    G1i=w*x.*J;
    G1=2*(2*pi)^(N/2)*G1i;
% if N is odd
else
    G1i=integral(@(b)dmd_G1Integral(N,x,b),0,3,'ArrayValue',true);
    G1=2*(2*pi)^(N/2)*w.*x.*G1i;
end

%% calculate G2
%calculate exponential interation of Tij
en_inT=ei(-1/2.*(Tij./(2*bmax^2)));
%calculate the sigma in G2
for m=1:N
    for i=1:L
        for j=1:L
            dif(i,j)=x(m,i)-x(m,j);
            sigma(i,j)=w*dif(i,j)*en_inT(i,j);
        end
    end
    tsigma=sigma';
    SIGMA(m,:)=w*sum(tsigma);
end
G2=pi^(N/2)/2.*SIGMA;

%% output
D=D1-2*D2+D3;
G=G1+G2;
G=reshape(G,[1,N*L]);