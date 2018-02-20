function IN=dmd_IN(N,b)
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
k=floor(N/2);
l=N-k*2;
I1=1/2*(b*(1+b^2)^(1/2)-asinh(b));
I2=1/2*(b^2-log(1+b^2));

%% calculate IN
%if N<02
if N==1
    IN=I1
else if N==2 
        IN=I2;
        %if N is even
    else if l==0
            if N==4
                IN=k*I2-(b^(2*k)*((1+b^2)^(1/2))^(2-2*k))/(2*k-2);
            else
                SIGMA=0;
                for i=2:k-1
                    SIGMA=SIGMA+(b^(2*i)*((1+b^2)^(1/2))^(2-2*i))/(2*i*(2*i-2));
                end
            IN=2*k*(I2-SIGMA)-(b^(2*k)*((1+b^2)^(1/2))^(2-2*k))/(2*k-2);
            end
            %if N is odd
        else
            SIGMA=0;
            for i=1:k-1
                SIGMA=SIGMA+(b^(2*i+1)*((1+b^2)^(1/2))^(1-2*i))/((2*i+1)*(2*i-1));
            end
            IN=(2*k+1)*(I1-SIGMA)-(b^(2*k+1)*((1+b^2)^(1/2))^(1-2*k))/(2*k-1);
        end
    end
end
end