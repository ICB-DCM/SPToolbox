function [f,dfdphi] = accuracyTestModel1(phi)
% accuracyTestModel1.m is a quadratic model used to test the accuracy of 
% different Dirac mixture distribution approximations.
%

f = [phi(1)^2,phi(2)^2];

dfdphi(:,:,1) = [2,0];
dfdphi(:,:,2) = [0,2];

end