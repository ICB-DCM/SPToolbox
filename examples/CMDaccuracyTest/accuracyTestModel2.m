function [f,dfdphi] = accuracyTestModel2(phi)

f = [phi(1)^4/(1+phi(1)^4),phi(2)^4/(1+phi(2)^4)];

dfdphi(:,:,1) = [4*phi(1)^3/((1+phi(1)^4)^2),0];
dfdphi(:,:,2) = [0,4*phi(2)^3/((1+phi(2)^4)^2)];

end