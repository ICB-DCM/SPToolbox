function [f,dfdphi] = accuracyTestModel1(phi)

f = [phi(1)^2,phi(2)^2];

dfdphi(:,:,1) = [2,0];
dfdphi(:,:,2) = [0,2];

end