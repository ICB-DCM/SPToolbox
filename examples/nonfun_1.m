function [f,dfdphi] = nonfun_1(phi)

t = [0,1,2,3];

f = [phi(1)*t(1),(phi(1)+phi(2))^2*t(1)^2,phi(3)^3*t(1)^2;
     phi(1)*t(2),(phi(1)+phi(2))^2*t(2)^2,phi(3)^3*t(2)^2;
     phi(1)*t(3),(phi(1)+phi(2))^2*t(3)^2,phi(3)^3*t(3)^2;
     phi(1)*t(4),(phi(1)+phi(2))^2*t(4)^2,phi(3)^3*t(4)^2];

dfdphi(:,:,1) = [t(1),2*(phi(1)+phi(2))*t(1)^2,0;
                 t(2),2*(phi(1)+phi(2))*t(2)^2,0;
                 t(3),2*(phi(1)+phi(2))*t(3)^2,0;
                 t(4),2*(phi(1)+phi(2))*t(4)^2,0];
dfdphi(:,:,2) = [0,2*(phi(1)+phi(2))*t(1)^2,0;
                 0,2*(phi(1)+phi(2))*t(2)^2,0;
                 0,2*(phi(1)+phi(2))*t(3)^2,0;
                 0,2*(phi(1)+phi(2))*t(4)^2,0];
dfdphi(:,:,3) = [0,0,3*phi(3)^2*t(1)^2;
                 0,0,3*phi(3)^2*t(2)^2;
                 0,0,3*phi(3)^2*t(3)^2;
                 0,0,3*phi(3)^2*t(4)^2];

end