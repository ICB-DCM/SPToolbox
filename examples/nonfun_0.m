function [f,dfdphi] = nonfun_0(phi)

t = [0,1,2,3];

f = [(phi(1)+phi(2)+phi(3))^2*t(1)^2;
     (phi(1)+phi(2)+phi(3))^2*t(2)^2;
     (phi(1)+phi(2)+phi(3))^2*t(3)^2;
     (phi(1)+phi(2)+phi(3))^2*t(4)^2];

dfdphi(:,:,1) = [2*(phi(1)+phi(2)+phi(3))*t(1)^2;
                 2*(phi(1)+phi(2)+phi(3))*t(2)^2;
                 2*(phi(1)+phi(2)+phi(3))*t(3)^2;
                 2*(phi(1)+phi(2)+phi(3))*t(4)^2];
dfdphi(:,:,2) = [2*(phi(1)+phi(2)+phi(3))*t(1)^2;
                 2*(phi(1)+phi(2)+phi(3))*t(2)^2;
                 2*(phi(1)+phi(2)+phi(3))*t(3)^2;
                 2*(phi(1)+phi(2)+phi(3))*t(4)^2];
dfdphi(:,:,3) = [2*(phi(1)+phi(2)+phi(3))*t(1)^2;
                 2*(phi(1)+phi(2)+phi(3))*t(2)^2;
                 2*(phi(1)+phi(2)+phi(3))*t(3)^2;
                 2*(phi(1)+phi(2)+phi(3))*t(4)^2];

% f(:,1) = [phi(1);
%           phi(2)^2;
%          (phi(1)+phi(2))^2;
%          phi(3)^3];
% 
% dfdphi = [1,0,0;
%           0,2*phi(2),0;
%           2*(phi(1)+phi(2)),2*(phi(1)+phi(2)),0
%           0,0,3*phi(3)^2];

end