function [f,dfdphi] = accuracyTestModel3(phi)
%% MODEL
% State vector
% x = 1/Omega*[[X1],[X2]]^T

% Initial condition
x0 = [2
      0
      0
      0
      0
      0];

% Parameter vector
theta = exp(phi);
         
% Reaction fluxes
F = @(t,x) [-theta(1)*x(1)+theta(2)*x(2);
             theta(1)*x(1)-theta(2)*x(2);
            -theta(1)*x(3)+theta(2)*x(4)-x(1);
             theta(1)*x(3)-theta(2)*x(4)+x(1);
            -theta(1)*x(5)+theta(2)*x(6)+x(2);
             theta(1)*x(5)-theta(2)*x(6)-x(2)];

%% SIMULATION WITH ODE15s
[~,x] = ode15s(F,[0:0.1:5],x0);
% f = [x(1:10:51,1),x(1:10:51,2)];
% dfdphi(:,:,1) = [x(1:10:51,3),x(1:10:51,4)];
% dfdphi(:,:,2) = [x(1:10:51,5),x(1:10:51,6)];

f = x(1:10:51,2);
dfdphi(:,:,1) = x(1:10:51,4);
dfdphi(:,:,2) = x(1:10:51,6);
end