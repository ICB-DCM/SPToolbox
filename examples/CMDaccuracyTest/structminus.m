function dApp_MCP = structminus(App, MCP)
%% 
% structminus.m computes absolute difference between approximation and true
% values by computing each fields respectively.
%
% Parameters:
%   App: approximated values
%   MCP: true values computed by Monte Carlo method
%
% Return values:
%   dApp_MC: absolute difference between approximation and true values
dApp_MCP.my = abs(App.my - MCP.my);
dApp_MCP.Cy = abs(App.Cy - MCP.Cy);
end