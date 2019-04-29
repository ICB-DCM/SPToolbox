function [my, Cy] = average_time(dApp_MCt)
my = sum(dApp_MCt.my, 1)/size(dApp_MCt.my, 1);
Cy = sum(dApp_MCt.Cy, 1)/size(dApp_MCt.Cy, 1);
end