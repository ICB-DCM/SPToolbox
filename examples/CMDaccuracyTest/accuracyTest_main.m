function accuracyTest_main(Model)
%% 
% accuracyTest_main.m computes the absulute errors to the true mean and 
% true covariance matrix, and generate a figure showing the trajectories.
%
% Parameters:
%   Model: 'model1' for \phi.^2
%          'model2' for \phi.^4 ./ (1 + \phi.^4)
%          'model3' for conversion reaction model
%
% Return values:
%   figure plots the comparison of Halton, CMD and SP methods
%
% History:
% * 2019/04/29 Dantong Wang
%% switch models
switch Model
    case 'model1'
        model = @(phi)accuracyTestModel1(phi);
    case 'model2'
        model = @(phi)accuracyTestModel2(phi);
    case 'model3'
        model = @(phi)accuracyTestModel3(phi);
end
%% Font
Fontsize = 10;
Marksize = 3;
Linewidth = 1;
%% color defination
color1 = [0.85,0.325,0.098];   %MC
color2 = [0,0.447,0.741];   %Halton
color3 = [0.929,0.694,0.125];   %Sobol
color4 = [0.494,0.184,0.556];   %SP1
color5 = [0.466,0.674,0.188];   %SP2
color6 = [0.301,0.745,0.933];   %SP3
color7 = [0.635,0.078,0.184];   %SP4
color8 = [0.88,0.5,0.7];   %SP5
color9 = [0.3,0.6,0.7];   %SP6
color10 = [0.8,0.1,0.3];   %DMD
%% prepare figure
subplot(1, 2, 1)
set(gca,'FontSize',Fontsize, 'YScale', 'log', 'XScale', 'log')
xlabel('Sample size','FontSize',Fontsize)
ylabel('Error in mean','FontSize',Fontsize)
hold on;
subplot(1, 2, 2)
set(gca,'FontSize',Fontsize, 'YScale', 'log', 'XScale', 'log')
xlabel('Sample size','FontSize',Fontsize)
ylabel('Error in covariance','FontSize',Fontsize)
hold on;
%% Example
% Parameterization
A = eye(2);
B = eye(2);

% xi
switch Model
    case 'model3'
        xi = [log(1);log(2);0;0];
        dim_out = 1;
    otherwise
        xi = [0;0;0;0];
        dim_out = 2;
end

% beta
beta_fun = @(xi) [xi(1);
                  xi(2)];
dbetadxi = [1,0,0,0;
            0,1,0,0];
        
%delta
delta_fun = @(xi) [xi(3);
                  xi(4)];
ddeltadxi = [0,0,1,0;
             0,0,0,1];

estruct.beta = beta_fun;
estruct.dbetadxi = @(xi) dbetadxi;
estruct.delta = delta_fun;
estruct.ddeltadxi = @(xi) ddeltadxi;
estruct.phi = @(beta,b) A*beta + B*b;
estruct.dphidbeta = @(beta,b) A;
estruct.dphidb = @(beta,b) B;
estruct.sigma_noise = @(phi) zeros(1,1);
estruct.dsigma_noisedphi = @(phi) zeros(1,1,2);
dim_in = 2;

%% Results from Monte Carlo
op_SP.nderiv = 0;
op_SP.req = [1,1,0,0,0];
op_SP.type_D = 'diag-matrix-logarithm';       
op_SP.approx = 'samples';
samples = mvnrnd(zeros(2,1),eye(2),100000);
op_SP.samples = samples;
MCP = getSigmaPointApp(model,xi,estruct,op_SP); 
%% Number of points
PointN = [2, 3, 4, 5, 6, 8, 10, 13, 16, 20, 25, 32, 40, 50, 63, 79, 100];
%% Results of Halton quasi Monte Carlo
color = color2;
op_SP.approx = 'halton';
[dHalton_MC, dHalton_MC_neg, dHalton_MC_pos] = compMeanErrorbar(model,xi,estruct,op_SP, PointN, dim_out, dim_in, MCP);
% plot
subplot(1, 2, 1)
errorbar(PointN, dHalton_MC.my, dHalton_MC_neg.my, dHalton_MC_pos.my, '-','Color', color, 'MarkerSize', Linewidth, 'LineWidth', Linewidth)
subplot(1, 2, 2)
errorbar(PointN, dHalton_MC.Cy, dHalton_MC_neg.Cy, dHalton_MC_pos.Cy, '-','Color', color, 'MarkerSize', Linewidth, 'LineWidth', Linewidth)
%% Results of CMD method
color = color10;
op_SP.approx = 'cmd';
[dCMD_MC, dCMD_MC_neg, dCMD_MC_pos] = compMeanErrorbar(model,xi,estruct,op_SP, PointN, dim_out, dim_in, MCP);
% plot
subplot(1, 2, 1)
errorbar(PointN, dCMD_MC.my, dCMD_MC_neg.my, dCMD_MC_pos.my, '-','Color', color, 'MarkerSize', Linewidth, 'LineWidth', Linewidth)
subplot(1, 2, 2)
errorbar(PointN, dCMD_MC.Cy, dCMD_MC_neg.Cy, dCMD_MC_pos.Cy, '-','Color', color, 'MarkerSize', Linewidth, 'LineWidth', Linewidth)
%% Results of non uniform weights CMD methods
color = color5;
op_SP.approx = 'cmd-non-uniform';
[dCMDoptw_MC, dCMDoptw_MC_neg, dCMDoptw_MC_pos] = compMeanErrorbar(model,xi,estruct,op_SP, PointN, dim_out, dim_in, MCP);
% plot
subplot(1, 2, 1)
errorbar(PointN, dCMDoptw_MC.my, dCMDoptw_MC_neg.my, dCMDoptw_MC_pos.my, '-','Color', color, 'MarkerSize', Linewidth, 'LineWidth', Linewidth)
subplot(1, 2, 2)
errorbar(PointN, dCMDoptw_MC.Cy, dCMDoptw_MC_neg.Cy, dCMDoptw_MC_pos.Cy, '-','Color', color, 'MarkerSize', Linewidth, 'LineWidth', Linewidth)
%% Results of sigma point methods
imethod = 1;
n_samples = 1;
while imethod <= 6
    % initialize matrices, for Charalampidis method, more than one points
    % are used
    if imethod == 5
        n_samples = n_samples + 1;
        if n_samples == 2
            dSP_MC.my = zeros(length(2:10), dim_out);
            dSP_MC.Cy = zeros(length(2:10), dim_out, dim_out);
        end
    else
        dSP_MC.my = zeros(1, dim_out);
        dSP_MC.Cy = zeros(1, dim_out, dim_out);
    end
    switch imethod
        case 1
            color = color4;
            marker = '+';
            xaxis = 2 * dim_in + 1;
            op_SP.approx = 'Julier1' ;
            imethod = imethod + 1;
        case 2
            color = color5;
            marker = 'x';
            xaxis = 2 * dim_in + 1;
            op_SP.approx = 'Julier2' ;
            imethod = imethod + 1;
        case 3
            color = color6;
            marker = 'd';
            xaxis = dim_in + 1;
            op_SP.approx = 'Menegaz' ;
            imethod = imethod + 1;
        case 4
            color = color7;
            marker = 'p';
            xaxis = 2 * dim_in ^2 + 1;
            op_SP.approx = 'Lerner' ;
            imethod = imethod + 1;
        case 5
            color = color9;
            marker = 's';
            xaxis = [2:10] .^ dim_in;
            op_SP.approx = 'Charalampidis';
            op_SP.n_samples = n_samples;
            if n_samples == 10
                imethod = imethod + 1;
            end
        case 6
            color = color8;
            marker = 'h';
            xaxis = 2 * dim_in + 1;
            op_SP.approx = 'sp' ;
            imethod = imethod + 1;
    end
    SP = getSigmaPointApp(model,xi,estruct,op_SP);
    dSP_MCt = structminus(SP, MCP);
    if strcmp(op_SP.approx, 'Charalampidis')
        [dSP_MC.my(n_samples - 1, :), dSP_MC.Cy(n_samples - 1, :, :)] =  average_time(dSP_MCt);
        if n_samples == 10
            dSP_MC = average_difference(dSP_MC, dim_out);
            subplot(1, 2, 1)
            plot(xaxis, dSP_MC.my, 'Marker', marker, 'Color', color, 'MarkerSize', Marksize*2, 'LineWidth', Marksize-1, 'LineStyle', 'none')
            subplot(1, 2, 2)
            plot(xaxis, dSP_MC.Cy, 'Marker', marker, 'Color', color, 'MarkerSize', Marksize*2, 'LineWidth', Marksize-1, 'LineStyle', 'none')
        end
    else
        [dSP_MC.my(1, :), dSP_MC.Cy(1, :, :)] =  average_time(dSP_MCt);
        dSP_MC = average_difference(dSP_MC, dim_out);
        subplot(1, 2, 1)
        plot(xaxis, dSP_MC.my, 'Marker', marker, 'Color', color, 'MarkerSize', Marksize*2, 'LineWidth', Marksize-1, 'LineStyle', 'none')
        subplot(1, 2, 2)
        plot(xaxis, dSP_MC.Cy, 'Marker', marker, 'Color', color, 'MarkerSize', Marksize*2, 'LineWidth', Marksize-1, 'LineStyle', 'none')
    end
end
%% legend
subplot(1, 2, 2)
legend('Halton, J. (1964)', 'CMD', 'non-uniform weight', 'Julier et al. (1995)', 'Julier and Uhlmann (2004)', 'Menegaz et al. (2011)', 'Lerner (2002)', 'Charalampidis and Papavassilopoulos (2011)', 'van der Merwe (2004)')