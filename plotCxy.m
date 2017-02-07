% function fh = plotSCTLstatstat(Data,Sim,fh,options)
function fh = plotCxy(varargin)

%% Check and assign inputs
if nargin >= 2
    CxySP = varargin{1};
    Cxytrue = varargin{2};
else
    error('Not enough inputs.')
end

% Figure handel
if nargin >= 3
    if ~isempty(varargin{3})
        fh = figure(varargin{3});
    else
        fh = figure;
    end
else
    fh = figure;
end

% Options
options.data.col = 'b';
options.data.ls = '-';
options.data.lw = 2;
options.sim.col = 'r';
options.sim.ls = '-';
options.sim.lw = 2;
options.error.col = 'b';
options.error.ls = '-';
options.error.lw = 1;
if nargin == 4
    options = setdefault(varargin{4},options);
end

%% Visualization: Data and Simulation
% Number of outputs
ny = 1;

% Covariances - Data
Cxy_abs_max = max([max(abs(CxySP(:))),max(abs(Cxytrue(:)))]);

subplot(1,2,1); hold off;
title('Cxy SP')
imagesc(CxySP,Cxy_abs_max*[-1,1]); hold on;
xlabel('column'); ylabel('row');
colorbar;
caxis([min(min(min(CxySP,Cxy_abs_max))),max(max(max(Data.SCTLstat.Cz,Cxy_abs_max)))])
axis square

subplot(1,2,2); hold off;
title('Cxy sampled')
imagesc(Cxytrue,Cxy_abs_max*[-1,1]); hold on;
xlabel('column'); ylabel('row');
colorbar;
caxis([min(min(min(CxySP,Cxy_abs_max))),max(max(max(Data.SCTLstat.Cz,Cxy_abs_max)))])
axis square

%%
drawnow

end

