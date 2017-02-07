% function fh = plotSCTLstatstat(Data,Sim,fh,options)
function fh = plotCz(varargin)

%% Check and assign inputs
if nargin >= 2
    CzSP = varargin{1};
    Cztrue = varargin{2};
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
Cz_abs_max = max([max(abs(CzSP(:))),max(abs(Cztrue(:)))]);

subplot(1,2,1); hold off;

imagesc(CzSP,Cz_abs_max*[-1,1]); hold on;
title('Cz Sigma Point')
xlabel('column'); ylabel('row');
colorbar;
caxis([min(min(min(CzSP,Cz_abs_max))),max(max(max(CzSP,Cz_abs_max)))])
axis square

subplot(1,2,2); hold off;

imagesc(Cztrue,Cz_abs_max*[-1,1]); hold on;
title('Cz sampled')
xlabel('column'); ylabel('row');
colorbar;
caxis([min(min(min(CzSP,Cz_abs_max))),max(max(max(CzSP,Cz_abs_max)))])
axis square

%%
drawnow

end

