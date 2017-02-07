% function fh = plotmz(t,mzSP,mztrue,fh,options)
function fh = plotmz(varargin)

%% Check and assign inputs
if nargin >= 2
    mzSP = varargin{1};
    mztrue = varargin{2};    
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
options.data.area_col = [0.7,0.7,1];
options.data.ls = '-';
options.data.mean_lw = 2;
options.data.bound_lw = 1;
options.sim.col = 'r';
options.sim.ls = '--';
options.sim.mean_lw = 2;
options.sim.bound_lw = 1;
options.error.col = 'b';
options.error.ls = '-';
options.error.lw = 1;
if nargin == 4
    options = setdefault(varargin{4},options);
end

%% Subplot dimensions
n_y = 1;
n_t = size(mztrue,1);
nr = 1;
nc = 2;

%% Visualization: Data and Simulation
% Data and simulation
subplot(nr,nc,1); hold off;

lh(1) = plot(1:n_t,mzSP(:,1),'-',...
    'linewidth',options.data.mean_lw,...
    'linestyle',options.data.ls,...
    'color',options.data.col); hold on;

lh(2) = plot(1:n_t,mztrue(:,1),'-',...
    'linewidth',options.sim.mean_lw,...
    'linestyle',options.sim.ls,...
    'color',options.sim.col); hold on;

xlabel('time'); ylabel('state');
xlim([1,n_t]);
legend(lh,{'mz SP','mz sampled'});


subplot(nr,nc,2); hold off;
plot(1:n_t,mzSP(:,1)-mztrue(:,1),'-',...
    'linewidth',options.error.lw,...
    'linestyle',options.error.ls,...
    'color',options.error.col); hold on;

xlabel('time'); ylabel('error');
xlim([1,n_t]);

%%
drawnow

end

