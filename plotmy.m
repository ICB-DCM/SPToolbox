% function fh = plotPA(t,mySP,mytrue,fh,options)
function fh = plotmy(varargin)

%% Check and assign inputs
if nargin >= 2
    mySP = varargin{1};
    mytrue = varargin{2};
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
n_y = size(mytrue,2);
n_t = size(mytrue,1);
nc = 2;
nr = n_y;

%% Visualization: Data and Simulation
% Loop: measurands
for j = 1:n_y
    % Data and simulation
    subplot(nr,nc,2*(j-1)+1); hold off;
    
    lh(1) = plot(1:n_t,mySP(:,j),'-',...
        'linewidth',options.data.mean_lw,...
        'linestyle',options.data.ls,...
        'color',options.data.col); hold on;
    
    lh(2) = plot(1:n_t,mytrue(:,j),'-',...
        'linewidth',options.sim.mean_lw,...
        'linestyle',options.sim.ls,...
        'color',options.sim.col); hold on;
    
    xlabel('time'); ylabel('state');
    xlim([1,max(n_t,2)]);
    if j == 1
        legend(lh,{'my SP','my sampled'});
    end
    
    subplot(nr,nc,2*(j-1)+2); hold off;
    plot(1:n_t,mySP(:,j)-mytrue(:,j),'-',...
        'linewidth',options.error.lw,...
        'linestyle',options.error.ls,...
        'color',options.error.col); hold on;
    
    xlabel('time'); ylabel('error');
    xlim([1,max(n_t,2)]);
end

%%
drawnow

end

