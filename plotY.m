% function fh = plotY(SP,fh,options)
function fh = plotY(varargin)

%% Check and assign inputs
if nargin >= 2
    SP = varargin{1};
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
n_y = size(SP.Y,2);
n_t = size(SP.Y,1);
nc = n_y;
nr = n_y;

%% Visualization: Data and Simulation
% Loop: measurands
for it = 1:n_t
    figure
    for iy = 1:n_y
        for jy = iy:n_y
            % Data and simulation
            subplot(nr,nc,n_y*(jy-1)+iy); hold off;
            
            if(iy == jy)
                histogram(squeeze(SP.Y_true(it,iy,:)),'Normalization','pdf');
%                 bar(xx,hh/trapz(xx,hh))
                hold on
                xl = get(gca,'XLim');
                xx = linspace(xl(1),xl(2),100);
                xlabel(['state ' num2str(iy)]);
                plot(xx,normpdf(xx,SP.my(it,iy),sqrt(SP.Cy(it,iy,iy))),'r-');
            else
                plot(squeeze(SP.Y_true(it,iy,:)),squeeze(SP.Y_true(it,jy,:)),'bx');
                hold on
                X = getEllipse([SP.my(it,iy),SP.my(it,jy)],squeeze(SP.Cy(it,[iy,jy],[iy,jy])),2);
                plot(SP.my(it,iy),SP.my(it,jy),'ro','linewidth',3);
                plot(X(1,:),X(2,:),'r-','linewidth',2);
                ylabel(['state ' num2str(jy)]);
                xlabel(['state ' num2str(iy)]);
            end
            
        end
    end
end

%%
drawnow

end

