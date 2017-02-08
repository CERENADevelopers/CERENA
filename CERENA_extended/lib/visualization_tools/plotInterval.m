% function plotInterval(x,varx,t,options)
function plotInterval(varargin)
if nargin >=2
    x = varargin{1};
    varx = varargin{2};
    if (nargin >= 3) && (~isempty(varargin{3}))
        t = varargin{3};
    else
        t = 1:size(x,1);
    end
else
    error('At least two input arguments are required!');
end

% Defaut options - plotting
options.lw = 2;
options.lc = [0.0,0.0,0.5];
options.m  = false;
options.ms = 10;
options.msym  = 'o';
options.fs = 16;

% Defaut options - labling
options.figure = true;
options.hold = false;
options.fig_title = 'Figure 1';
options.xlabel = 'time';
options.ylabel = cell(1,size(x,2));
for i = 1:size(x,2)
    options.ylabel{i} = ['x_',num2str(i)];
end

% Defaut options - saving
options.save = false;
options.save_path = [pwd,'/figures'];

if nargin >= 4
    options = setdefault(varargin{4},options);
    if nargin > 4
        warning('Only the first four input arguments are used!')
    end
end
if ~isfield(options,'save_name') || isempty(options.save_name)
    options.save_name = options.fig_title{1};
end

lw = options.lw;
lc = options.lc;
ms = options.ms;
msym = options.msym;
fs = options.fs;

if options.figure
    figure
elseif isfield(options,'fh')
    figure(options.fh)
end
if options.hold
    hold on
end

nSubFig = size(x,2);
nCols = ceil(sqrt(nSubFig));
nRows = ceil(nSubFig/nCols);
for i = 1:nSubFig
    subplot(nRows,nCols,i)
    if options.m
        plot(t,x(:,i),['''-',msym,''''],'color',lc,'linewidth',lw,'MarkerSize',ms); hold on
    else
        plot(t,x(:,i),'-','color',lc,'linewidth',lw); hold on
    end
    if ~isempty(varx)
        plot(t,x(:,i)+sqrt(varx(:,i)),'--','color',lc,'linewidth',lw)
        plot(t,x(:,i)-sqrt(varx(:,i)),'--','color',lc,'linewidth',lw)
    end
    xlabel(strrep(char(options.xlabel),'_','\_'),'fontsize',fs)
    ylabel(strrep(char(options.ylabel(i)),'_','\_'),'fontsize',fs)
    set(gca,'fontsize',fs)
end
set(gcf,'nextPlot','add')
axes;
ht = title(options.fig_title,'fontsize',fs);
ht = strrep(ht,'_','\_');
set(gca,'Visible','off')
set(ht,'Visible','on')
if options.save
    if ~exist(options.save_path,'dir')
        mkdir(options.save_path)
    end
    savefig([options.save_path,'/',options.save_name])
    print('-depsc2','-r1000',[options.save_path,'/',options.save_name])
end