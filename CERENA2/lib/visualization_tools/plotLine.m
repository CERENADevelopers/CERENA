% function plotLine(x,t,options)
function plotLine(varargin)
if nargin >=1
    x = varargin{1};
    if (nargin >= 2) && (~isempty(varargin{2}))
        t = varargin{2};
    else
        t = 1:size(x,1);
    end
else
    error('At least one input argumetn is required!');
end

% Default colors

% Defaut options - plotting
options.lw = 2;
options.lc = [0.0,0.0,0.5];
options.m  = false;
options.ms = 10;
options.msym  = 'o';
options.fs = 16;
options.leg = false;

% Defaut options - labling
options.figure = true;
options.hold = false;
options.fig_title = 'Figure_1';
options.xlabel = 'time';
% options.ylabel = cell(1,size(x,2));
% for i = 1:size(x,2)
%     options.ylabel(i) = ['x'];
% end

% Defaut options - saving
options.save = false;
options.save_path = [pwd,'/figures'];

if nargin >= 3
    options = setdefault(varargin{3},options);
    if nargin > 3
        warning('Only the first three input arguments are used!')
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
        strplot = ['plot(t,x(:,i),''-',msym,''',''color'',lc,''linewidth'',lw,''MarkerSize'',ms)'];
        eval(strplot)
%         plot(t,x(:,i),['''-',msym,''''],'color',lc,'linewidth',lw,'MarkerSize',ms)
    else
        plot(t,x(:,i),'-','color',lc,'linewidth',lw)
    end
    if options.hold
        hold on
    end
    if options.leg && i==1
        leg = legend(options.legstr);
        set(leg,'fontsize',fs)
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