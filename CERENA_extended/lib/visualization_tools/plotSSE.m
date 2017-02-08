% function plotMM(System,x,mx,y,t,options)
% function plotMM(System,options)
function plotSSE(varargin)
if nargin >= 1
    System = varargin{1};
else
    error('At least one input argument is required!')
end
x = System.sol.x;
mx = System.sol.mx;
y = System.sol.y;
t = System.sol.t;

options.save = false;
options.plot_xsse = true;
options.plot_x = true;
options.plot_y = true;
options.plot_xo = false;
options.plot_yo = false;
options.compare = false;
if nargin >= 2
    options = setdefault(varargin{2},options);
end
%% States
if options.plot_xsse
    options_xsse = options;
    options_xsse.fig_title = {'SSE states'};
    options_xsse.ylabel = System.states;
    plotLine(x,t,options_xsse);
end
%% Species moments
if options.plot_x
    ind_mean_x = find(sum(System.state.order>=1,2) == 1);
    if ~strcmp(System.name,'RRE')
        ind_var_x = find((System.state.order(:,end-1)== System.state.order(:,end)).*(sum(System.state.order~=0,2)==2));
    end
    mean_x = mx(:,ind_mean_x);
    if ~strcmp(System.name,'RRE')
        var_x = mx(:,ind_var_x);
    else
        var_x = [];
    end
    
    options_x = options;
    options_x.ylabel = System.state.name;
    if options.compare
        options_x.fig_title = {options.fig_title_x{1}};
        options_x.fh = options.fhx{1};
        plotLine(mean_x,t,options_x);
        options_x.fig_title = {options.fig_title_x{2}};
        options_x.fh = options.fhx{2};
        plotLine(var_x,t,options_x);
    else
        options_x.fig_title = {'SSE species moments'};
        plotInterval(mean_x,var_x,t,options_x);
    end
end
%% Outputs
if isfield(System,'output')
    if options.plot_y
        ind_mean_y = find(sum(System.output.order>=1,2) == 1);
        if ~strcmp(System.name,'RRE')
            ind_var_y = find((System.output.order(:,end-1)== System.output.order(:,end)).*(sum(System.output.order~=0,2)==2));
        end
        mean_y = y(:,ind_mean_y);
        if ~strcmp(System.name,'RRE')
            var_y = y(:,ind_var_y);
        else
            var_y = [];
        end
        
        options_y = options;
        options_y.ylabel = System.output.name;
        if options.compare
            options_y.fig_title = {options.fig_title_y{1}};
            options_y.fh = options.fhy{1};
            plotLine(mean_y,t,options_y);
            options_y.fig_title = {options.fig_title_y{2}};
            options_y.fh = options.fhy{2};
            plotLine(var_y,t,options_y);
        else
            options_y.fig_title = {'SSE output moments'};
            plotInterval(mean_y,var_y,t,options_y);
        end
    end
end
%% Higher-order moments
if ~strcmp(System.name,'RRE')
    if options.plot_xo
        options_xo = options;
        ls_species = [];
        for i = 1:System.state.number-1
            ls_species = [ls_species,num2str(i),': ',System.state.name{i},',   '];
        end
        ls_species = [ls_species,num2str(System.state.number),': ',System.state.name{end}];
        options_xo.fs = 12;
        if any([strcmp(System.name,'LNA'),strcmp(System.name,'EMRE')])
            xo = 2;
        elseif strcmp(System.name,'IOS')
            xo = 3;
        end
        for i = 2:xo
            options_xo.fig_title = {['Species moments of order ',num2str(i)],ls_species};
            ind_xo = find((sum(System.state.order~=0,2)==i));
            for j = 1:length(ind_xo)
                options_xo.ylabel(j,:) = {['C_{',num2str(System.state.order(ind_xo(j),System.state.order(ind_xo(j),:)~=0)),'}']};
            end
            plotLine(mx(:,ind_xo),t,options_xo)
        end
    end
end

if ~strcmp(System.name,'RRE')
    if isfield(System,'output')
        if options.plot_yo
            options_yo = options;
            ls_outputs = [];
            for i = 1:System.output.number-1
                ls_outputs = [ls_outputs,num2str(i),': ',System.output.name{i},',   '];
            end
            ls_outputs = [ls_outputs,num2str(System.output.number),': ',System.output.name{end}];
            options_yo.fs = 12;
            if any([strcmp(System.name,'LNA'),strcmp(System.name,'EMRE')])
                yo = 2;
            elseif strcmp(System.name,'IOS')
                yo = 3;
            end
            for i = 2:yo
                options_yo.fig_title = {['Outputs moments of order ',num2str(i)],ls_outputs};
                ind_yo = find((sum(System.output.order~=0,2)==i));
                for j = 1:length(ind_yo)
                    options_yo.ylabel(j,:) = {['C_{',num2str(System.output.order(ind_yo(j),System.output.order(ind_yo(j),:)~=0)),'}']};
                end
                plotLine(y(:,ind_yo),t,options_yo)
            end
        end
    end
end
