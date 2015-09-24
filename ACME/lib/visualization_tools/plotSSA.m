% function plotSSA(System,x,y,mean_x,var_x,mean_y,var_y,t,options)
% function plotSSA(System,options)
function plotSSA(varargin)
if nargin >= 1
    System = varargin{1};
else
    error('At least one input argument is required!')
end
x = System.sol.x;
y = System.sol.y;
mean_x = System.sol.mean_x;
var_x = System.sol.var_x;
mean_y = System.sol.mean_y;
var_y = System.sol.var_y;
t = System.sol.t;

options.save = false;
options.mean_realization = true;
options.plot_xssa = true;
options.plot_x = true;
options.plot_y = true;
options.plot_xo = true;
options.plot_yo = true;
options.compare = false;
if nargin >= 2
    options = setdefault(varargin{2},options);
end
%% Realization
if options.plot_xssa
    if isfield(options,'n_realization')
        nr = options.n_realization;
    else
        nr = min(10,size(x,3));
    end
    ir = randperm(size(x,3),nr);
    options_r = options;
    options_r.ylabel = System.state.name;
    options_r.figure = false;
    options_r.hold = true;
    options_r.lc = [0.0,0.7,0.0];
    options_r.lw = 0.25;
    options_r.fig_title = {'SSA - realizations'};
    figure;
    for i=1:nr
        plotLine(x(:,:,ir(i)),t,options_r)
    end
    if options_r.mean_realization
        options_r.lw = 1;
        options_r.lc = [0.0,0.0,0.5];
        plotInterval(mean_x,var_x,t,options_r)
    end
end
%% States
if options.plot_x
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
        options_x.fig_title = {'SSA - Moments of species'};
        plotInterval(mean_x,var_x,t,options_x);
    end
end
%% Outputs
if options.plot_y
    % Indices of means and variances
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
        options_y.fig_title = {'SSA - output moments'};
        plotInterval(mean_y,var_y,t,options_y);
    end
end
%% Higher-order moments
if options.plot_xo
    if isfield(options,'state_order')
        ls_species = [];
        for i = 1:System.state.number-1
            ls_species = [ls_species,num2str(i),': ',System.state.name{i},',   '];
        end
        ls_species = [ls_species,num2str(System.state.number),': ',System.state.name{end}];
        xo = options.state_order;
        options_xo = options;
        options_xo.fs = 12;
        for i = 2:xo
            options_xo.fig_title = {['SSA - Species moments of order ',num2str(i)],ls_species};
            [ind_xo,~] = getMomentIndexSet(System.state.number,i);
            ind_xo = ind_xo(sum(ind_xo~=0,2)==i,:);
            n_M_xo = size(ind_xo,1);
            M_xo = zeros(size(x,1),n_M_xo);
            for j = 1:n_M_xo
                I_x = ind_xo(j,:);
                alpha_x = convertI2alpha(I_x,System.state.number);
                M_xo(:,j) = mean(prod(bsxfun(@power,bsxfun(@minus,x,mean_x),alpha_x),2),3);
                options_xo.ylabel(j,:) = {['C_{',num2str(I_x),'}']};
            end
            plotLine(M_xo,t,options_xo)
        end
    end
end
if options.plot_yo
    if isfield(options,'output_order')
        ls_outputs = [];
        for i = 1:System.output.number-1
            ls_outputs = [ls_outputs,num2str(i),': ',System.output.name{i},',   '];
        end
        ls_outputs = [ls_outputs,num2str(System.output.number),': ',System.output.name{end}];
        yo = options.output_order;
        options_yo = options;
        options_yo.fs = 12;
        for i = 2:yo
            options_yo.fig_title = {['SSA - Outputs moments of order ',num2str(i)],ls_outputs};
            [ind_yo,~] = getMomentIndexSet(System.output.number,i);
            ind_yo = ind_yo(sum(ind_yo~=0,2)==i,:);
            n_M_yo = size(ind_yo,1);
            M_yo = zeros(size(y,1),n_M_yo);
            for j = 1:n_M_yo
                I_y = ind_yo(j,:);
                alpha_y = convertI2alpha(I_y,System.output.number);
                M_yo(:,j) = mean(prod(bsxfun(@power,bsxfun(@minus,y,mean_y),alpha_y),2),3);
                options_yo.ylabel(j,:) = {['C_{',num2str(I_y),'}']};
            end
            plotLine(M_yo,t,options_yo)
        end
    end
end