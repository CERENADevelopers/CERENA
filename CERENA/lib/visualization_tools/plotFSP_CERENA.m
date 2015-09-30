% function plotFSP(System,x,y,t,options)
% function plotFSP(System,options)
function plotFSP_ACME(varargin)
if nargin >= 1
    System = varargin{1};
else
    error('At least one input arguments is required!')
end

x = System.sol.x;
if isfield(System.sol,'Mx')
    M_FSP = System.sol.Mx;
else
    M_FSP = [];
end
if isfield(System.sol,'My')
    y = System.sol.My;
else
    y = [];
end
t = System.sol.t;

options.save = false;
options.plot_xfsp = true;
options.plot_x = true;
options.plot_y = true;
options.plot_xo = false;
options.plot_yo = false;
options.calculate_x = false;
options.calculate_y = false;
options.compare = false;
if nargin >= 2
    options = setdefault(varargin{2},options);
end
%% Probability distribution over FSP states
if options.plot_xfsp
    options.species = System.state.name;
    for i =1:System.state.number
        plotFSP(x,t,System.index,i,options)
    end
end

%% States
if options.plot_x
    if isempty(M_FSP) || options.calculate_x
        [M_FSP, M_ind_FSP] = getMomentsFSP_centered(x,System.index,2);
    else
%         M_FSP = System.state.moment;
        M_ind_FSP = System.state.order;
    end
    ind_mean_x = find(sum(M_ind_FSP>=1,2) == 1);
    ind_var_x = find((M_ind_FSP(:,end-1)== M_ind_FSP(:,end)).*(sum(M_ind_FSP~=0,2)==2));
    mean_x = M_FSP(:,ind_mean_x);
    var_x = M_FSP(:,ind_var_x);
    % Indices of variances
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
        options_x.fig_title = {'FSP - Moments of species'};
        plotInterval(mean_x,var_x,t,options_x);
    end
end
%% Outputs
if options.plot_y
    % Indices of means and variances
    if isempty(y)
        if options.calculate_y
            options.output.calculate = 1;
            options.moment_order = 2;
            options.moment_order_output = 2;
            [fun_y,H,My,Iy,n_Iy] = getOutputMoments(System,M_ind_FSP,M_FSP,[],[],[],[],[],options);
            y = fun_y(M_FSP,System.sol.theta);
            System.output.order = Iy;
        end
    end
    ind_mean_y = find(sum(System.output.order>=1,2) == 1);
    ind_var_y = find((System.output.order(:,end-1)== System.output.order(:,end)).*(sum(System.output.order~=0,2)==2));
    mean_y = y(:,ind_mean_y);
    var_y = y(:,ind_var_y);
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
        options_y.fig_title = {'FSP - output moments'};
        plotInterval(mean_y,var_y,t,options_y);
    end
end
%% Higher-order moments
if options.plot_xo
    if isfield(options,'state_order')
        xo = options.state_order;
        options_xo = options;
        ls_species = [];
        for i = 1:System.state.number-1
            ls_species = [ls_species,num2str(i),': ',System.state.name{i},',   '];
        end
        ls_species = [ls_species,num2str(System.state.number),': ',System.state.name{end}];
        options_xo.fs = 12;
        if options.calculate_xo
            [M_FSP, M_ind_FSP] = getMomentsFSP_centered(x,System.index,xo);
        end
        for i = 2:xo
            options_xo.fig_title = {['FSP - Species moments of order ',num2str(i)],ls_species};
            ind_xo = find((sum(M_ind_FSP~=0,2)==i));
            for j = 1:length(ind_xo)
                options_xo.ylabel(j,:) = {['C_{',num2str(M_ind_FSP(ind_xo(j),M_ind_FSP(ind_xo(j),:)~=0)),'}']};
            end
            plotLine(M_FSP(:,ind_xo),t,options_xo)
        end
    end
end
if options.plot_yo
    if isfield(options,'output_order')
        yo = options.output_order;
        options_yo = options;
        options_yo.fs = 12;
        ls_outputs = [];
        for i = 1:System.output.number-1
            ls_outputs = [ls_outputs,num2str(i),': ',System.output.name{i},',   '];
        end
        ls_outputs = [ls_outputs,num2str(System.output.number),': ',System.output.name{end}];
        if options.calculate_yo
            options.output.calculate = 1;
            options.moment_order = xo;
            options.moment_order_output = yo;
            [fun_y,H,My,Iy,n_Iy] = getOutputMoments(System,M_ind_FSP,M_FSP,[],[],[],[],[],options);
            y = fun_y(M_FSP,System.sol.theta);
            System.output.order = Iy;
        end
        for i = 2:yo
            options_yo.fig_title = {['FSP - Outputs moments of order ',num2str(i)],ls_outputs};
            ind_yo = find((sum(System.output.order~=0,2)==i));
            for j = 1:length(ind_yo)
                options_yo.ylabel(j,:) = {['C_{',num2str(System.output.order(ind_yo(j),System.output.order(ind_yo(j),:)~=0)),'}']};
            end
            plotLine(y(:,ind_yo),t,options_yo)
        end
    end
end