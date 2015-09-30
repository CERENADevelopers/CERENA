% function plotMM(System,x,y,t,options)
% function plotMM(System,options)
function plotMCM(varargin)
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
options.plot_xmcm = true;
options.plot_x = true;
options.plot_y = true;
options.plot_xo = false;
options.plot_yo = false;
options.compare = false;
if nargin >= 2
    options = setdefault(varargin{2},options);
end
%% Reorder
if options.plot_xmcm
    n_y = size(System.CMM.state.stochatic.FSP_index,1);
    n_z = length(System.CMM.state.expectation.state_index);
    n_C = size(System.CMM.state.expectation.C_index,1);
    for iy = 1:n_y
        p_CMM{iy}   = x(:,iy);
        cmu_CMM{iy} = x(:,(n_z+n_C)*(iy-1)+n_y    +[1:n_z]);
        cC_CMM{iy}  = x(:,(n_z+n_C)*(iy-1)+n_y+n_z+[1:n_C]);
    end
    %% Marginal probabilities
    p_x = x(:,1:n_y);
    
    options_p = options;
    options_p.fig_title = {'MCM marginal probabilities'};
    options_p.ylabel = System.CMM.state.sym.p;
    plotLine(p_x,t,options_p);
    %% States
    % Indices of variances
    options_x = options;
    options_x.ylabel = System.CMM.state.expectation.name;
    if System.CMM.order >= 2
        ind_var_x = find((System.CMM.state.expectation.C_index(:,end-1)== System.CMM.state.expectation.C_index(:,end)).*(sum(System.CMM.state.expectation.C_index~=0,2)==2));
    end
    
    for iy=1:n_y
        options_x.fig_title = {['MCM conditional moments - stochastic state y_',num2str(iy)]};
        mean_x = cmu_CMM{iy};
        if System.CMM.order >= 2
            var_x = cC_CMM{iy}(:,ind_var_x);
        else
            var_x = [];
        end
        plotInterval(mean_x,var_x,t,options_x);
    end
end
%% Overall moments of species
if options.plot_x
    % Indices of means and variances
    ind_mean_x = find(sum(System.CMM.state_moments.order>=1,2) == 1);
    if System.CMM.order >= 2
        ind_var_x = find((System.CMM.state_moments.order(:,end-1)== System.CMM.state_moments.order(:,end)).*(sum(System.CMM.state_moments.order~=0,2)==2));
    end
    
    mean_x = mx(:,ind_mean_x);
    if System.CMM.order >= 2
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
        options_x.fig_title = {'MCM species overall moments'};
        plotInterval(mean_x,var_x,t,options_x);
    end
end

%% Outputs
if isfield(System,'output')
    if options.plot_y
        % Indices of means and variances
        ind_mean_y = find(sum(System.CMM.output.order>=1,2) == 1);
        if System.CMM.output_order >= 2
            ind_var_y = find((System.CMM.output.order(:,end-1)== System.CMM.output.order(:,end)).*(sum(System.CMM.output.order~=0,2)==2));
        end
        
        mean_y = y(:,ind_mean_y);
        if System.CMM.output_order >= 2
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
            options_y.fig_title = {'MCM output moments'};
            plotInterval(mean_y,var_y,t,options_y);
        end
    end
end
%% Higher-order moments
if System.CMM.order >= 2
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
            for iy=1:n_y
                for i = 2:xo
                    options_xo.fig_title = {['Conditional moments of order ',num2str(i),' - stochastic state y_',num2str(iy)],ls_species};
                    ind_xo = find((sum(System.CMM.state.expectation.C_index~=0,2)==i));
                    options_xo.ylabel = System.CMM.state.sym.C{iy}(ind_xo);
                    plotLine(cC_CMM{iy}(:,ind_xo),t,options_xo)
                end
            end
        end
    end
end
if System.CMM.order >= 2
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
            for i = 2:xo
                options_xo.fig_title = {['Overall moments of order ',num2str(i)],ls_species};
                ind_xo = find((sum(System.CMM.state_moments.order~=0,2)==i));
                for j = 1:length(ind_xo)
                    options_xo.ylabel(j,:) = {['C_{',num2str(System.CMM.state_moments.order(ind_xo(j),System.CMM.state_moments.order(ind_xo(j),:)~=0)),'}']};
                end
                plotLine(mx(:,ind_xo),t,options_xo)
            end
        end
    end
end
if isfield(System,'output')
    if System.CMM.output_order >= 2
        if options.plot_yo
            if isfield(options,'output_order')
                yo = options.output_order;
                options_yo = options;
                ls_outputs = [];
                for i = 1:System.output.number-1
                    ls_outputs = [ls_outputs,num2str(i),': ',System.output.name{i},',   '];
                end
                ls_outputs = [ls_outputs,num2str(System.output.number),': ',System.output.name{end}];
                options_yo.fs = 12;
                for i = 2:yo
                    options_yo.fig_title = {['Outputs moments of order ',num2str(i)],ls_outputs};
                    ind_yo = find((sum(System.CMM.output.order~=0,2)==i));
                    for j = 1:length(ind_yo)
                        options_yo.ylabel(j,:) = {['C_{',num2str(System.CMM.output.order(ind_yo(j),System.CMM.output.order(ind_yo(j),:)~=0)),'}']};
                    end
                    plotLine(y(:,ind_yo),t,options_yo)
                end
            end
        end
    end
end