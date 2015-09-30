% function plotCompare(System, options)
function plotCompare(varargin)
if nargin >=1
    System = varargin{1};
else
    error('At least one input argument is required!')
end
n_System = length(System);
for i = 1:n_System
    method{i} = System{i}.name;
end

options.save = false;
options.m = false;
options.plot_x = 1;
options.plot_y = 1;
if nargin >= 2
    options = setdefault(varargin{2},options);
end
%% Line colors
lc  =  [0.8,0.0,0.0;
    0.0,0.0,0.8;
    0.0,0.8,0.0;
    0.5,0.3,0.9;
    0.7,0.8,1.0;
    1.0,0.6,0.6;
    0.6,0.6,0.6;
    0.4,0.3,0.7;
    0.7,0.6,1.0;
    1.0,0.2,0.2;
    0.5,0.0,0.0;
    0.0,0.0,0.5;
    0.0,0.5,0.0;
    0.3,0.0,0.0;
    0.0,0.0,0.3;
    0.0,0.3,0.0;
    0.4,0.4,0.4;
    1.0,0.4,0.4;
    0.0,0.0,1.0;
    0.6,0.6,1.0;
    0.4,0.4,1.0;
    0.2,0.2,1.0;
    1.0,0.0,0.0;
    0.5,0.5,0.5;
    0.6,0.6,0.6;
    0.9,0.9,0.9];
msym = {'*','s','o','^','v','+','.'};
%% Mean and variance of species
optionsAll = options;
if options.plot_x
    h_mean_x = figure;
    h_var_x = figure;
    optionsAll.plot_x = true;
    optionsAll.fhx = {h_mean_x,h_var_x};
    optionsAll.fig_title_x = {'Mean of species','Variance of species'};
else
    optionsAll.plot_x = false;
end
if options.plot_y
    h_mean_y = figure;
    h_var_y = figure;
    optionsAll.plot_y = true;
    optionsAll.fhy = {h_mean_y,h_var_y};
    optionsAll.fig_title_y = {'Mean of output variables','Variance of output variables'};
else
    optionsAll.plot_y = false;
end
optionsAll.plot_xo = false;
optionsAll.plot_yo = false;
optionsAll.figure = false;
optionsAll.hold = true;
optionsAll.save = options.save;
optionsAll.compare = 1;
optionsAll.leg = true;
for i = 1:n_System
    System_i = System{i};
    method_i = method{i};
    optionsAll.legstr{i} = method_i;
    if isfield(options,'lc') && size(options.lc,1) >= i
        optionsAll.lc = options.lc(i,:);
    else
        optionsAll.lc = lc(i,:);
    end
    if options.m
        if isfield(options,'msym') && length(options.msym) >= i
            optionsAll.msym = options.msym{i};
        else
            optionsAll.msym = msym{i};
        end
    end
    switch method_i
        case 'MM'
            options_MM = optionsAll;
%             plotMM(System_i,System_i.sol.x,System_i.sol.y,System_i.sol.t,options_MM)
            plotMM(System_i,options_MM)
        case 'MCM'
            options_MCM = optionsAll;
            options_MCM.plot_xmcm = false;
%             plotMCM(System_i,System_i.sol.x,System_i.sol.mx,System_i.sol.y,System_i.sol.t,options_MCM)
            plotMCM(System_i,options_MCM)
        case {'RRE','EMRE','LNA','IOS'}
            options_SSE = optionsAll;
            options_SSE.plot_xsse = false;
%             plotSSE(System_i,System_i.sol.x,System_i.sol.mx,System_i.sol.y,System_i.sol.t,options_SSE)
            plotSSE(System_i,options_SSE)
        case 'FSP'
            options_FSP = optionsAll;
            options_FSP.plot_xfsp = false;
%             plotFSP_ACME(System_i,System_i.sol.x,[],System_i.sol.t,options_FSP)
            plotFSP_ACME(System_i,options_FSP)
        case 'SSA'
            options_SSA = optionsAll;
            options_SSA.plot_xssa = false;
%             plotSSA(System_i,System_i.sol.x,System_i.sol.y,System_i.sol.mean_x,System_i.sol.var_x,System_i.sol.mean_y,System_i.sol.var_y,System_i.sol.t,options_SSA)
            plotSSA(System_i,options_SSA)
    end
end