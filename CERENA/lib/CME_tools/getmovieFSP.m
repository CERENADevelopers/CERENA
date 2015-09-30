% getmovieFSP.m
% computes the solution of the CME.
%
% USAGE:
% ======
% [mov] = getmovieFSP(P,t,index,dim,options)
%
% INPUT:
% ======
% P ... n x N matrix, where n is the number of states of the FSP and N is
%       the number of timepoints. Each column provides the state of the FSP
%       at one particular time instance.
% t ... N x 1 vector containing the time instances at which the FSP has
%       been evaluated during the simulation.
% index ... n x m matrix describing the mapping between states of the FSP 
%       and molecule numbers. The i.th row of 'index' contains the molecule
%       numbers associated to the i.th state of the FSP. The system has m
%       species.
% dim ... vector of dimensions which are kept after the marginalization.
%       This variable is assumed to have unique entries.
% options ... options of the algorithm:
%   .species ... names of the species.
%   .view ... point of view (default = [140,50]).
%   .shading ... visualization type (default = 'interp').
%   .open_figure ... open figure (default = 1).
%
% OUTPUT:
% =======
% mov ... movie.
%       
% 27/01/2011 ... Jan Hasenauer

% function [mov] = getmovieFSP(P,t,index,dim,options)
function [mov] = getmovieFSP(varargin)

%% CHECK / ASSIGN INPUTS
% Check required inputs
if nargin >= 3
    if (~isempty(varargin{1}) && ...
        ~isempty(varargin{2}) && ...
        ~isempty(varargin{3}))
        % Assign inputs
        P = varargin{1};
        t = varargin{2};
        index = varargin{3};
    else
        % Error message:
        error('This routine requires three non-empty inputs!');
    end
else
    % Error message:
    error('This routine requires three inputs!');
end

% Optional input
dim = 1:min(size(index,2),2);
if nargin >= 4
    if ~isempty(varargin{4})
        dim = setdefault(varargin{4},dim);
    end
end
% Check 'dim'
if min(size(dim)) > 1
    % Error message:
    error('The variable dim must be a vector.');
end 
% Check that 'dim' has unique entries:
dim = sort(unique(dim));
if (max(dim) > size(index,2)) || (min(dim) < 1)
    % Error message:
    error('Dim exceeds dimensions of index.');
end
if length(dim) >= 3
    % Error message:
    error('The maximal possible dimension of dim is two.');
end

% Check dimensions
if size(P,1) ~= size(index,1)
    % Error message:
    error('Dimension of P and index does not agree!');
end
if size(P,2) ~= length(t)
    % Error message:
    error('Dimension of P and t does not agree!');
end

% Options
for i = 1:size(P,1)
    options.species{i} = sprintf('s_%d',i);
end
options.view = [140,50];
options.shading = 'interp';
options.open_figure = 1;
if nargin == 5
    options = setdefault(varargin{5},options);
end

%% OPEN FIGURE
if options.open_figure == 1
    figure;
end
options.open_figure = 0;

%% GENERATE FSP
% Plot FSP
switch length(dim)
    % 1-dimensional system
    case 1
        % Loop: timepoints
        for i = 1:length(t)
            % Limit
            options.ylim = [0,max(P(:))];
            % Plot
            plotFSP(P(:,i),t(i),index,dim,options);
            % Label
            title(sprintf('t = %.2f',t(i)));
            % View
            view([0,90]);
            % Plot properties
            drawnow;
            hold off;
            % Get frame
            mov(i) = getframe;
       end
    % 2-dimensional system
    case 2
        % Define Grid
        for i = 1:length(t)
            % Limit
            options.zlim = [0,max(P(:))];
            % Plot
            plotFSP(P(:,i),t(i),index,dim,options);
            % Label
            title(sprintf('t = %.2f',t(i)));
            % Plot properties
            drawnow;
            hold off;
            % Get frame
            mov(i) = getframe;
       end
    % n-dimensional system with n >= 3
    otherwise
        % Error message:
        error('This routine only allows for 1- and 2-dimensional data.')    
end