% plotFSPindex.m
% plots the indexset 'index' obtained for instance via a noise-free flow
% cyometry measurement. 'index' does here not contain the mapping between
% states of the FSP and phyical species number, but a colection of observed
% species numbers. If the system is high-dimensional, a marginalization 
% can be performed such that only dimension 'dim' are kept.
%
% USAGE:
% ======
% [fig] = plotFSPindex(t,index,dim,options)
%
% INPUT:
% ======
% t ... N x 1 vector containing the time instances at which the FSP has
%       been evaluated during the simulation.
% index{i} ... n x m matrix describing the mapping between states of the FSP 
%       and molecule numbers. The i.th row of 'index' contains the molecule
%       numbers associated to the i.th data point. The system has m
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
% fig ... figure handle.
%       
% 27/01/2011 ... Jan Hasenauer

% function [fig] = plotFSPindx(t,index,dim,options)
function [fig] = plotFSPindx(varargin)

%% CHECK / ASSIGN INPUTS
% Check required inputs
if nargin >= 2
    if (~isempty(varargin{1}) && ...
        ~isempty(varargin{2}))
        % Assign inputs
        t = varargin{1};
        index = varargin{2};
    else
        % Error message:
        error('This routine requires two non-empty inputs!');
    end
else
    % Error message:
    error('This routine requires two inputs!');
end

% Optional input
dim = 1:min(size(index{1},2),2);
if nargin >= 3
    if ~isempty(varargin{3})
        dim = setdefault(varargin{3},dim);
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

% Options
options.open_figure = 1;
if nargin == 4
    options = setdefault(varargin{4},options);
end

%% GENERATE STATE INDEX SET:
% Lower and upper bounds:
xmin =  inf*ones(size(index{1},2),1);
xmax = -inf*ones(size(index{1},2),1);
for i = 1:length(t)
    xmin = min((min(index{i},[],1))',xmin);
    xmax = max((max(index{i},[],1))',xmax);
end
% Initialize index matrix:
nx = size(index{1},2);
max_index = prod(xmax-xmin+1);
index_full = zeros(max_index,nx);
% Initialize state
ind = xmin;
% Construct index matrix
for i = 1:prod(xmax-xmin+1)
    % Assign index
    index_full(i,:) = ind';
    % Update current index:
    j = min(setdiff(1:nx,find(ind == xmax))); % index which is updated
    ind(j) = ind(j) + 1; 
    ind(1:j-1) = xmin(1:j-1); % reset of all previous indexes
end

%% CONSTRUCTION OF P
% Initialization
P = zeros(length(t),size(index_full,1));
% Loop: time instances
for i = 1:length(t)
    % Loop: elements of index set 'index'
    for j = 1:size(index{i},1)
        % Determine match of 'index' and 'index_full'
        ind = (index{i}(j,1) == index_full(:,1));
        for k = 2:size(index{i},2)
            ind = ind .* (index{i}(j,k) == index_full(:,k));
        end
        ind_j = find(ind);
        % Update P
        P(i,ind_j) = P(i,ind_j) + 1;
    end
    % Normaliztion
    P(i,:) = P(i,:)./sum(P(i,:));    
end

%% CALL plotFSP.m TO PLOT P
[fig] = plotFSP(P,t,index_full,dim,options);
