% getMarginalization.m
% computes the marginalization of states of an FSP along one or more
% dimensions. This function does not assume anythin about the structure of
% the FSP.
%
% USAGE:
% ======
% [P_red,index_red] = getMarginalization(P,index,dim)
%
% INPUT:
% ======
% P ... n x N matrix, where n is the number of states of the FSP and N is
%       the number of timepoints. Each column provides the state of the FSP
%       at one particular time instance.
% index ... n x m matrix describing the mapping between states of the FSP 
%       and molecule numbers. The i.th row of 'index' contains the molecule
%       numbers associated to the i.th state of the FSP. The system has m
%       species.
% dim ... vector of dimensions which are kept after the marginalization.
%       This variable is assumed to have unique entries.
%
% OUTPUT:
% =======
% P_red ... marginalized probability distribution.
% index_red ... reduced index matrix only containing the dimensions dim
%       over which no marginalization is performed.
%       
% 27/01/2011 ... Jan Hasenauer

% function [P_red,index_red] = getMarginalization(P,index,dim)
function [P_red,index_red] = getMarginalization(varargin)

%% CHECK / ASSIGN INPUTS
% Check required inputs
if nargin >= 2
    if (~isempty(varargin{1}) && ...
        ~isempty(varargin{2}))
        % Assign inputs
        P = varargin{1};
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
dim = 1;
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

% Check dimensions
if size(P,1) ~= size(index,1)
    % Error message:
    error('Dimensions of P and index do not agree!');
end

%% PERMFORM MARGINALIZATION
% Check if marginalization is required
if length(dim) < size(index,2)
    % Define reduced index set:
    index_red = zeros(size(index,1),length(dim));
    % Ensure that index_red is unique
    l = 1:size(index,1);
    i = 1;
    while ~isempty(l)
        % Compare index_red(i) to all elements of index_red and find
        % those which are identical:
        z = (index(l,dim(1)) == index(l(1),dim(1)));
        for j = 2:length(dim)
            z = z .* (index(l,dim(j)) == index(l(1),dim(j)));
        end
        z = find(z);
        % Assign the mapping and reduce l:
        index_red(i,:) = index(l(1),dim);
        mapping{i} = l(z);
        l = setdiff(l,l(z));
        % Increase counter
        i = i + 1;
    end
    % Truncate index_red:
    index_red = index_red(1:length(mapping),:);
    
    % Perform marginalization
    P_red = zeros(size(index_red,1),size(P,2)); 
    for i = 1:size(index_red,1)
        P_red(i,:) = sum(P(mapping{i},:),1);
    end
    
else
    P_red = P;
    index_red = index;
end
