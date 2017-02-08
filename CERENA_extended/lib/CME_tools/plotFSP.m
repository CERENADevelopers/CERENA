% plotFSP.m
% plots the solution P of a simulation of a FSP with a FSP state to species
% number mapping stored in 'index'. If the system is high-dimensional,
% a marginalization can be performed such that only dimension 'dim'
% are kept.
%
% USAGE:
% ======
% [fig] = plotFSP(P,t,index)
% [fig] = plotFSP(P,t,index,dim)
% [fig] = plotFSP(P,t,index,dim,options)
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
%   .xlim ... xlim of figure.
%   .ylim ... ylim of figure.
%   .zlim ... zlim of figure.
%
% OUTPUT:
% =======
% fig ... figure handle.
%       
% 27/01/2011 ... Jan Hasenauer

% function [fig] = plotFSP(P,t,index,dim,options)
function [fig] = plotFSP(varargin)

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
options.xlim = [];
options.ylim = [];
options.zlim = [];
if length(dim) == 1
    options.view = [0,90];
else
    options.view = [140,50];
end
options.shading = 'interp';
options.open_figure = 1;
if nargin == 5
    options = setdefault(varargin{5},options);
end
% check xlim
if    (~isempty(options.xlim) || length(options.xlim)) ...
   && (~isempty(options.ylim) || length(options.ylim)) ...
   && (~isempty(options.zlim) || length(options.zlim))
    error('Dimension of xlim, ylim or zlim does not agree.');
end

%% MARGINALIZATION
% Check if marginalization is required
if length(dim) < size(index,2)
    [P,index] = getMarginalization(P,index,dim);
    % Sort
    M = max(index(:));
    w = M.^[size(index,2)-1:-1:0];
    W = zeros(size(index,1),1);
    for i = 1:size(index,1)
        W(i) = (M-index(i,:))*w';
    end
    [~,W] = sort(W,1,'descend');
    index = index(W,:);
    P = P(W,:);
end

% Adapt dimensions of other variables:
options.species = options.species(dim);

%% OPEN FIGURE
if options.open_figure == 1
    fig = figure;
else
    fig = [];
end

%% PLOT FSP
switch size(index,2)
    % 1-dimensional system
    case 1
        % Distigish between length(t) = 1 and length(t) > 1
        if length(t) == 1
            % Define Grid
            x = min(index):max(index);
            % Plot
            plot(x,P);
            % Label
            xlabel(char(options.species(1)));
            ylabel('probability');
            % View
            view([0,90]);
            xlim([min(index),max(index)]);
            ylim([0,multimax(P)]);
            grid on;
        else
%             % Define Grid
%             [X1,X2] = meshgrid(t(:),min(index):max(index));
%             % Plot
%             C = sqrt((P+multimin(nonzeros(P(:)))));
%             surface(X1,X2,P,C);
%             % Label
%             xlabel('time');
%             ylabel(char(options.species(1)));
%             zlabel('probability');
%             % View
%             view(options.view);
%             xlim([min(t),max(t)]);
%             ylim([min(index(:,1)),max(index(:,1))]);
%             zlim([0,multimax(P)]);
%             shading('interp')
%             grid on;

            X = min(index):max(index);
            for i = 1:length(X)
                k = X(i);
                % Vertex, face and color Assignment
                x = t(round(0.5:0.5:length(t)))';
                y = repmat(k + 0.5*[1;-1],length(t),1);
                verts = [x,y];
                faces = bsxfun(@plus,[1,2,4,3],2*[0:length(t)-2]');
                cdata = P(i,round(0.5:0.5:size(P,2)))';
                % Patch visualization
                p = patch('Faces',faces,'Vertices',verts); hold on;
                set(p,'FaceColor','interp','FaceVertexCData',cdata,...
                      'EdgeColor','interp','LineWidth',5);
            end
            colorbar;
            % Label
            xlabel('time');
            ylabel(char(options.species(1)));
            zlabel('probability');
            % View
            view(options.view);
            xlim([min(t),max(t)]);
            ylim([min(index(:,1))-0.5,max(index(:,1))+0.5]);

        end
    % 2-dimensional system
    case 2
        % Determine number of subfigures
        l = ceil(sqrt(length(t)));
        if l*(l-1) >= length(t)
            l1 = l-1;
        else
            l1 = l;
        end
        % Define Grid
        tri = delaunay(index(:,1),index(:,2));
        % Loop: timepoints
        for i = 1:length(t)
            % Open subplot
            subplot(l1,l,i);
            % Plot
	        trisurf(tri,index(:,1),index(:,2),P(:,i));            
            % Label
            xlabel(char(options.species(1)));
            ylabel(char(options.species(2)));
            zlabel('probability');
            title(sprintf('t = %.2f',t(i)));
            % View
            view(options.view);
            xlim([multimin(index(:,1)),multimax(index(:,1))]);
            ylim([multimin(index(:,2)),multimax(index(:,2))]);
            zlim([0,max(P(:,i))])
            shading('interp')
            grid on;
       end
    % n-dimensional system with n >= 3
    otherwise
        % Error message:
        error('This routine only allows for 1- and 2-dimensional data.')    
end

if ~isempty(options.xlim)
    xlim(options.xlim);
end
if ~isempty(options.ylim)
    ylim(options.ylim);
end
if ~isempty(options.zlim)
    zlim(options.zlim);
end

drawnow;