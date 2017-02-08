% bar_vw.m plots a bar diagram width variable bar width.
%
% USAGE:
% ======
% fig = bar_vw(h,x)
% fig = bar_vw(h,x,options)
%
% INPUT:
% ======
% h ... height of individual bars.
% x ... bounds for bars. The lower bound of bar i is x(i) and the upper
%       bound of bar i is x(i+1).
% options ... options of the routine:
%   .color ... color of the bars.
%   .linestyle ... linestyle of the edge of the bar. (default = 'k')
%                       Note: If no line is wanted enter 'none'. 
%   .open_figure ... = 1 (default) => new figure is opened 
%   .normalize ... = 1 (default) => histogram is normalized to have
%                       integral one.
%   .reduce_width ... the width of each bar is reduced by this value.
%                       (default = 0).
% 
% OUTPUT:
% =======
% fig ... figure handle.
%
% 16/01/2011 - Jan Hasenauer

% function fig = bar_vw(h,x,options)
function [varargout] = bar_vw(varargin)

%% CHECK / ASSIGN INPUTS
% Data
if nargin >= 2
	if isempty(varargin{1}) || isempty(varargin{2})
		error('First two inputs must be non-empty!');
	else
		h = varargin{1};
		x = varargin{2};
	end
else
	error('Not enough inputs.');
end

% Check inputs
if length(x) ~= (length(h)+1)
    error('Length second input must be length of first input + 1.');
end

% Set default options
options.color = 'b';
options.linestyle = '-';
options.open_figure = 1;
options.normalize = 1;
options.reduce_width = 0;
if nargin >= 3
	options = setdefault(varargin{3},options);
end

%% OPEN FIGURE
if options.open_figure == 1
    fig = figure;
else
    fig = [];
end

%% NORMALIZE HISTORGAM
if options.normalize == 1
    s = diff(x(:))'*h(:);
    h = h/s;
end

%% PLOT HISTOGRAM
for i = 1:length(h)
    fill([x(i)+options.reduce_width;x(i)+options.reduce_width;...
          x(i+1)-options.reduce_width;x(i+1)-options.reduce_width],...
          [0;h(i);h(i);0],options.color,...
        'LineStyle',options.linestyle); hold on;
end

%% ASSIGN OUTPUT
switch nargout
    case 0
        varargout{1} = [];
    case 1
        varargout{1} = fig;
    otherwise
        error('Too many output arguments.')
end
