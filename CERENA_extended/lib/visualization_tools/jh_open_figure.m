% jh_open_figure.m
% can be used to open a formatted or unformatted figure.
%
% USAGE:
% ======
% fig_h = jh_open_figure(options)
%
% INPUT:
% ======
% options ... options for figure:
%   .type ... type of figure:
%       = 'formatted' ... formatted figure is opened.
%       = 'unformatted' (default) ... formatted figure is opened.
%   .size ... size of formatted figure (default = [10,8]).
%   .offset ... offset of formatted figure (default = [1.8,1.4]).
%
% OUPUT:
% ======
% fig_h ... figure handle.
%
% 17/02/2011 - Jan Hasenauer

% function fig_h = jh_open_figure(options)
function fig_h = jh_open_figure(varargin)

%% CHECK / ASSIGN INPUT
if nargin == 1
    options.type = 'unformatted';
    options.size = [10,8];
    options.offset = [1.8,1.4];
    options = setdefault(varargin{1},options);
else
    error('This function requires precisely one input')
end

%% CHOOSE AND OPEN FIGURE TYPE
switch options.type
    % Formatted figure
    case 'formatted'
        [fig_h, axes_h] = jh_figures(options.size,[1,1],options.offset);
        axes(axes_h);
    % Unformatted figure
    otherwise
        fig_h = figure;
end