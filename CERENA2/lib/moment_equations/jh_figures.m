function [handle, axes_h] = jh_figures(Size,n,Axpos,FontSize)
% JH_FIGURE Generates a useful figure
%   [handle, axes_h] = jh_figures([width,height],[numCol,numRow],[xmargin,ymargin],FontSize) 
%   opens a new figure with the specified size in centimeters, 
%   TimesNewRoman as 'DefaultAxesFontName', 
%   and FontSize as DefaultAxesFontSize.
%   numCol and numRow give the number of subplots axes, with corresponding
%   handles, e.g. axes_h(1,2) for first column and second row.

%--------------------------------------------------------------------------
if nargout ~= 2
    error('jh_figures: exactly two outputs needed');
end

%--------------------------------------------------------------------------
% Default values
if nargin < 2
    % default: one subplot
    n = [1 1];
elseif numel(n) == 1
    n = [n n];
end

if nargin < 3
    % made with trial and error for nifty plots
    Axpos = [1.75 1.4 0 0]; %[2.40 1.8];
else
    if length(Axpos) == 2
        Axpos = [rowvector(Axpos),0,0];
    end
end

if nargin < 4
    FontSize = 7;
end

%--------------------------------------------------------------------------
% Einheiten einstellen -> macht manchmal hinterher Probleme!
set(0,   'DefaultFigureUnits',    'centimeters'  );
% Figure �ffnen
handle = figure(  'Position',     [5 5 Size(1) Size(2)] );
set(gcf, 'PaperUnits',            'centimeters'  );
set(gcf, 'PaperPosition',         [0 0 Size(1) Size(2)] );
set(gcf, 'PaperSize',             [Size(1) Size(2)] );
%set(gcf, 'DefaultAxesFontName',   'TimesNewRoman');
%set(gcf, 'DefaultAxesFontSize',   FontSize       );
%set(gcf, 'DefaultTextFontName',   'TimesNewRoman');
%set(gcf, 'DefaultTextFontSize',   FontSize       );
set(gcf, 'DefaultLineLineWidth',  1.0            );
set(gcf, 'DefaultAxesLineWidth',  0.5            );
set(gcf, 'DefaultAxesBox',        'On');
set(gcf, 'DefaultAxesUnits',      'centimeters'  );
% set(gcf, 'DefaultAxesUnits',      'normalized');

%--------------------------------------------------------------------------
if( n(1)*n(2) ~= 0)
    % wenn Anzahl an Subplots gegeben

    % Groesse f�r Subplots ausrechnen
    subaxsize  = Size./n - Axpos(1:2) - [0.4 0.114].*Axpos(1:2)/n - Axpos(3:4);

    % Subplot Axes erzeugen und handle an axes_h �bergeben
    for col = 1:n(1)
        for row = 1:n(2)
            subaxpos = Axpos(1:2) + ([1.0 1.0].*Axpos(1:2)+subaxsize).*[col-1 row-1];
            eval(['axes_h(',num2str(col),',',num2str(row),') = jh_axis(subaxpos,subaxsize);']);
        end
    end
else
    axes_h = [];
end % if 

%end % function jh_figures
%--------------------------------------------------------------------------

function handle = jh_axis(axpos,axsize)
% function handle = jh_axes(axpos,axsize, Box);
% Build axes for figures.
% axpos gives position of current axis
% axsize gives size of axis
% Box 'on' or 'off'

Pos = [axpos axsize];

handle = axes('Position',Pos,'Box','on');
%end % function jh_axes
%--------------------------------------------------------------------------