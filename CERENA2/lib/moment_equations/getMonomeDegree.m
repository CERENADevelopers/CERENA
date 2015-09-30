% determine_monomial_degree_symbolic.m
% determines the degree of a monomial, and the 
% variables the monomial is build off.
% If a vector of variables is provided the degree
% of the monomial with respect to these varaibles
% is determined
%
% USAGE:
% ======
% [total_degree] = determine_monomial_degree(monomial)
% [degree] = determine_monomial_degree(monomial,xi)
%
% INPUTS:
% =======
% monomial ... monomial as symbolic expression.
%
% OUTPUTS:
% ========
% total_degree ... total degree of monomial.
% degree ... degree of monomial with respect to the
% 	suppplied variables.
%
% 11/09/2010 Jan Hasenauer

% function [degree] = determine_monomial_degree_symbolic(monomial,xi)
function [varargout] = determine_monomial_degree_symbolic(varargin)
%% CHECK INPUTS AND ASSIGN DEFAULTS
if nargin >= 1
    % Check that input is non-empty:
    if ~isempty(varargin{1})
        % Assign inputs
        monomial = varargin{1};
    else
        % Error message:
        error('The first input has to be non-empty.');
    end
else
    % Error message:
    error('Not the correct number of inputs.');
end

%% DETERMINE MONOMES WITH RESPECT TO WITH DEGREE IS CALCULATED
if nargin == 2
    xi = varargin{2};
else
    xi = symvar(monomial);
end

%% PREPROCESSING OF STRING
% Initialize degree
degree = zeros(length(xi),1);
% Loop: entries of xi
for i = 1:length(xi)
    while diff(monomial,xi(i)) ~= 0
        % Update monomial
        monomial = monomial/xi(i);
        % Update degree
        degree(i) = degree(i) + 1;
    end
end

%% ASSIGN OUTPUTS
% Select output type
switch nargin
    % One output:
    case 1
        % Assign degree:
        varargout{1} = sum(degree);
    % Two outputs:
    case 2
        % Assign degree:
        varargout{1} = degree;
end

%% End
end