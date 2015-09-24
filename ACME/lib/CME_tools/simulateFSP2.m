% simulateFSP.m
% simulates a ODE system obtained by finite state projection.
%
% USAGE:
% ======
% [P] = simulateFSP(model,t)
% [P] = simulateFSP(model,t,solver_options)
%
% IPUTS:
% ======
% model ... FSP model:
% 	.A ... A matrix of FSP (dp/dt = A p, p(0) = p0).
% 	.p0 ... initial condition of FSP.
% 	.index ... index mapping of FSP.
% t ... time instances then ODE is evaluated.
% solver_options ... solver options 
%           (default = odeset('RelTol',1e-7,'AbsTol',1e-7)).
%
% OUTPUT:
% =======
% P ... n x N matrix, where n is the number of states of the FSP and N is
%       the number of timepoints. Each column provides the state of the FSP
%       at one particular time instance.
% If no output is defined the simulated reponse is plotted.
%
% 28/01/2011 - Jan Hasenauer

% function P = simulateFSP(model,t,solver_options)
function [vargout] = simulateFSP(varargin)

%% CHECK / ASSIGN INPUTS
% Check required inputs
if nargin >= 2
    if (~isempty(varargin{1}) && ...
        ~isempty(varargin{2}))
        % Assign inputs
        model = varargin{1};
        t = rowvector(varargin{2});
    else
        % Error message:
        error('This routine requires two non-empty inputs!');
    end
else
    % Error message:
    error('This routine requires two inputs!');
end
% Check model
if (~isfield(model,'A')  || ...
    ~isfield(model,'p0'))
    % Error message:
    error('Model does not contain all required information.');
else
    if ((size(model.A,1) ~= size(model.A,2)) || ...
        (size(model.A,1) ~= length(model.p0)))
        % Error message:
        error('Dimension disagreement in model.');
    end
end


% Solver options:
solver_options = odeset('RelTol',1e-7,...
                        'AbsTol',1e-7,...
                        'NonNegative',1:length(model.p0));
if nargin == 3
    solver_options = odeset(solver_options,varargin{3});
end
% Set Jacobian
solver_options = odeset(solver_options,'Jacobian',model.A);

%% ASSIGN SYSTEM PROPERTIES
A = model.A;
p0 = model.p0;

%% SIMULATE SYSTEM
if t(1) == 0
    % Simulate system:
    [T,P] = ode45(@(t,p) A*p,t,p0,solver_options);
    P = P';
    % Check for special case that length(t) = 2
    if length(t) == 2
        P = P([1,end],:);
    end    
else
    % Simulate system:
    [T,P] = ode15s(@(t,p) A*p,[0,t],p0,solver_options);
    P = P(2:end,:)';
    % Check for special case that length(t) = 1
    if length(t) == 1
        P = P(end,:);
    end
end

%% WARNING
if sum(P(:,end))/sum(p0) < 0.99
    warning(['The probability within the region of interest ' ...
        'decreased by more than 1 percent!']);
end

%% PLOT SIMULATION IF NO OUTPUT IS DEFINED
switch nargout
    case 0
        if isfield(model,'index')
            plotFSP(P,t,model.index);
        end
    case 1
        vargout = P;
    otherwise
        error('Too many output arguments.');
end

    

end