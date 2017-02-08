% function System = simulateFSP_matlab(modelDefName,t,theta,kappa,System_MCM,options)
function System = simulateFSP_matlab(varargin)
if nargin >= 3
    modelDefName = varargin{1};
    t = varargin{2};
    theta = varargin{3};
    if nargin >= 4
        kappa = varargin{4};
    else
        kappa = [];
    end
else
    error('at least three input arguments are required!')
end
eval(modelDefName)
System.state.mu0 = subs(System.state.mu0,System.parameter.variable,theta);
if isfield(System,'kappa')
    System.state.mu0 = subs(System.state.mu0,System.kappa.variable,kappa);
end
System = completeSystemFSP(System);
theta_FSP = [theta;kappa];
if nargin >= 5
    System_MCM = varargin{5};
else
    System_MCM = [];
end
if nargin >= 6
    options = varargin{6};
    System = simulate_FSP(System,t,theta_FSP,System_MCM,options);
else
    System = simulate_FSP(System,t,theta_FSP,System_MCM);
end
    
