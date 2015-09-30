% function System = simulateSSA_matlab(modelDefName,t,theta,kappa,Nssa,options)
function System = simulateSSA_matlab(varargin)
if nargin >= 4
    modelDefName = varargin{1};
    t = varargin{2};
    theta = varargin{3};
    if nargin >= 4
        kappa = varargin{4};
        if nargin >= 5
            Nssa = varargin{5};
        else
            Nssa = 10;
        end
    else
        kappa = [];
    end
else
    error('At least three input arguments are required!')
end

options.mode = 'constant';
if nargin >= 6
    options = setdefault(varargin{6},options);
end
eval(modelDefName);
System = completeSystem(System);
System = completeSystemSSA(System);
System.sol = simulate_SSA(System,t,theta,kappa,Nssa,options);
System.sol.theta = theta;
System.sol.kappa = kappa;
