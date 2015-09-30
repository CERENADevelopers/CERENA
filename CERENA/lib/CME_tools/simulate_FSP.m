% function [x,Mx,My,cMx] = simulate_FSP(System,t,theta,System_MCM,options)
function varargout = simulate_FSP(varargin)
if nargin >= 3
    System = varargin{1};
    t = varargin{2};
    theta = varargin{3};
    if nargin >=4
        System_MCM = varargin{4};
    else
        System_MCM = [];
    end
else
    error('At least three input arguments are required!')
end

options.ode15s.reltol = 1e-8;
options.ode15s.abstol = 1e-8;
options.moment_order = 2;
options.moment_order_output = 2;
options.output.calculate = 1;
if nargin >=5
    options = setdefault(varargin{5},options);
end

model.p0 = System.x0;
model.A = System.A(theta);
reltol = options.ode15s.reltol;
abstol = options.ode15s.abstol;
options_ode15s = odeset('RelTol',reltol,'AbsTol',abstol);
x = simulateFSP(model,t,options_ode15s); % Probabilities of FSP states

%% Calculating the overall moments of species and outputs
xo = options.moment_order;
[Mx,Mx_ind] = getMomentsFSP_centered(x,System.index,xo);
yo = options.moment_order_output;
[fun_y,H,My_sym,Iy,n_Iy] = getOutputMoments(System,Mx_ind,[],[],[],[],[],[],options);
My = fun_y(Mx,theta,[]); %fun_y(Mx,theta,kappa)
System.output.order = Iy;
System.state.order = Mx_ind;
%% Calculating the conditional moments of species
if ~isempty(System_MCM)
    [p,cmu,cC] = getMofFSPsol(x',System.index,...
        System_MCM.CMM.state.stochatic.state_index,...
        System_MCM.CMM.state.stochatic.FSP_index,...
        System_MCM.CMM.state.expectation.state_index,...
        System_MCM.CMM.state.expectation.C_index,'central');
    
    % Construct conditional moments
    cMx = [];
    for iy = 1:size(p,2)
        cMx = [cMx,p{iy}];
    end
    for iy = 1:size(p,2)
        cMx = [cMx,cmu{iy},cC{iy}];
    end
end
%% Assembling output
if nargout == 1
    System.sol.t = t;
    System.sol.x = x;
    System.sol.Mx = Mx;
    System.sol.My = My;
    System.sol.theta = theta;
    System.output.order = Iy;
    System.sym.output.function = H;
    if ~isempty(System_MCM)
        System.sol.cMx = cMx;
    end
    varargout{1} = System;
elseif nargout > 1
    varargout{1} = x;
    varargout{2} = Mx;
    varargout{3} = My;
    if nargout > 3
        varargout{4} = cMx;
    end
end
