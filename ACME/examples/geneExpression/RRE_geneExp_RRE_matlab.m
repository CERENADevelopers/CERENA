% function [tout,X,MX,Y] = RRE_geneExp_RRE_matlab(t,theta,kappa) 
function varargout = RRE_geneExp_RRE_matlab(varargin) 

t = varargin{1};
theta = varargin{2};
if(nargin>=3)
    kappa=varargin{3};
   if(length(kappa)==1)
    kappa(2:29)=0;
   end
else
    kappa = zeros(1,29);
end
% Initial conditions
x0 = x0fun(theta,kappa);

% Simulation
[tout,X] = ode15s(@(t,x) rhs(t,x,theta,kappa),t,x0);
Y = rhsO(t,X,theta,kappa);

% Assign output
varargout{1} = tout;
if nargout >= 2
varargout{2} = X;
end
if nargout >= 3
    % Moments of species
    varargout{3} = Y(:,1:4);
end
if nargout >= 4
    % Moments of output variables
    varargout{4} = Y(:,5:end);
end
if nargout >= 5
    error('Too many output arguments.');
end


%% RIGHT-HAND SIDE
function [dxdt] = rhs(t,x,theta,kappa) 

dxdt = [theta(2)*x(2) - theta(1)*x(1) - kappa(1)*theta(7)*x(1)*x(4);...
         theta(1)*x(1) - theta(2)*x(2) + kappa(1)*theta(7)*x(1)*x(4);...
         theta(3)*x(2) - theta(4)*x(3);...
         theta(5)*x(3) - theta(6)*x(4)];

%% OUTPUT MAP
function y = rhsO(t,x,theta,kappa) 

y = [x(:,1),...
         x(:,2),...
         x(:,3),...
         x(:,4),...
         theta(9) + theta(8).*x(:,4)];


%% INITIAL CONDITIONS FOR STATE
function x0 = x0fun(theta,kappa) 

x0 = [(kappa(2)*kappa(16) - kappa(2) + 1)/kappa(1);...
       (kappa(3)*kappa(17))/kappa(1);...
       (kappa(4)*kappa(18) - theta(10)*(kappa(4) - 1))/kappa(1);...
       (kappa(5)*kappa(19))/kappa(1)];


