% function [X,Y,mX,mY,CX,CY] = simulate_SSA(system,t,theta,kappa,Nssa,options)
function [varargout] = simulate_SSA(varargin)

if nargin >=3
    system = varargin{1};
    t = varargin{2};
    theta = varargin{3};
    if nargin >=4
        kappa = varargin{4};
        if length(kappa) == system.kappa.nk1
            kappa(system.kappa.nk1+1:length(system.kappa.variable)) = 0;
        end
    else
        kappa = zeros(length(system.kappa.variable),1);
    end
    if nargin >=5 && ~isempty(varargin{5})
        Nssa = varargin{5};
    else
        Nssa = 10;
    end
else
    error('At least three input arguments are required!')
end

options.mode = 'constant';
options.scale = 'absolute';
if nargin >=6
    options = setdefault(varargin{6},options);
end
v = system.v;
S = system.S;
H = system.H;

X = zeros(length(t),system.state.number,Nssa);
if isfield(system,'output')
    Y = zeros(length(t),system.output.number,Nssa);
end
x0 = system.state.mu0;
if isfield(system,'input')
    x0 = mysubs(x0,system.input.variable,system.input.function);
end
x0 = mysubs(x0,system.state.fmu0_sym,system.state.fmu0);
x0 = mysubs(x0,system.parameter.variable,theta);
x0 = mysubs(x0,system.kappa.variable,kappa);
x0 = floor(double(x0));

dnssa = floor(Nssa/10);

for i = 1:Nssa
    if (mod(i,dnssa)==1) && (i~=1)
        disp([num2str(i-1),' SSA runs completed.'])
    end
    switch options.mode
        case 'constant'
            X(:,:,i) = simulateSSA(@(x) v(t,x,theta),S,t,x0)'; 
            if isfield(system,'output')
                Y(:,:,i) = H(t,squeeze(X(:,:,i)),theta);
            else
                Y = [];
            end
        case 'time-dependent'
            X(:,:,i) = simulateSSA_t(@(t,x) v(t,x,theta),S,t,x0)';
            if isfield(system,'output')
                Y(:,:,i) = H(t,squeeze(X(:,:,i)),theta);
            else
                Y = [];
            end
    end
    if (mod(i,100)==0)
        save([num2str(i),'SSA'])
    end
end
if strcmp(options.scale,'concentration')
    volumes = mysubs(system.state.volume,system.parameter.variable,theta);
    volumes = double(mysubs(volumes,system.kappa.variable,kappa));
    for i=1:system.state.number
        X(:,i,:) = X(:,i,:)./volumes(i);
    end
end

if nargout == 1 || nargout >= 3
    mX = mean(X,3);
end
if nargout == 1 || nargout >= 4
    mY = mean(Y,3);
end
if nargout == 1 || nargout >= 5
    CX = var(X,[],3);
end
if nargout == 1 || nargout >= 6
    CY = var(Y,[],3);
end
if nargout == 1
    sol.x = X;
    sol.y = Y;
    sol.mean_x = mX;
    sol.mean_y = mY;
    sol.var_x = CX;
    sol.var_y = CY;
    sol.t = t;
    varargout{1} = sol;
elseif nargout > 1   
    varargout{1} = X;
    varargout{2} = Y;
    if nargout >=3 
        varargout{3} = mX;
        if nargout >=4
            varargout{4} = mY;
            if nargout >=5
                varargout{5} = CX;
                if nargout >=6
                    varargout{6} = CY;
                end
            end
        end
    end
end

end

% better subs
function out = mysubs(in, old, new)
    if(~isnumeric(in) && ~isempty(old) && ~isempty(findsym(in)))
        matVer = ver('MATLAB');
        if(str2double(matVer.Version)>=8.1)
            out = subs(in, old(:), new(:));
        else
            out = subs(in, old(:), new(:), 0);
        end
    else
        out = in;
    end
end