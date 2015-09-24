% function logL = logLikelihood_FACS_FSP(theta,system,histogram,options)
function varargout = logLikelihood_FACS_FSP(varargin)

%% CHECK/ASSIGN INPUTS:
if nargin >= 3
	theta = varargin{1};
	system = varargin{2};
	histogram = varargin{3};
else
	error('Not enough inputs!');
end

% Assign defaults for options
options.solver_options = [];
options.sign = 'positive';
options.scale = 'lin';
options.grad_ind = [1:length(theta)]';
if nargin == 4
	options = setdefault(varargin{4},options);
end

%% Simulation and marginalization
if nargout == 1
    % Initialization of model
    system.A = system.A(theta);
    options_cvode = CVodeSetOptions('JacobianFn', @Jac,...
                                    'RelTol',1.e-10,...
                                    'AbsTol',1.e-10,...
                                    'MaxNumSteps', 200000,...
                                    'LinearSolver','GMRES',...
                                    'UserData',system);
    CVodeInit(@rhsfn,'BDF','Newton',0,system.p0,options_cvode);

    % Simulation
    P = zeros(length(system.p0),length(histogram.time));
    P(:,1) = system.p0;
    [~,~,P(:,2:end)] = CVode(histogram.time(2:end), 'Normal');

    % Marginalization
    dim = getindstr(system.state.name,histogram.data.measurands);
    P = getMarginalization(P,system.index,dim);
else
    % Initialization of model
    system.A = system.A(theta);
    options_cvode = CVodeSetOptions('JacobianFn', @Jac,...
                                    'RelTol',1.e-10,...
                                    'AbsTol',1.e-10,...
                                    'MaxNumSteps', 200000,...
                                    'LinearSolver','GMRES',...
                                    'UserData',system);
    CVodeInit(@rhsfn,'BDF','Newton',0,system.p0,options_cvode);

    % Initialization of sensitvity model
    FSAoptions = CVodeSensSetOptions('method','Simultaneous',...
                                     'ErrControl',true,...
                                     'ParamScales',1:length(theta));
    CVodeSensInit(length(theta),@rhsSfn,zeros(length(system.p0),length(theta)),FSAoptions);

    % Simulation
    P = zeros(length(system.p0),length(histogram.time));
    PS = zeros(length(system.p0),length(theta),length(histogram.time));
    P(:,1) = system.p0;
    [~,~,P(:,2:end),PS(:,:,2:end)] = CVode(histogram.time(2:end), 'Normal');

    % Marginalization
    dim = getindstr(system.state.name,histogram.data.measurands);
    [P,index] = getMarginalization(P,system.index,dim);
    PSr = zeros(size(index,1),length(theta),length(histogram.time));
    for i = 1:length(theta)
        PSr(:,i,:) = getMarginalization(squeeze(PS(:,i,:)),system.index,dim);
    end
end

%% Evaluation of log-likelihood function
logL = 0;
if nargout == 2
    grad = zeros(length(theta),1);
end
% Loop: time points
for k = 1:length(histogram.time)
    % Reassignment
    h = histogram.data.values(:,k);
    % Constant
    I = log(1:histogram.data.cellsMeasured(k));
    for i = 1:length(h)
        logL = logL - sum(I(1:h(i)));
    end
    % Model dependent
    ind = find(h);    
    logL = logL + h(ind)'*log(P(ind,k));
    if nargout == 2
        for i = 1:length(theta)
            grad(i) = grad(i) + (h(ind)'*(PSr(ind,i,k)./P(ind,k)))';
        end
    end
end

%% ASSIGN OUTPUT
if nargout == 2
    if strcmp(options.scale,'log')
        grad = grad(options.grad_ind).*theta(options.grad_ind);
    else
        grad = grad(options.grad_ind);
    end
end

switch nargout
    % One output
    case {0,1}
        switch  options.sign
            case 'positive'
                varargout{1} =  logL;
            case 'negative'
                varargout{1} = -logL;
        end
    % Two outputs
    case 2
        switch  options.sign
            case 'positive'
                varargout{1} =  logL;
                varargout{2} =  grad;
            case 'negative'
                varargout{1} = -logL;
                varargout{2} = -grad;
        end
end


end


    
function [J,flag,new_data] = Jac(~,~,~,y,system)
    J = system.A*y;
    flag = 0;
    new_data = [];
end

function [yd,flag,new_data] = rhsfn(~,y,system)
    yd = system.A*y;
    flag = 0;
    new_data = [];
end

function [ySd, flag, new_data] = rhsSfn(t,y,yd,yS,system)
    ySd = system.A*yS;
    for i = 1:length(system.A_i)
        ySd(:,i) = ySd(:,i) + system.A_i{i}*y;
    end
    flag = 0;
    new_data = [];
end
