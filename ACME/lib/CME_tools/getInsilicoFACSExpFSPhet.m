% getInsilicoFACSExpFSP.m
% generate a flow cyometry measurement by simulation the FSP model 'model'
% and sampling the solution. Also the corruption of this measurement with
% noise is possible.
%
% USAGE:
% ======
% [measurement] = getInsilicoFACSExpFSP(model,timevector, mu, sigma)
% [measurement] = getInsilicoFACSExpFSP(model,timevector, mu, sigma, measurands)
% [measurement] = getInsilicoFACSExpFSP(model,timevector,measurands,options)
%
% INPUTS:
% ======
% model ... FSP model used to generate FACS experiments:
% 	.A ... A matrix of FSP (dp/dt = A(p1,theta) p, p(0) = p0). p1 is the
% 	distributrd parameter
% 	.p0 ... initial condition of FSP.
% 	.species ... initial condition of FSP.
%   .index ... matrix describing the mapping between states of the FSP and
%       molecule numbers. The i.th row of 'index' contains the molecule
%       numbers associated to the i.th state of the FSP.
% timevector ... time instances then measurements are taken.
% measurands ... species which are observed.
% options ... options of algorithm
%   .Nt ... number of cells measured at the individual time instances.
%   .ode_solver ... solver options 
%           (default = odeset('RelTol',1e-7,'AbsTol',1e-7)).
%   .noise ... type of measurement noise (default = 'none'),
%          possible choices: 'none', 'mult log'
% mu ... Parameter distribution mean
% sigma ... Parameter distribution standard deviation
% theta ... Parameter vector containing the nondistributed parameters
%
% OUTPUT:
% =======
% measurement ... information about measurement:
%   .time ... vector of time instances then measurements are
%           performed.
%   .data ... sturcture containing the experimental data:
%       .measurands ... names of measurands.
%       .values{i} ... measurement data at time instance i:
%
%                   y     
%                ------    
%               |      |    
%         cell  |      |     
%               |      |      
%                ------
%
% 28/01/2011 - Jan Hasenauer
% modified 01/02/2011 - Jan Hasenauer

% function measurement = getInsilicoFACSExpFSP(model,timevector, mu, sigma, theta, measurands,options)
function measurement = getInsilicoFACSExpFSPhet(varargin)

%% CHECK / ASSIGN INPUTS
% Check required inputs
if nargin >= 5
    if (~isempty(varargin{1}) && ...
        ~isempty(varargin{2}))
        % Assign inputs
        model = varargin{1};
        timevector = rowvector(varargin{2});
        mu=varargin{3};
        sigma=varargin{4};
        theta=varargin{5};
    else
        % Error message:
        error('This routine requires two non-empty inputs!');
    end
else
    % Error message:
    error('This routine requires five inputs!');
end
% Check model
if (~isfield(model,'A')  || ...
    ~isfield(model,'p0') || ...
    ~isfield(model,'species') || ...
    ~isfield(model,'index'))
    % Error message:
    error('Model does not contain all required information.');
else
    if (size(model.index,2) ~= length(model.species))
        % Error message:
        error('Dimension disagreement in model.');
    end
end

% Optional input
measurands = model.species;
if nargin > 5
	measurands = setdefault(varargin{6},measurands);

% Check 'species'
if max(isnan(getindstr(model.species,measurands)))
    % Error message:
    error('Not all species that shall be measured are contained in the model.');
end 
end
% Options
options.Nt = 1000;
options.ode_solver = odeset('RelTol',1e-7,'AbsTol',1e-7);
options.noise.type = 'none';
if nargin == 7
    options = setdefault(varargin{7},options);
end

% Check options.Nt
if length(options.Nt) == 1
    options.Nt = ones(length(timevector),1)*options.Nt;
else
    if length(options.Nt) ~= length(timevector)
        % Error message:
        error('Dimension of options.Nt does not fit to dimension of timevector.');
    end
end

%% SIMULATE FSP
P = getProbInt(model,timevector,theta,mu,sigma,1000);

%% MARGINALIZE DISTRIBUTION
% Determine index of measurands
dim = getindstr(model.species,measurands);
% Perform marginalization
[P,index] = getMarginalization(P,model.index,dim);

%% CONSTRUCT MEASUREMENT
measurement.time = timevector;
measurement.data.measurands = measurands;
measurement.index=index;
%% SAMPLE FSP SOLUTION
% Scaling of FSP for sampling:
P = P*diag(1./sum(P,1));
% Loop: time points
for i = 1:length(timevector)
    % Construct cumulative probability distribution
    cpdf_i = cumsum(P(:,i));
    % Assign measurements to states of the FSP
    index_meas_i = zeros(length(options.Nt(i)),length(measurands));
    for j = 1:options.Nt(i)
        k = max(find(rand >= [0;cpdf_i]));
        index_meas_i(j,:) = index(k,:);
    end
    measurement.data.values{i} = index_meas_i;
    % Add measurement noise
    switch options.noise.type
        case 'none'
                measurement.data.values{i} = index_meas_i;
        case 'mult log'
            for l=1:length(measurement.data.values{i})
                for m=1:length(measurement.data.values{i}(l,:))
            measurement.data.values{i}(l,m) = poissrnd(measurement.data.values{i}(l,m));
                end
            end
        otherwise
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % include measurement noise here! %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end