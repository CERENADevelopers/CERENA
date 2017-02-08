% getInsilicoWesternblotFSP.m
% generate a Westernblot measurement by simulation the FSP model 'model'
% and sampling the solution. Also the corruption of this measurement with
% noise is possible.
%
% USAGE:
% ======
% [measurement] = getInsilicoWesternblotFSP(model,timevector)
% [measurement] = getInsilicoWesternblotFSP(model,timevector,measurands)
% [measurement] = getInsilicoWesternblotFSP(model,timevector,measurands,options)
%
% INPUTS:
% ======
% model ... FSP model used to generate FACS experiments:
% 	.A ... A matrix of FSP (dp/dt = A p, p(0) = p0).
% 	.p0 ... initial condition of FSP.
% 	.species ... initial condition of FSP.
%   .index ... matrix describing the mapping between states of the FSP and
%       molecule numbers. The i.th row of 'index' contains the molecule
%       numbers associated to the i.th state of the FSP.
% timevector ... time instances then measurements are taken.
% measurands ... species which are observed.
% options ... options of algorithm
%   .ode_solver ... solver options 
%           (default = odeset('RelTol',1e-7,'AbsTol',1e-7)).
%   .noise ... type of measurement noise (default = 'none'),
%          possible choices: 'none', 'mult log'
%   .nt ... Number of experiments
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
% 07/11/2011 - Jan Hasenauer


% function measurement = getInsilicoWesternblot(model,timevector,measurands,options)
function measurement = getInsilicoWesternblotFSP(varargin)

%% CHECK / ASSIGN INPUTS
% Check required inputs
if nargin >= 2
    if (~isempty(varargin{1}) && ...
        ~isempty(varargin{2}))
        % Assign inputs
        model = varargin{1};
        timevector = rowvector(varargin{2});
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
    ~isfield(model,'p0') || ...
    ~isfield(model,'species') || ...
    ~isfield(model,'index'))
    % Error message:
    error('Model does not contain all required information.');
else
    if ((size(model.A,1) ~= size(model.A,2)) || ...
        (size(model.A,1) ~= length(model.p0)) || ...
        (size(model.index,2) ~= length(model.species)))
        % Error message:
        error('Dimension disagreement in model.');
    end
end

% Optional input
measurands = model.species;
% if nargin >= 3
% 	measurands = setdefault(varargin{3},measurands);
%     
% 
% % Check 'species'
% if max(isnan(getindstr(model.species,measurands)))
%     % Error message:
%     error('Not all species that shall be measured are contained in the model.');
% end 
% end


% Options
options.ode_solver = odeset('RelTol',1e-7,'AbsTol',1e-7);
options.noise.type = 'none';
options.nt=1000;
if nargin == 4
    options = setdefault(varargin{4},options);
end


%% SIMULATE FSP
P = simulateFSP(model,timevector,options.ode_solver);

%% MARGINALIZE DISTRIBUTION
% Determine index of measurands
dim = getindstr(model.species,measurands);
% Perform marginalization
[P,index] = getMarginalization(P,model.index,dim);

%% CONSTRUCT MEASUREMENT
measurement.data.time = timevector;
measurement.measurands = measurands;

%% GENERATE WESTERNBLOT
for i=1:options.nt
measurement.data.values{i}=getExpectationFSP(P,index);
    % Add measurement noise
    
    switch options.noise.type
        case 'none'
            
        case 'mult log'
            for j=1:size(measurement.data.values{i},2)
            measurement.data.values{i}(:,j) = ...
                measurement.data.values{i}(:,j).*exp(options.noise.sigma.*randn(size(measurement.data.values{i}(:,j)))...
                    +options.noise.mu);  
            end  
        otherwise
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % include measurement noise here! %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    end

end