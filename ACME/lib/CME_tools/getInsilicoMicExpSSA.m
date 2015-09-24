% getInsilicoMicExpSSA.m
% generate a fluorescense microscopy measurement by simulation the
%stochastic model consisting of stoichiometric matrix S and propensity
%matrix v and sampling the solution. 
%Also the corruption of this measurement with noise is possible.
%
% USAGE:
% ======
% [measurement] = getInsilicoMicExpSSA(model,timevector)
% [measurement] = getInsilicoMicExpSSA(model,timevector,measurands)
% [measurement] = getInsilicoMicExpSSA(model,timevector,measurands,options)
%
% INPUTS:
% ======
% model ... model of the system used to generate Experiments
%   .propensities ... Propensity matrix.
%   .S ... stoichiometric matrix.
%   .x0 ... initial condition of the System.
%   .species ... Cell array containing the symbolic expressions of the species
% timevector ... time instances the measurements are taken.
% measurands ... symbolic expressions of the species which are observed.
% options ... options of algorithm
%   .Nt ... Number of cells measured
%   .noise ... type of measurement noise (default = 'none'),
%          possible choices: 'none', 'mult log'
%
% OUTPUT:
% =======
% measurement ... information about measurement:
% .data ... structure containing the experimental data:
%   .time ... vector of time instances the measurements are
%             performed.
%   .values{i} ... measurement data of experiment i:
%       
% .species ... names of measured species.
%       
%
% 28/10/2011 - Sebastian Waider


% function measurement = getInsilicoMicExpSSA(model,timevector,measurands,options)
function measurement = getInsilicoMicExpSSA(varargin)

%% CHECK / ASSIGN INPUTS
% Check required inputs
if nargin >= 2
    if (~isempty(varargin{1}) && ...
        ~isempty(varargin{2}))
        % Assign inputs
        model=varargin{1};
        v = model.propensities;
        S=model.S;
        timevector = rowvector(varargin{2});
        x0= model.x0;
        species=model.species;
    else
        % Error message:
        error('This routine requires four non-empty inputs!');
    end
else
    % Error message:
    error('This routine requires four inputs!');
end
if(nargin>=3)
    measurands=varargin{3};
else
    measurands=species;
end
% Check 'species'
 if max(isnan(getindstr(model.species,measurands)))
     % Error message:
error('Not all species that shall be measured are contained in the model.');
end 

% Options
options.noise.type = 'none';
options.Nt=100;
if nargin == 4
    options = setdefault(varargin{4},options);
end


%% SIMULATE SSA
for i=1:options.Nt
X = simulateSSAcontinuous(v,S,timevector,x0);

%% Only keep the species given in measurands
if(length(species)~=length(measurands))
    Xm{i}=X.data.values(getindstr(species,measurands),:)';
else Xm{i}=X.data.values';
end
end
%% CONSTRUCT MEASUREMENT
measurement.data.time = timevector;
measurement.data.values = Xm;
measurement.species=measurands;

%% CORRUPT SOLUTION WITH NOISE
% Loop: time points
for i = 1:length(options.Nt)
    for j=1:length(timevector)
    % Add measurement noise
    switch options.noise.type
        case 'none'
              
        case 'mult log'
            for q=1:length(measurement.data.values{i}(j,:))
            measurement.data.values{i}(j,q) = ...
                poissrnd(measurement.data.values{i}(j,q));
            end
        otherwise
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % include measurement noise here! %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    end
end