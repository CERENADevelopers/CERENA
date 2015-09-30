% Definition of symbolic variables
syms   S1 S2 S3 X S4
syms   k s1 s2
syms   time
% Define state vector
system.time = time;
system.compartments = {'compartment'};
system.volumes      = [1];
system.state.variable    = [S1; S2; S3; X; S4];
system.state.compartment = {'compartment'; 'compartment'; 'compartment'; 'compartment'; 'compartment'};
system.state.number      = length(system.state.variable);
system.state.type        = {'moment'; 'moment'; 'moment'; 'moment'; 'moment'; };
system.state.name        = {'S1'; 'S2'; 'S3'; 'X'; 'S4'};
system.state.xmin        = [0  0  0  0  0]';
system.state.xmax        = [Inf  Inf  Inf  Inf  Inf]';
system.state.mu0         = [0.03           0           0           0           0]';
system.state.C0          = zeros(system.state.number*(system.state.number+1)/2,1);
system.parameter.variable = [[k; s1; s2]];
system.parameter.name     = {'k'; 's1'; 's2'};

% Define reactions
% (R1)
system.reaction(1).educt      = [S1];
system.reaction(1).product      = [S2];
system.reaction(1).propensity      = S1*k;

% (R2)
system.reaction(2).educt      = [S2];
system.reaction(2).product      = [S4];
system.reaction(2).propensity      = S2*k;

system.output.variable = [S1, S2, S3, X, S4];
system.output.function = [S1, S2, S3, X, S4];
system.output.number   = length(system.output.variable);
system.output.type     = {'moment'; 'moment'; 'moment'; 'moment'; 'moment'; };
system.output.name     = {'S1'; 'S2'; 'S3'; 'X'; 'S4'};

system.input.function = [];
system.input.variable = [];
system.input.number = 0;
system.input.type     = {};
system.input.name     = {};

system = completeSystem(system);
