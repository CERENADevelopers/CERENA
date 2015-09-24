% Definition of symbolic variables
syms   S1 S2 S3 X S4
syms   k k1 s1 s2 s3 c c1
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
system.state.mu0         = [3  6  0  0  0]';
system.state.C0          = zeros(system.state.number*(system.state.number+1)/2,1);
system.parameter.variable = [[k; k1; s1; s2; s3; c; c1]];
system.parameter.name     = {'k'; 'k1'; 's1'; 's2'; 's3'; 'c'; 'c1'};

% Define reactions
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
