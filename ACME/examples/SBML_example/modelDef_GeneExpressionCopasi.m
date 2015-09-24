% Definition of symbolic variables
syms   DNA_off DNA_on mRNA Protein
syms   tau_on k_m gamma_m k_p tau_on_p gamma_p tau_off
syms   time
% Define state vector
system.time = time;
system.compartments = {'cyt'};
system.volumes      = [1];
system.state.variable    = [DNA_off; DNA_on; mRNA; Protein];
system.state.compartment = {'cyt'; 'cyt'; 'cyt'; 'cyt'};
system.state.number      = length(system.state.variable);
system.state.type        = {'moment'; 'moment'; 'moment'; 'moment'; };
system.state.name        = {'DNA_off'; 'DNA_on'; 'mRNA'; 'Protein'};
system.state.xmin        = transpose([0  0  0  0]);
system.state.xmax        = transpose([10  10  10  10]);
system.state.mu0         = transpose([1  1  1  1]);
system.state.C0          = zeros(system.state.number*(system.state.number+1)/2,1);
system.parameter.variable = [tau_on; k_m; gamma_m; k_p; tau_on_p; gamma_p; tau_off];
system.parameter.name     = {'tau_on'; 'k_m'; 'gamma_m'; 'k_p'; 'tau_on_p'; 'gamma_p'; 'tau_off'};

% Define reactions
% (R1)
system.reaction(1).educt      = [DNA_off];
system.reaction(1).product      = [DNA_on];
system.reaction(1).MacroscopicPropensity      = DNA_off*cyt*tau_on;

% (R2)
system.reaction(2).educt      = [DNA_on];
system.reaction(2).product      = [DNA_off];
system.reaction(2).MacroscopicPropensity      = DNA_on*cyt*tau_off;

% (R3)
system.reaction(3).educt      = [];
system.reaction(3).product      = [mRNA];
system.reaction(3).MacroscopicPropensity      = cyt*k_m;

% (R4)
system.reaction(4).educt      = [mRNA];
system.reaction(4).product      = [];
system.reaction(4).MacroscopicPropensity      = cyt*gamma_m*mRNA;

% (R5)
system.reaction(5).educt      = [mRNA];
system.reaction(5).product      = [mRNA, Protein];
system.reaction(5).MacroscopicPropensity      = cyt*k_p*mRNA;

% (R6)
system.reaction(6).educt      = [Protein];
system.reaction(6).product      = [];
system.reaction(6).MacroscopicPropensity      = Protein*cyt*gamma_p;

% (R7)
system.reaction(7).educt      = [Protein, DNA_off];
system.reaction(7).product      = [Protein, DNA_on];
system.reaction(7).MacroscopicPropensity      = DNA_off*Protein*cyt*tau_on_p;

system.output.variable = [DNA_off, DNA_on, mRNA, Protein];
system.output.function = [DNA_off, DNA_on, mRNA, Protein];
system.output.number   = length(system.output.variable);
system.output.name     = {'DNA_off'; 'DNA_on'; 'mRNA'; 'Protein'};

system.input.function = [];
system.input.variable = [];
system.input.number = 0;
system.input.name     = {};

