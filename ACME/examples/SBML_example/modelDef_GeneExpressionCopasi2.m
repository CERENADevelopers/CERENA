% Definition of symbolic variables
syms   DNA_off DNA_on mRNA Protein
syms   tau_off tau_on k_m gamma_m k_p gamma_p tau_on_p
syms   time
% Define state vector
System.time = time;
System.compartments = {'cyt'};
System.volumes      = [1];
System.state.variable    = [DNA_off; DNA_on; mRNA; Protein];
System.state.compartment = {'cyt'; 'cyt'; 'cyt'; 'cyt'};
System.state.number      = length(System.state.variable);
System.state.type        = {'moment'; 'moment'; 'moment'; 'moment'; };
System.state.name        = {'DNA_off'; 'DNA_on'; 'mRNA'; 'Protein'};
System.state.xmin        = transpose([0  0  0  0]);
System.state.xmax        = transpose([10    0   40  100]);
System.state.mu0         = sym(transpose([1   0   4  10]));
System.state.C0          = sym(zeros(System.state.number*(System.state.number+1)/2,1));
System.parameter.variable = [tau_off; tau_on; k_m; gamma_m; k_p; gamma_p; tau_on_p];
System.parameter.name     = {'tau_off'; 'tau_on'; 'k_m'; 'gamma_m'; 'k_p'; 'gamma_p'; 'tau_on_p'};

% Define reactions
% (R1)
System.reaction(1).educt      = [DNA_off];
System.reaction(1).product      = [DNA_on];
System.reaction(1).MacroscopicPropensity      = DNA_off*tau_on_p;

% (R2)
System.reaction(2).educt      = [DNA_on];
System.reaction(2).product      = [DNA_off];
System.reaction(2).MacroscopicPropensity      = DNA_on*tau_off;

% (R3)
System.reaction(3).educt      = [DNA_on];
System.reaction(3).product      = [mRNA, DNA_on];
System.reaction(3).MacroscopicPropensity      = DNA_on*k_m;

% (R4)
System.reaction(4).educt      = [mRNA];
System.reaction(4).product      = [];
System.reaction(4).MacroscopicPropensity      = gamma_m*mRNA;

% (R5)
System.reaction(5).educt      = [mRNA];
System.reaction(5).product      = [mRNA, Protein];
System.reaction(5).MacroscopicPropensity      = k_p*mRNA;

% (R6)
System.reaction(6).educt      = [Protein];
System.reaction(6).product      = [];
System.reaction(6).MacroscopicPropensity      = Protein*gamma_p;

% (R7)
System.reaction(7).educt      = [DNA_off, Protein];
System.reaction(7).product      = [DNA_on, Protein];
System.reaction(7).MacroscopicPropensity      = (DNA_off*Protein*tau_on_p)/1;

System.output.variable = [DNA_off, DNA_on, mRNA, Protein];
System.output.function = [DNA_off, DNA_on, mRNA, Protein];
System.output.number   = length(System.output.variable);
System.output.name     = {'DNA_off'; 'DNA_on'; 'mRNA'; 'Protein'};

System.input.function = [];
System.input.variable = [];
System.input.number = 0;
System.input.name     = {};

