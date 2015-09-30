% Model definition of the reaction network.
% This routine should return a struct, called System, with the following
% fields. Some fields are optional. Consult the ACME documentation
% (ACME/doc/ACME_doc.pdf) for detailed instructions.
%
%       System.time 
%       System.compartments
%       System.volumes
%       System.state.variable
%       System.state.compartment
%       System.state.type
%       System.state.name
%       System.state.xmin
%       System.state.xmax
%       System.state.mu0 
%       System.state.C0  
%       System.state.constraint 
%       System.parameter.variable 
%       System.parameter.name 
%       System.scaleIndicator
%       System.reaction.educt  
%       System.reaction.product
%       System.reaction.propensity  
%       System.output.variable
%       System.output.function
%       System.output.name     = {'Protein'};
%       System.input.function = [];
%       System.input.variable = [];
%       System.input.name     = {};
%% SYSTEM / SETUP
% (R1) DNA_off -> DNA_on  , tau_on*DNA_off
% (R2) DNA_on  -> DNA_off , tau_off*DNA_on
% (R3)       0 -> mRNA    , k_m*DNA_on
% (R4)    mRNA -> 0       , gamma_m*mRNA
% (R5)       0 -> Protein , k_p*mRNA
% (R6) Protein -> 0       , gamma_p*Protein
% (R7) Protein + DNA_off -> Protein + DNA_on , tau_on_p*Protein*DNA_off

%% MODEL
% Definition of symbolic variables:
syms DNA_off DNA_on mRNA Protein %  species 
syms scaledProtein % observables
syms tau_on tau_off k_m gamma_m k_p gamma_p tau_on_p % kinetic rates
syms Omega % volume
syms scaleP offsetP % scaling parameters for the observable 
syms r0 % initial count of mRNA
syms time % time

% Defininition of general fileds:
System.time = time;
System.compartments = {'cell'};
System.volumes = [Omega];

% Defininition of state fileds:
System.state.variable = [ DNA_off    ;  DNA_on    ; mRNA       ; Protein ];
System.state.name = { 'DNA_{off}'    ;  'DNA_{on}'    ; 'mRNA'       ; 'Protein' };
System.state.compartment = { 'cell';'cell'  ;'cell';'cell'};
System.state.type     = {'stochastic';'stochastic';'moment'    ;'moment' };
System.state.xmin     = [      0     ;      0     ;    0       ;    0    ];
System.state.xmax     = [      1     ;      1     ;    40      ;   150   ];
System.state.mu0      = [      1     ;      0     ;    r0      ;    0   ];
System.state.C0       = zeros(length(System.state.variable)*(length(System.state.variable)+1)/2,1);
System.state.constraint = @(x) ((x(1)+x(2)) == 1);

% Definition of parameters field:
System.parameter.variable = [ tau_on     ; tau_off     ; k_m ; gamma_m   ; k_p ; gamma_p   ; tau_on_p; scaleP; offsetP; r0];
% Definition of constant parameters filed:
System.kappa.variable = [Omega];

% Specifying the scale of propensities and parameters
System.scaleIndicator = 'microscopic'; % options are 'microscopic' and 'macroscopic'

% Define propensities:
% (R1)
System.reaction(1).educt      = DNA_off;
System.reaction(1).product    = DNA_on;
System.reaction(1).propensity = tau_on*DNA_off;
% (R2)
System.reaction(2).educt      = DNA_on;
System.reaction(2).product    = DNA_off;
System.reaction(2).propensity = tau_off*DNA_on;
% (R3)
System.reaction(3).educt      = DNA_on;
System.reaction(3).product    = [mRNA,DNA_on];
System.reaction(3).propensity = k_m*DNA_on;
% (R4)
System.reaction(4).educt      = mRNA;
System.reaction(4).product    = [];
System.reaction(4).propensity = gamma_m*mRNA;
% (R5)
System.reaction(5).educt      = mRNA;
System.reaction(5).product    = [Protein,mRNA];
System.reaction(5).propensity = k_p*mRNA;
% (R6)
System.reaction(6).educt      = Protein;
System.reaction(6).product    = [];
System.reaction(6).propensity = gamma_p*Protein;
% (R7)
System.reaction(7).educt      = [DNA_off,Protein];
System.reaction(7).product    = [DNA_on ,Protein];
System.reaction(7).propensity = tau_on_p*Protein*DNA_off;

System.output.variable = [ scaledProtein ];
System.output.name     = {'scaledProtein'};
System.output.function = [ offsetP + scaleP * Protein ];

