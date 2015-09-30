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
%% Initialization:
% Define states as symbolic variables:
syms g1 g2 g3 g4 LRP r;
% Define parameters as symbolic variables:
syms par1 par2 par3 par4 par5 par6 par7 par8 par9 par10 par11 par12 par13 par14;
syms time
%% Model:
% Define state vector:
System.time = time;
System.compartments = {'cell'};
System.volumes = 1;
System.state.variable = [     g1     ;     g2     ;     g3     ;     g4     ;   LRP  ;    r   ;];
System.state.compartment = { 'cell';'cell'  ;'cell';'cell';'cell';'cell'};
System.state.type     = {'stochastic';'stochastic';'stochastic';'stochastic';'moment';'moment';};
System.state.name     = {    'g1'    ;    'g2'    ;    'g3'    ;    'g4'    ;  'LRP' ;   'r'  ;};
System.state.xmin     = [      0     ;      0     ;      0     ;      0     ;    98  ;    0   ;];
System.state.xmax     = [      1     ;      1     ;      1     ;      1     ;   100  ;   50  ;];
System.state.mu0      = [      1     ;      0     ;      0     ;      0     ;   100  ;    0    ];
System.state.C0       = zeros(length(System.state.variable)*(length(System.state.variable)+1)/2,1);

System.parameter.variable = [  par1  ;   par2   ;   par3   ;  par4  ;   par5   ;   par6   ;  par7  ;   par8   ;   par9   ;...
    par10  ;   par11   ;   par12   ;  par13  ;   par14  ;];
System.parameter.name     = {'\c_1';'\c_2';'\c_3';'\c_4';'\c_5';'\c_6';'\c_7';'\c_8';'\c_9';...
    '\c_{10}';'\c_{11}';'\c_{12}';'\c_{13}';'\c_{14}';};

System.state.constraint = @(X) ((X(1)+X(2)+X(3)+X(4)) == 1);

% Specifying the scale of propensities and parameters
System.scaleIndicator = 'microscopic'; % options are 'microscopic' and 'macroscopic'

% Define propensities:
% (R1)  
System.reaction(1).educt      = [g1,LRP];
System.reaction(1).product    = g2;
System.reaction(1).propensity = par1 * g1 * LRP;
% (R2)  
System.reaction(2).educt      = g2;
System.reaction(2).product    = [g1,LRP];
System.reaction(2).propensity =  g2 * (par2+par3* r/(r+1));
% (R3)  
System.reaction(3).educt      = [g1,LRP];
System.reaction(3).product    = g3;
System.reaction(3).propensity = par4 * g1 * LRP;
% (R4)  
System.reaction(4).educt      = g3;
System.reaction(4).product    = [g1,LRP];
System.reaction(4).propensity = (par5 +par6 * r/(r+1))* g3;
% (R5)  
System.reaction(5).educt      = [g2,LRP];
System.reaction(5).product    = g4;
System.reaction(5).propensity = par7 * g2 * LRP;
% (R6) 
System.reaction(6).educt      = g4;
System.reaction(6).product    = [g2,LRP];
System.reaction(6).propensity = (par8+par9 * r/(r+1)) * g4;
% (R7) 
System.reaction(7).educt      = [g3,LRP];
System.reaction(7).product    = g4;
System.reaction(7).propensity = par10 * g3 * LRP;
% (R8)  
System.reaction(8).educt      = g4;
System.reaction(8).product    = [g3,LRP];
System.reaction(8).propensity = (par11+par12 * r/(r+1)) * g4;
% (R9)  
System.reaction(9).educt      = g2;
System.reaction(9).product    = [g2,r];
System.reaction(9).propensity = par13 * g2;
% (R10) 
System.reaction(10).educt      = r;
System.reaction(10).product    = [];
System.reaction(10).propensity = par14 * r;

% system.output.variable = [ r];
% system.output.function = [ r];
% system.output.name     = {'PapI'};
