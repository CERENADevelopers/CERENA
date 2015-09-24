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
% (R1) 0       -> 40*X1            , c1
% (R2) X1 + X1 -> 0                , c2*X1*(X1-1)/2
% (R3) X1 + X1 -> X1 + X1 + 15*X2  , c3*X1*(X1-1)/2
% (R4)    X2   -> 0                , c4*X2

%% MODEL
% Definition of symbolic variables:
syms X1 X2
syms c1 c2 c3 c4
syms time V

% Define state vector:
system.time = time;
system.compartments = {'V'};
system.volumes = [V];
system.state.variable = [ X1    ;  X2 ];
system.state.compartment = { 'V';'V' };
system.state.type     = {'moment'    ;'moment' };
system.state.name     = {'X1';'X2'};
system.state.xmin     = [      0     ;      0     ];
system.state.xmax     = [   200      ;   200   ];
system.state.mu0      = sym([      10     ;      25     ]);
system.state.C0       = sym([5;0.0001;12.5]);
% Define parameter vector:
system.parameter.variable = [ c1     ; c2     ; c3 ; c4 ;V];
system.parameter.name     = {'c1';'c2';'c3';'c4';'V'};

% Specifying the scale of propensities and parameters
System.scaleIndicator = 'microscopic'; % options are 'microscopic' and 'macroscopic'

% Define propensities:
% (R1)
system.reaction(1).educt      = [];
system.reaction(1).product    = X1*ones(1,40);
system.reaction(1).propensity = c1;
% (R2)
system.reaction(2).educt      = [X1,X1];
system.reaction(2).product    = [];
system.reaction(2).propensity = c2*X1*(X1-1)/2;
% (R3)
system.reaction(3).educt      = [X1,X1];
system.reaction(3).product    = [X1,X1,X2*ones(1,15)];
system.reaction(3).propensity = c3*X1*(X1-1)/2;
% (R4)
system.reaction(4).educt      = [X2];
system.reaction(4).product    = [];
system.reaction(4).propensity = c4*X2;

system.output.variable = [ X2 ];
system.output.function = [ X2 ];
system.output.name     = {'X2'};
