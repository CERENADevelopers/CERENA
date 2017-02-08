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
%% REACTIONS
%          ->  X1    ,    "k0"
% X1       ->  X2    ,    "k1*X1"
% X2       ->  X1    ,    "k1_1*X2"
% X2       ->  X3    ,    "k2*X2"
% X3       ->  X2    ,    "k2_1*X3"
% X3       ->  X4    ,    "k3*X3"
% X4       ->  X3    ,    "k3_1*X4"
% X4       ->  X5    ,    "k4*X4"
% X5       ->  X4    ,    "k4_1*X5"
% X5       ->  X6    ,    "k5*X5"
% X6       ->  X5    ,    "k5_1*X6"
% X6       ->  X7    ,    "k6*X6"
% X7       ->  X6    ,    "k6_1*X7"
% X7       ->  X8    ,    "k7*X7"
% X8       ->  X7    ,    "k7_1*X8"
% X8       ->  X9    ,    "k8*X8"
% X9       ->  X8    ,    "k8_1*X9"
% X9       ->  X10   ,    "k9*X9"
% X10      ->  X9    ,    "k9_1*X10"
% X10      ->        ,    "k10*X10"
%% Initialization:
% Define states as symbolic variables:
syms X1 X2 X3 X4 X5 X6 X7 X8 X9 X10
% Define parameters as symbolic variables:
syms k0 k1 k1_1 k2 k2_1 k3 k3_1 k4 k4_1 k5 k5_1 k6 k6_1 k7 k7_1 k8 k8_1 k9 k9_1 k10;
% Define output as symbolic variable
% General stuff
syms time

%% Model:
System.time = time;
System.compartments = {'cell'};
System.volumes = 1;
% Define state vector:
System.state.variable = [ X1; X2; X3; X4; X5; X6; X7; X8; X9; X10];
System.state.compartment = {'cell';'cell';'cell';'cell';'cell';'cell';'cell';'cell';'cell';'cell'};
System.state.number   = length(System.state.variable);
System.state.type     = cell(System.state.number,1);
for i=1:System.state.number
    System.state.type(i) = {'moment'};
end
System.state.name     = { 'X1'; 'X2'; 'X3'; 'X4'; 'X5'; 'X6'; 'X7'; 'X8'; 'X9'; 'X10'};
System.state.mu0      = [10;0;0;0;0;0;0;0;0;0];
System.state.C0       = zeros(length(System.state.variable)*(length(System.state.variable)+1)/2,1);

% Define parameter vector:
System.parameter.name = {'k0';'k1';'k1_1'; 'k2';'k2_1'; 'k3';'k3_1';'k4';'k4_1';'k5';'k5_1';'k6';'k6_1';'k7';'k7_1';'k8';'k8_1';'k9';'k9_1';'k10'};
System.parameter.variable = [k0;k1;k1_1; k2;k2_1; k3;k3_1;k4; k4_1; k5; k5_1; k6; k6_1; k7; k7_1; k8; k8_1; k9; k9_1; k10];
% Definition of constant parameters filed:
System.kappa.variable = [];

% Specifying the scale of propensities and parameters
System.scaleIndicator = 'microscopic'; % options are 'microscopic' and 'macroscopic'

% Define propensities:
% (R1)
%          ->  A    ,    "k0"
System.reaction(1).educt      =  [];
System.reaction(1).product    = X1;
System.reaction(1).propensity = k0;
% (R2)
% A        ->  B    ,    "k1*A"
System.reaction(2).educt      =  X1;
System.reaction(2).product    = X2;
System.reaction(2).propensity = k1*X1;
% (R3)
% B       ->  A    ,    "k1_1*B"
System.reaction(3).educt      = X2;
System.reaction(3).product    = X1;
System.reaction(3).propensity = k1_1*X2;
% (R4)
% B       ->  C     ,    "k2*B"
System.reaction(4).educt      = X2;
System.reaction(4).product    = X3;
System.reaction(4).propensity = k2*X2;
% (R5)
% C       ->  B    ,    "k2_1*C"
System.reaction(5).educt      = X3;
System.reaction(5).product    = X2;
System.reaction(5).propensity = k2_1*X3;
% (R6)
% C       ->  D  ,    "k3*C"
System.reaction(6).educt      = X3;
System.reaction(6).product    = X4;
System.reaction(6).propensity = k3*X3;
% (R7)
% D       ->  C    ,    "k3_1*D"
System.reaction(7).educt      = X4;
System.reaction(7).product    = X3;
System.reaction(7).propensity = k3_1*X4;
% (R8)
% D       ->  E  ,    "k4*D"
System.reaction(8).educt      = X4;
System.reaction(8).product    = X5;
System.reaction(8).propensity = k4*X4;
% (R9)
% E       ->  D    ,    "k4_1*E"
System.reaction(9).educt      = X5;
System.reaction(9).product    = X4;
System.reaction(9).propensity = k4_1*X5;
% (R10)
% E       ->  F  ,    "k5*E"
System.reaction(10).educt      = X5;
System.reaction(10).product    = X6;
System.reaction(10).propensity = k5*X5;
% (R11)
% F       ->  E    ,    "k5_1*F"
System.reaction(11).educt      = X6;
System.reaction(11).product    = X5;
System.reaction(11).propensity = k5_1*X6;
% (R12)
% F       ->  G  ,    "k6*F"
System.reaction(12).educt      = X6;
System.reaction(12).product    = X7;
System.reaction(12).propensity = k6*X6;
% (R13)
% G       ->  F    ,    "k6_1*G"
System.reaction(13).educt      = X7;
System.reaction(13).product    = X6;
System.reaction(13).propensity = k6_1*X7;
% (R14)
% G       ->  H  ,    "k7*G"
System.reaction(14).educt      = X7;
System.reaction(14).product    = X8;
System.reaction(14).propensity = k7*X7;
% (R15)
% H       ->  G    ,    "k7_1*H"
System.reaction(15).educt      = X8;
System.reaction(15).product    = X7;
System.reaction(15).propensity = k7_1*X8;
% (R16)
% H       ->  I  ,    "k8*H"
System.reaction(16).educt      = X8;
System.reaction(16).product    = X9;
System.reaction(16).propensity = k8*X8;
% (R17)
% I       ->  H    ,    "k8_1*I"
System.reaction(17).educt      = X9;
System.reaction(17).product    = X8;
System.reaction(17).propensity = k8_1*X9;
% (R18)
% I       ->  J  ,    "k9*I"
System.reaction(18).educt      = X9;
System.reaction(18).product    = X10;
System.reaction(18).propensity = k9*X9;
% (R19)
% J       ->  I    ,    "k9_1*J"
System.reaction(19).educt      = X10;
System.reaction(19).product    = X9;
System.reaction(19).propensity = k9_1*X10;
% (R20)
% J       ->       ,    "k10*J"
System.reaction(20).educt      = X10;
System.reaction(20).product    = [];
System.reaction(20).propensity = k10*X10;

System.output.variable = [X10];
System.output.function = [X10];
System.output.number   = length(System.output.variable);
System.output.name     = {'X10'};

%% Network Structure
System.reactionDistanceInd = [1  2;1  3;1  4;1  5;1  6;1  7;1  8;1  9;1  10;...
                              2  3;2  4;2  5;2  6;2  7;2  8;2  9;2  10;...
                              3  4;3  5;3  6;3  7;3  8;3  9;3  10;...
                              4  5;4  6;4  7;4  8;4  9;4  10;...
                              5  6;5  7;5  8;5  9;5  10;...
                              6  7;6  8;6  9;6  10;...
                              7  8;7  9;7  10;...
                              8  9;8  10;...
                              9  10];
System.reactionDistance = [1;2;3;4;5;6;7;8;9;...
                           1;2;3;4;5;6;7;8;...
                           1;2;3;4;5;6;7;...
                           1;2;3;4;5;6;...
                           1;2;3;4;5;...
                           1;2;3;4;...
                           1;2;3;...
                           1;2;...
                           1];
                       
System.reactionGraph  = [0 1 0 0 0 0 0 0 0 0;1 0 1 0 0 0 0 0 0 0;0 1 0 1 0 0 0 0 0 0;0 0 1 0 1 0 0 0 0 0;...
                         0 0 0 1 0 1 0 0 0 0;0 0 0 0 1 0 1 0 0 0;0 0 0 0 0 1 0 1 0 0;0 0 0 0 0 0 1 0 1 0;...
                         0 0 0 0 0 0 0 1 0 1;0 0 0 0 0 0 0 0 1 0;];      