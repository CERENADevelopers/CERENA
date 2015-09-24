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
%% MODEL
% Definition of symbolic variables:
syms STAT pSTAT pSTAT_pSTAT npSTAT_npSTAT nSTAT1 nSTAT2 nSTAT3 nSTAT4 nSTAT5;
syms k_1 k_2 k_3 k_4 STAT_total_conc u Omega_cyt Omega_nuc 
syms totpSTAT totSTAT offset_pSTAT scale_pSTAT offset_tSTAT scale_tSTAT
syms y_1 y_2
syms pEpoR 
syms pEpoR sp1 sp2 sp3 sp4 sp5
syms a b
syms time

% Define state vector:
System.time = time;
System.compartments = {'cyt','nuc'};
System.volumes = [Omega_cyt,Omega_nuc];
System.state.variable = [STAT; pSTAT; pSTAT_pSTAT; npSTAT_npSTAT; nSTAT1; nSTAT2; nSTAT3; nSTAT4; nSTAT5];
System.state.compartment = { 'cyt';'cyt'  ;'cyt';'nuc';'nuc';'nuc';'nuc';'nuc';'nuc'};
System.state.number   = length(System.state.variable);
System.state.type     = {'moment';'moment';'moment';'moment';'moment';'moment';'moment';'moment'; 'moment'};
System.state.name     = {'STAT';'pSTAT';'pSTAT:pSTAT';'npSTAT:npSTAT';'nSTAT_1';'nSTAT_2';'nSTAT_3';'nSTAT_4';'nSTAT_5'};
System.state.xmin     = [0;0;0;0;0;0;0;0;0];
System.state.mu0      = [ STAT_total_conc * Omega_cyt;0;0;0;0;0;0;0;0];
System.state.C0       = sym(zeros(System.state.number*(System.state.number+1)/2,1));
% Define parameter vector:
System.parameter.variable = [ k_1 ; k_2 ; k_3 ; k_4 ; STAT_total_conc; Omega_cyt; Omega_nuc; ...
                              offset_pSTAT; offset_tSTAT; scale_pSTAT; scale_tSTAT; sp1; sp2; sp3; sp4; sp5;a;b];
System.parameter.name     = {'k_1';'k_2';'k_3';'k_4';'STAT_{total, conc}'; 'Omega_{cyt}'; 'Omega_{nuc}';...
                             'offset_pSTAT'; 'offset_tSTAT'; 'scale_pSTAT'; 'scale_tSTAT';'sp1'; 'sp2'; 'sp3'; 'sp4'; 'sp5';'a';'b'};
System.scaleIndicator = 'macroscopic';
% Define propensities:
% (R1)
System.reaction(1).educt      =  STAT;
System.reaction(1).product    = pSTAT;
System.reaction(1).propensity = k_1*STAT*pEpoR;
% (R2)
System.reaction(2).educt      = [pSTAT,pSTAT];
System.reaction(2).product    = pSTAT_pSTAT;
System.reaction(2).propensity = k_2/STAT_total_conc*pSTAT^2;
% (R3)
System.reaction(3).educt      = pSTAT_pSTAT;
System.reaction(3).product    = npSTAT_npSTAT;
System.reaction(3).propensity = k_3*pSTAT_pSTAT;
% (R4)
System.reaction(4).educt      = npSTAT_npSTAT;
System.reaction(4).product    = [nSTAT1,nSTAT1];
System.reaction(4).propensity = k_4*npSTAT_npSTAT;
% (R5)
System.reaction(5).educt      = nSTAT1;
System.reaction(5).product    = nSTAT2;
System.reaction(5).propensity = k_4*nSTAT1;
% (R6)
System.reaction(6).educt      = nSTAT2;
System.reaction(6).product    = nSTAT3;
System.reaction(6).propensity = k_4*nSTAT2;
% (R7)
System.reaction(7).educt      = nSTAT3;
System.reaction(7).product    = nSTAT4;
System.reaction(7).propensity = k_4*nSTAT3;
% (R8)
System.reaction(8).educt      = nSTAT4;
System.reaction(8).product    = nSTAT5;
System.reaction(8).propensity = k_4*nSTAT4;
% (R9)
System.reaction(9).educt      = nSTAT5;
System.reaction(9).product    = STAT;
System.reaction(9).propensity = k_4*nSTAT5;

System.output.variable = [ totpSTAT ; totSTAT; ];
System.output.function = [offset_pSTAT + scale_pSTAT/STAT_total_conc *(pSTAT+2*pSTAT_pSTAT);offset_tSTAT + scale_tSTAT/STAT_total_conc*(STAT+pSTAT+2*pSTAT_pSTAT)];
System.output.name     = {'totpSTAT';'totSTAT'};

System.input.variable = [  pEpoR ];
System.input.function = [ (b*time^2)*exp(-a*time) ];
% system.input.function = [spline_pos5(time, 0.0, sp1, 5.0, sp2, 10.0, sp3, 20.0, sp4, 60.0, sp5, 0, 0.0)];
System.input.name     = {'pEpoR'};
