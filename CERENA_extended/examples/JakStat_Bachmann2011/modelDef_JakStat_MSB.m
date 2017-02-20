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
% Define states as symbolic variables:
syms EpoRJAK2 EpoRpJAK2 p1EpoRpJAK2 p2EpoRpJAK2 p12EpoRpJAK2 ...
     EpoRJAK2_CIS SHP1 SHP1Act STAT5 pSTAT5 npSTAT5 CISnRNA1 ...
     CISnRNA2 CISnRNA3 CISnRNA4 CISnRNA5 CISRNA CIS SOCS3nRNA1 ...
     SOCS3nRNA2 SOCS3nRNA3 SOCS3nRNA4 SOCS3nRNA5 SOCS3RNA SOCS3;
% Define parameters as symbolic variables:
syms JAK2ActEpo SOCS3Inh JAK2EpoRDeaSHP1 EpoRActJAK2 EpoRCISInh EpoRCISRemove SHP1ActEpoR ...
     SHP1Dea STAT5ActJAK2 STAT5ActEpoR CISInh STAT5Imp STAT5Exp CISRNAEqc CISRNATurn ActD ...
     CISRNADelay CISEqc CISTurn CISoe CISEqcOE SOCS3RNAEqc SOCS3RNATurn SOCS3RNADelay ...
     SOCS3Eqc SOCS3Turn SOCS3oe SOCS3EqcOE Omega_cyt Omega_nuc;
% Define output as symbolic variable
syms y_1 y_2 y_3
% Define initial conditions
syms init_EpoRJAK2 init_SHP1 init_STAT5
% Input
syms Epo epo_level
% General stuff
syms time

% Define state vector:
System.time = time;
System.compartments = {'cyt','nuc','frac'};
System.volumes = [Omega_cyt, Omega_nuc, 1];
% Define state vector:
System.state.variable = [ EpoRJAK2; EpoRpJAK2; p1EpoRpJAK2; p2EpoRpJAK2; p12EpoRpJAK2; ...
     EpoRJAK2_CIS; SHP1; SHP1Act; STAT5; pSTAT5; npSTAT5; CISnRNA1; ...
     CISnRNA2; CISnRNA3; CISnRNA4; CISnRNA5; CISRNA; CIS; SOCS3nRNA1; ...
     SOCS3nRNA2; SOCS3nRNA3; SOCS3nRNA4; SOCS3nRNA5; SOCS3RNA; SOCS3;];
System.state.compartment = {'cyt';'cyt';'cyt';'cyt';'cyt';'frac';'cyt';'cyt';...
     'cyt';'cyt';'nuc';'nuc';'nuc';'nuc';'nuc';'nuc';'cyt';'cyt';'nuc';'nuc';...
     'nuc';'nuc';'nuc';'cyt';'cyt'};
System.state.name     = {'EpoRJAK2'; 'EpoRpJAK2'; 'p1EpoRpJAK2';'p2EpoRpJAK2';'p12EpoRpJAK2';...
     'EpoRJAK2\_CIS';'SHP1';'SHP1Act';'STAT5';'pSTAT5';'npSTAT5';'CISnRNA1';...
     'CISnRNA2';'CISnRNA3';'CISnRNA4';'CISnRNA5';'CISRNA';'CIS';'SOCS3nRNA1';...
     'SOCS3nRNA2';'SOCS3nRNA3';'SOCS3nRNA4';'SOCS3nRNA5';'SOCS3RNA';'SOCS3';};
System.state.number   = length(System.state.variable);
System.state.mu0      = sym(zeros(System.state.number,1));
System.state.mu0(1)   = init_EpoRJAK2 * Omega_cyt;
System.state.mu0(7)   = init_SHP1 * Omega_cyt;
System.state.mu0(9)   = init_STAT5 * Omega_cyt;
% Define parameter vector:
dynParameter.variable = [JAK2ActEpo; SOCS3Inh; JAK2EpoRDeaSHP1; EpoRActJAK2; EpoRCISInh; EpoRCISRemove; SHP1ActEpoR; ...
     SHP1Dea; STAT5ActJAK2; STAT5ActEpoR; CISInh; STAT5Imp; STAT5Exp; CISRNAEqc; CISRNATurn; ActD; ...
     CISRNADelay; CISEqc; CISTurn; CISoe; CISEqcOE; SOCS3RNAEqc; SOCS3RNATurn; SOCS3RNADelay; ...
     SOCS3Eqc; SOCS3Turn; SOCS3oe; SOCS3EqcOE];
dynParameter.name     = {'JAK2ActEpo'; 'SOCS3Inh'; 'JAK2EpoRDeaSHP1'; 'EpoRActJAK2'; 'EpoRCISInh'; 'EpoRCISRemove'; 'SHP1ActEpoR'; ...
     'SHP1Dea'; 'STAT5ActJAK2'; 'STAT5ActEpoR'; 'CISInh'; 'STAT5Imp'; 'STAT5Exp'; 'CISRNAEqc'; 'CISRNATurn'; 'ActD'; ...
     'CISRNADelay'; 'CISEqc'; 'CISTurn'; 'CISoe'; 'CISEqcOE'; 'SOCS3RNAEqc'; 'SOCS3RNATurn'; 'SOCS3RNADelay'; ...
     'SOCS3Eqc'; 'SOCS3Turn'; 'SOCS3oe'; 'SOCS3EqcOE'};
[dynParameter.name,sort_ind] = sort(dynParameter.name);
System.parameter.name = [dynParameter.name;'init_EpoRJAK2'; 'init_SHP1'; 'init_STAT5'; 'Omega_cyt'; 'Omega_nuc'; 'epo_level'];
System.parameter.variable = [dynParameter.variable(sort_ind);init_EpoRJAK2; init_SHP1; init_STAT5; Omega_cyt; Omega_nuc; epo_level];
System.parameter.kinetic  = ones(length(System.parameter.variable),1);
System.parameter.kinetic(32:end-1) = 0;
System.parameter.kinetic(end) = 1;

System.scaleIndicator = 'macroscopic';
% Define propensities:
% (R1)
% EpoRJAK2        ->  EpoRpJAK2    ,    "JAK2ActEpo * EpoRJAK2 * Epo / (1 + SOCS3Inh * SOCS3)"
System.reaction(1).educt      =  EpoRJAK2;
System.reaction(1).product    = EpoRpJAK2;
System.reaction(1).propensity = JAK2ActEpo * EpoRJAK2 * epo_level / (1 + SOCS3Inh * SOCS3);%JAK2ActEpo * EpoRJAK2 * Epo / (1 + SOCS3Inh * SOCS3);
% (R2)
% EpoRpJAK2       ->  EpoRJAK2     ,    "JAK2EpoRDeaSHP1 * EpoRpJAK2 * SHP1Act"
System.reaction(2).educt      = EpoRpJAK2;
System.reaction(2).product    = EpoRJAK2;
System.reaction(2).propensity = JAK2EpoRDeaSHP1 * EpoRpJAK2 * SHP1Act;
% (R3)
% EpoRpJAK2       ->  p1EpoRpJAK2  ,    "EpoRActJAK2 * EpoRpJAK2 / (1 + SOCS3Inh * SOCS3)"
System.reaction(3).educt      = EpoRpJAK2;
System.reaction(3).product    = p1EpoRpJAK2;
System.reaction(3).propensity = EpoRActJAK2 * EpoRpJAK2 / (1 + SOCS3Inh * SOCS3);
% (R4)
% EpoRpJAK2       ->  p2EpoRpJAK2  ,    "3*EpoRActJAK2 * EpoRpJAK2 / (1 + SOCS3Inh * SOCS3) / (EpoRCISInh * EpoRJAK2_CIS + 1)"
System.reaction(4).educt      = EpoRpJAK2;
System.reaction(4).product    = p2EpoRpJAK2;
% system.reaction(4).propensity = 3*EpoRActJAK2 * EpoRpJAK2 / (1 + SOCS3Inh * SOCS3/Av/Omega_cyt) / (EpoRCISInh * EpoRJAK2_CIS/Av/Omega_cyt + 1);
System.reaction(4).propensity = 3*EpoRActJAK2 * EpoRpJAK2 / (1 + SOCS3Inh * SOCS3) / (EpoRCISInh * EpoRJAK2_CIS + 1);
% (R5)
% p1EpoRpJAK2     ->  p12EpoRpJAK2 ,    "3*EpoRActJAK2 * p1EpoRpJAK2 / (1 + SOCS3Inh * SOCS3) / (EpoRCISInh * EpoRJAK2_CIS + 1)"
System.reaction(5).educt      = p1EpoRpJAK2;
System.reaction(5).product    = p12EpoRpJAK2;
% system.reaction(5).propensity = 3*EpoRActJAK2 * p1EpoRpJAK2 / (1 + SOCS3Inh * SOCS3/Av/Omega_cyt) / (EpoRCISInh * EpoRJAK2_CIS/Av/Omega_cyt + 1);
System.reaction(5).propensity = 3*EpoRActJAK2 * p1EpoRpJAK2 / (1 + SOCS3Inh * SOCS3) / (EpoRCISInh * EpoRJAK2_CIS + 1);
% (R6)
% p2EpoRpJAK2     ->  p12EpoRpJAK2 ,    "EpoRActJAK2 * p2EpoRpJAK2 / (1 + SOCS3Inh * SOCS3)"
System.reaction(6).educt      = p2EpoRpJAK2;
System.reaction(6).product    = p12EpoRpJAK2;
System.reaction(6).propensity = EpoRActJAK2 * p2EpoRpJAK2 / (1 + SOCS3Inh * SOCS3);
% (R7)
% p1EpoRpJAK2     ->  EpoRJAK2     ,    "JAK2EpoRDeaSHP1 * p1EpoRpJAK2 * SHP1Act"
System.reaction(7).educt      = p1EpoRpJAK2;
System.reaction(7).product    = EpoRJAK2;
System.reaction(7).propensity = JAK2EpoRDeaSHP1 * p1EpoRpJAK2 * SHP1Act;
% (R8)
% p2EpoRpJAK2     ->  EpoRJAK2     ,    "JAK2EpoRDeaSHP1 * p2EpoRpJAK2 * SHP1Act"
System.reaction(8).educt      = p2EpoRpJAK2;
System.reaction(8).product    = EpoRJAK2;
System.reaction(8).propensity = JAK2EpoRDeaSHP1 * p2EpoRpJAK2 * SHP1Act;
% (R9)
% p12EpoRpJAK2    ->  EpoRJAK2     ,    "JAK2EpoRDeaSHP1 * p12EpoRpJAK2 * SHP1Act"
System.reaction(9).educt      = p12EpoRpJAK2;
System.reaction(9).product    = EpoRJAK2;
System.reaction(9).propensity = JAK2EpoRDeaSHP1 * p12EpoRpJAK2 * SHP1Act;
% (R10)
% EpoRJAK2_CIS    ->               ,    "EpoRCISRemove * EpoRJAK2_CIS * (p1EpoRpJAK2 + p12EpoRpJAK2)"
System.reaction(10).educt      = EpoRJAK2_CIS;
System.reaction(10).product    = [];
System.reaction(10).propensity = EpoRCISRemove * EpoRJAK2_CIS * (p1EpoRpJAK2 + p12EpoRpJAK2);

% system.reaction(10).educt      = EpoRJAK2_CIS;
% system.reaction(10).product    = [];
% system.reaction(10).propensity = EpoRCISRemove * EpoRJAK2_CIS * p1EpoRpJAK2 /Omega_cyt;
% 
% system.reaction(11).educt      = EpoRJAK2_CIS;
% system.reaction(11).product    = [];
% system.reaction(11).propensity = EpoRCISRemove * EpoRJAK2_CIS * p12EpoRpJAK2/Omega_cyt;
% system.reaction(10).rate       = EpoRCISRemove * (p1EpoRpJAK2 + p12EpoRpJAK2)/Omega_cyt;
% system.reaction(10).propensity = EpoRCISRemove * EpoRJAK2_CIS * (p1EpoRpJAK2 + p12EpoRpJAK2);
% system.reaction(10).parameter  = [EpoRCISRemove,Av,Omega_cyt];
% (R11)
% SHP1            ->  SHP1Act      ,    "SHP1ActEpoR * SHP1 * (EpoRpJAK2 + p1EpoRpJAK2 + p2EpoRpJAK2 + p12EpoRpJAK2)"
System.reaction(11).educt      = SHP1;
System.reaction(11).product    = SHP1Act;
System.reaction(11).propensity = SHP1ActEpoR * SHP1 * (EpoRpJAK2 + p1EpoRpJAK2 + p2EpoRpJAK2 + p12EpoRpJAK2);

% system.reaction(12).educt      = SHP1;
% system.reaction(12).product    = SHP1Act;
% system.reaction(12).propensity = SHP1ActEpoR * SHP1 * EpoRpJAK2/Omega_cyt;
% 
% system.reaction(13).educt      = SHP1;
% system.reaction(13).product    = SHP1Act;
% system.reaction(13).propensity = SHP1ActEpoR * SHP1 * p1EpoRpJAK2/Omega_cyt;
% 
% system.reaction(14).educt      = SHP1;
% system.reaction(14).product    = SHP1Act;
% system.reaction(14).propensity = SHP1ActEpoR * SHP1 * p2EpoRpJAK2/Omega_cyt;
% 
% system.reaction(15).educt      = SHP1;
% system.reaction(15).product    = SHP1Act;
% system.reaction(15).propensity = SHP1ActEpoR * SHP1 * p12EpoRpJAK2/Omega_cyt;
% system.reaction(11).rate       = SHP1ActEpoR * (EpoRpJAK2 + p1EpoRpJAK2 + p2EpoRpJAK2 + p12EpoRpJAK2)/Omega_cyt;
% system.reaction(11).parameter  = [SHP1ActEpoR,Av,Omega_cyt];
% (R12)
% SHP1Act         ->  SHP1         ,    "SHP1Dea * SHP1Act"
System.reaction(12).educt      = SHP1Act;
System.reaction(12).product    = SHP1;
System.reaction(12).propensity = SHP1Dea * SHP1Act;
% (R13)
% STAT5           ->  pSTAT5       ,    "STAT5ActJAK2 * STAT5 * (EpoRpJAK2 + p1EpoRpJAK2 + p2EpoRpJAK2 + p12EpoRpJAK2) / (1 + SOCS3Inh * SOCS3)"
System.reaction(13).educt      = STAT5;
System.reaction(13).product    = pSTAT5;
System.reaction(13).propensity = STAT5ActJAK2 * STAT5 * (EpoRpJAK2 + p1EpoRpJAK2 + p2EpoRpJAK2 + p12EpoRpJAK2)/ (1 + SOCS3Inh * SOCS3);

% system.reaction(17).educt      = STAT5;
% system.reaction(17).product    = pSTAT5;
% system.reaction(17).propensity = STAT5ActJAK2 * STAT5 * EpoRpJAK2 /Omega_cyt;%/ (1 + SOCS3Inh * SOCS3/Omega_cyt);
% 
% 
% system.reaction(18).educt      = STAT5;
% system.reaction(18).product    = pSTAT5;
% system.reaction(18).propensity = STAT5ActJAK2 * STAT5 * p1EpoRpJAK2 /Omega_cyt;%/ (1 + SOCS3Inh * SOCS3/Omega_cyt);
% 
% 
% system.reaction(19).educt      = STAT5;
% system.reaction(19).product    = pSTAT5;
% system.reaction(19).propensity = STAT5ActJAK2 * STAT5 * p2EpoRpJAK2 /Omega_cyt;%/ (1 + SOCS3Inh * SOCS3/Omega_cyt);
% 
% 
% system.reaction(20).educt      = STAT5;
% system.reaction(20).product    = pSTAT5;
% system.reaction(20).propensity = STAT5ActJAK2 * STAT5 * p12EpoRpJAK2 /Omega_cyt;%/ (1 + SOCS3Inh * SOCS3/Omega_cyt);

% system.reaction(13).rate       = STAT5ActJAK2 * (EpoRpJAK2 + p1EpoRpJAK2 + p2EpoRpJAK2 + p12EpoRpJAK2) /Omega_cyt/ (1 + SOCS3Inh * SOCS3/Omega_cyt);
% system.reaction(13).parameter  = [STAT5ActJAK2,SOCS3Inh,Av,Omega_cyt];
% (R14)
% STAT5           ->  pSTAT5       ,    "STAT5ActEpoR * STAT5 * (p1EpoRpJAK2 + p12EpoRpJAK2)^2 / (1 + SOCS3Inh * SOCS3) / (1 + CISInh * CIS)"
System.reaction(14).educt      = STAT5;
System.reaction(14).product    = pSTAT5;
System.reaction(14).propensity = STAT5ActEpoR * STAT5 * (p1EpoRpJAK2 + p12EpoRpJAK2)^2/ (1 + SOCS3Inh * SOCS3) / (1 + CISInh * CIS);
% (R15)
% pSTAT5          ->  npSTAT5      ,    "STAT5Imp * pSTAT5"
System.reaction(15).educt      = pSTAT5;
System.reaction(15).product    = npSTAT5;
System.reaction(15).propensity = STAT5Imp * pSTAT5;
% (R16)
% npSTAT5         ->  STAT5        ,    "STAT5Exp * npSTAT5"
System.reaction(16).educt      = npSTAT5;
System.reaction(16).product    = STAT5;
System.reaction(16).propensity = STAT5Exp * npSTAT5;
% (R17)
%                 ->  CISnRNA1     ,    "CISRNAEqc * CISRNATurn * npSTAT5 * (1-ActD)"
System.reaction(17).educt      = [];
System.reaction(17).product    = CISnRNA1;
System.reaction(17).propensity = CISRNAEqc * CISRNATurn * npSTAT5 * (1-ActD);
% (R18)
% CISnRNA1        ->  CISnRNA2     ,    "CISRNADelay * CISnRNA1"
System.reaction(18).educt      = CISnRNA1;
System.reaction(18).product    = CISnRNA2;
System.reaction(18).propensity = CISRNADelay * CISnRNA1;
% (R19)
% CISnRNA2        ->  CISnRNA3     ,    "CISRNADelay * CISnRNA2"
System.reaction(19).educt      = CISnRNA2;
System.reaction(19).product    = CISnRNA3;
System.reaction(19).propensity = CISRNADelay * CISnRNA2;
% (R20)
% CISnRNA3        ->  CISnRNA4     ,    "CISRNADelay * CISnRNA3"
System.reaction(20).educt      = CISnRNA3;
System.reaction(20).product    = CISnRNA4;
System.reaction(20).propensity = CISRNADelay * CISnRNA3;
% (R21)
% CISnRNA4        ->  CISnRNA5     ,    "CISRNADelay * CISnRNA4"
System.reaction(21).educt      = CISnRNA4;
System.reaction(21).product    = CISnRNA5;
System.reaction(21).propensity = CISRNADelay * CISnRNA4;
% (R22)
% CISnRNA5        ->  CISRNA       ,    "CISRNADelay * CISnRNA5"
System.reaction(22).educt      = CISnRNA5;
System.reaction(22).product    = CISRNA;
System.reaction(22).propensity = CISRNADelay * CISnRNA5;
% (R23)
% CISRNA          ->               ,    "CISRNATurn * CISRNA"
System.reaction(23).educt      = CISRNA;
System.reaction(23).product    = [];
System.reaction(23).propensity = CISRNATurn * CISRNA;
% (R24)
%                 ->  CIS          ,    "CISEqc * CISTurn * CISRNA"
System.reaction(24).educt      = [];
System.reaction(24).product    = CIS;
System.reaction(24).propensity = CISEqc * CISTurn * CISRNA;
% (R25)
% CIS             ->               ,    "CISTurn * CIS"
System.reaction(25).educt      = CIS;
System.reaction(25).product    = [];
System.reaction(25).propensity = CISTurn * CIS;
% (R26)
%                 ->  CIS          ,    "CISoe * CISEqcOE * CISTurn"
System.reaction(26).educt      = [];
System.reaction(26).product    = CIS;
System.reaction(26).propensity = CISoe * CISEqcOE * CISTurn;
% (R27)
%                 ->  SOCS3nRNA1   ,    "SOCS3RNAEqc * SOCS3RNATurn * npSTAT5 * (1-ActD)"
System.reaction(27).educt      = [];
System.reaction(27).product    = SOCS3nRNA1;
System.reaction(27).propensity = SOCS3RNAEqc * SOCS3RNATurn * npSTAT5 * (1-ActD);
% (R28)
% SOCS3nRNA1      ->  SOCS3nRNA2   ,    "SOCS3RNADelay * SOCS3nRNA1"
System.reaction(28).educt      = SOCS3nRNA1;
System.reaction(28).product    = SOCS3nRNA2;
System.reaction(28).propensity = SOCS3RNADelay * SOCS3nRNA1;
% (R29)
% SOCS3nRNA2      ->  SOCS3nRNA3   ,    "SOCS3RNADelay * SOCS3nRNA2"
System.reaction(29).educt      = SOCS3nRNA2;
System.reaction(29).product    = SOCS3nRNA3;
System.reaction(29).propensity = SOCS3RNADelay * SOCS3nRNA2;
% (R30)
% SOCS3nRNA3      ->  SOCS3nRNA4   ,    "SOCS3RNADelay * SOCS3nRNA3"
System.reaction(30).educt      = SOCS3nRNA3;
System.reaction(30).product    = SOCS3nRNA4;
System.reaction(30).propensity = SOCS3RNADelay * SOCS3nRNA3;
% (R31)
% SOCS3nRNA4      ->  SOCS3nRNA5   ,    "SOCS3RNADelay * SOCS3nRNA4"
System.reaction(31).educt      = SOCS3nRNA4;
System.reaction(31).product    = SOCS3nRNA5;
System.reaction(31).propensity = SOCS3RNADelay * SOCS3nRNA4;
% (R32)
% SOCS3nRNA5      ->  SOCS3RNA     ,    "SOCS3RNADelay * SOCS3nRNA5"
System.reaction(32).educt      = SOCS3nRNA5;
System.reaction(32).product    = SOCS3RNA;
System.reaction(32).propensity = SOCS3RNADelay * SOCS3nRNA5;
% (R33)
% SOCS3RNA        ->               ,    "SOCS3RNATurn * SOCS3RNA"
System.reaction(33).educt      = SOCS3RNA;
System.reaction(33).product    = [];
System.reaction(33).propensity = SOCS3RNATurn * SOCS3RNA;
% (R34)
%                 ->  SOCS3        ,    "SOCS3Eqc * SOCS3Turn * SOCS3RNA"
System.reaction(34).educt      = [];
System.reaction(34).product    = SOCS3;
System.reaction(34).propensity = SOCS3Eqc * SOCS3Turn * SOCS3RNA;
% (R35)
% SOCS3           ->               ,    "SOCS3Turn * SOCS3"
System.reaction(35).educt      = SOCS3;
System.reaction(35).product    = [];
System.reaction(35).propensity = SOCS3Turn * SOCS3;
% (R36)
%                 ->  SOCS3        ,    "SOCS3oe * SOCS3EqcOE * SOCS3Turn"
System.reaction(36).educt      = [];
System.reaction(36).product    = SOCS3;
System.reaction(36).propensity = SOCS3oe * SOCS3EqcOE * SOCS3Turn;


System.output.variable = [ y_1 ; y_2 ; y_3 ];
System.output.function = [2 * (EpoRpJAK2 + p1EpoRpJAK2 + p2EpoRpJAK2 + p12EpoRpJAK2);...
                          16 * (p1EpoRpJAK2 + p2EpoRpJAK2 + p12EpoRpJAK2);...
                          (STAT5+pSTAT5)];
System.output.name     = {'pJAK2';'pEpoR';'tSTAT5'};
% 
% System.input.variable = [ Epo ];
% System.input.function = [ epo_level ];
% System.input.name     = {'Epo'};
