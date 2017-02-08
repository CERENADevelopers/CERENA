% Definition of symbolic variables
syms   species_1 species_2 species_3 species_4 species_5 species_6 species_7 species_8 species_9 species_10 species_11 species_12 species_13 species_14 species_15 species_16 species_17 species_18 species_19 species_20 species_21 species_22 species_23 species_24 species_25 species_26 species_27 species_28 species_29 species_30 species_31 species_32 species_33 species_34 species_35 species_36 species_37 species_38 species_39
syms   parameter_1 parameter_2 0 0 parameter_5 parameter_6 0 parameter_8 parameter_9 parameter_10 parameter_11 parameter_12 parameter_13 parameter_14 parameter_15 0 parameter_17 parameter_18 parameter_19 parameter_20 parameter_21 parameter_22 parameter_23 parameter_24 parameter_25 parameter_26 0 0 0 parameter_30 parameter_31 parameter_32 parameter_33 parameter_34 parameter_35 parameter_36 parameter_37 parameter_38 0 0 0 0 0 0 0 0 parameter_47 parameter_48
syms   time
% Define state vector
System.time = time;
System.compartments = {'compartment_1'};
System.volumes      = [1];
System.state.variable    = [species_1; species_2; species_3; species_4; species_5; species_6; species_7; species_8; species_9; species_10; species_11; species_12; species_13; species_14; species_15; species_16; species_17; species_18; species_19; species_20; species_21; species_22; species_23; species_24; species_25; species_26; species_27; species_28; species_29; species_30; species_31; species_32; species_33; species_34; species_35; species_36; species_37; species_38; species_39];
System.state.compartment = {'compartment_1'; 'compartment_1'; 'compartment_1'; 'compartment_1'; 'compartment_1'; 'compartment_1'; 'compartment_1'; 'compartment_1'; 'compartment_1'; 'compartment_1'; 'compartment_1'; 'compartment_1'; 'compartment_1'; 'compartment_1'; 'compartment_1'; 'compartment_1'; 'compartment_1'; 'compartment_1'; 'compartment_1'; 'compartment_1'; 'compartment_1'; 'compartment_1'; 'compartment_1'; 'compartment_1'; 'compartment_1'; 'compartment_1'; 'compartment_1'; 'compartment_1'; 'compartment_1'; 'compartment_1'; 'compartment_1'; 'compartment_1'; 'compartment_1'; 'compartment_1'; 'compartment_1'; 'compartment_1'; 'compartment_1'; 'compartment_1'; 'compartment_1'};
System.state.number      = length(System.state.variable);
System.state.type        = {'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; };
System.state.name        = {'Bhi'; 'VattHi'; 'Vex'; 'Blo'; 'VattLo'; 'Ven'; 'Vfus'; 'VpCyt'; 'VpNuc'; 'Rc'; 'P_Rdrp'; 'RcRdrp'; 'P_Np'; 'P_M1'; 'VpNucM1'; 'VpCytM1'; 'Cp'; 'Rv'; 'RvRdrp'; 'Rm1'; 'Rm2'; 'Rm3'; 'Rm4'; 'Rm5'; 'Rm6'; 'Rm7'; 'Rm8'; 'P_Pb1'; 'P_Pb2'; 'P_Pa'; 'P_Nep'; 'P_Ha'; 'P_Na'; 'P_M2'; 'Vrel'; 'total cRNA'; 'total cRNA of a segment'; 'total vRNA'; 'total vRNA of a segment'};
System.state.xmin        = transpose([0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0]);
System.state.xmax        = transpose([1500      0    100  10000      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0]);
System.state.mu0         = transpose([150     0    10  1000     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0]);
System.state.C0          = zeros(System.state.number*(System.state.number+1)/2,1);
System.parameter.variable = [parameter_1; parameter_2; 0; 0; parameter_5; parameter_6; 0; parameter_8; parameter_9; parameter_10; parameter_11; parameter_12; parameter_13; parameter_14; parameter_15; 0; parameter_17; parameter_18; parameter_19; parameter_20; parameter_21; parameter_22; parameter_23; parameter_24; parameter_25; parameter_26; 0; 0; 0; parameter_30; parameter_31; parameter_32; parameter_33; parameter_34; parameter_35; parameter_36; parameter_37; parameter_38; 0; 0; 0; 0; 0; 0; 0; 0; parameter_47; parameter_48];
System.parameter.name     = {'kAttHi'; 'kAttLo'; 'kEn'; 'kFus'; 'kImp'; 'kSynP'; 'kSynV'; 'kSynC'; 'kBindRdrp'; 'kBindNp'; 'kBindM1'; 'kExp'; 'kRdrp'; 'kRel'; 'kDegR'; 'kDegRnp'; 'kDegM'; 'kDegRrdrp'; 'Ffus'; 'Drib'; 'Fspl7'; 'Fspl8'; 'L1'; 'L2'; 'L3'; 'L4'; 'L5'; 'L6'; 'L7'; 'L8'; 'kSynM'; 'kEqHi'; 'kEqLo'};

% Define reactions
System.scaleIndicator = 'macroscopic';
% (R1)
System.reaction(1).educt      = [species_3, species_1];
System.reaction(1).product      = [species_2];
System.reaction(1).propensity      = -1*((parameter_1*species_2)/parameter_47 - parameter_1*species_1*species_3);

% (R2)
System.reaction(2).educt      = [species_3, species_4];
System.reaction(2).product      = [species_5];
System.reaction(2).propensity      = -1*((parameter_2*species_5)/parameter_48 - parameter_2*species_3*species_4);

% (R3)
System.reaction(3).educt      = [species_2];
System.reaction(3).product      = [species_6, species_1];
System.reaction(3).propensity      = 1*parameter_5*species_2;

% (R4)
System.reaction(4).educt      = [species_5];
System.reaction(4).product      = [species_6, species_4];
System.reaction(4).propensity      = 1*parameter_5*species_5;

% (R5)
System.reaction(5).educt      = [species_6];
System.reaction(5).product      = [species_7, species_8];
System.reaction(5).propensity      = 1*parameter_6*species_6;

% (R6)
System.reaction(6).educt      = [species_6];
System.reaction(6).product      = [];
System.reaction(6).propensity      = -(1*parameter_6*species_6*(parameter_23 - 1))/parameter_23;

% (R7)
System.reaction(7).educt      = [species_8];
System.reaction(7).product      = [species_9];
System.reaction(7).propensity      = 1*parameter_8*species_8;

% (R8)
System.reaction(8).educt      = [species_9];
System.reaction(8).product      = [species_9, species_10];
System.reaction(8).propensity      = 1*parameter_11*species_9;

% (R9)
System.reaction(9).educt      = [species_10, species_11];
System.reaction(9).product      = [species_12];
System.reaction(9).propensity      = 1*parameter_12*species_10*species_11;

% (R10)
System.reaction(10).educt      = [species_12, species_13];
System.reaction(10).product      = [species_17];
System.reaction(10).propensity      = 1*parameter_13*species_12*species_13;

% (R11)
System.reaction(11).educt      = [species_9, species_14];
System.reaction(11).product      = [species_15];
System.reaction(11).propensity      = 1*parameter_14*species_9*species_14;

% (R12)
System.reaction(12).educt      = [species_15, species_31];
System.reaction(12).product      = [species_16];
System.reaction(12).propensity      = 1*parameter_15*species_15*species_31;

% (R13)
System.reaction(13).educt      = [species_17];
System.reaction(13).product      = [species_17, species_18];
System.reaction(13).propensity      = 1*parameter_10*species_17;

% (R14)
System.reaction(14).educt      = [species_18, species_11];
System.reaction(14).product      = [species_19];
System.reaction(14).propensity      = 1*parameter_12*species_11*species_18;

% (R15)
System.reaction(15).educt      = [species_19, species_13];
System.reaction(15).product      = [species_9];
System.reaction(15).propensity      = 1*parameter_13*species_13*species_19;

% (R16)
System.reaction(16).educt      = [species_9];
System.reaction(16).product      = [species_9, species_20];
System.reaction(16).propensity      = (1*parameter_38*species_9)/(8*parameter_30);

% (R17)
System.reaction(17).educt      = [species_9];
System.reaction(17).product      = [species_9, species_21];
System.reaction(17).propensity      = (1*parameter_38*species_9)/(8*parameter_31);

% (R18)
System.reaction(18).educt      = [species_9];
System.reaction(18).product      = [species_9, species_22];
System.reaction(18).propensity      = (1*parameter_38*species_9)/(8*parameter_32);

% (R19)
System.reaction(19).educt      = [species_9];
System.reaction(19).product      = [species_9, species_23];
System.reaction(19).propensity      = (1*parameter_38*species_9)/(8*parameter_33);

% (R20)
System.reaction(20).educt      = [species_9];
System.reaction(20).product      = [species_9, species_24];
System.reaction(20).propensity      = (1*parameter_38*species_9)/(8*parameter_34);

% (R21)
System.reaction(21).educt      = [species_9];
System.reaction(21).product      = [species_9, species_25];
System.reaction(21).propensity      = (1*parameter_38*species_9)/(8*parameter_35);

% (R22)
System.reaction(22).educt      = [species_9];
System.reaction(22).product      = [species_9, species_26];
System.reaction(22).propensity      = (1*parameter_38*species_9)/(8*parameter_36);

% (R23)
System.reaction(23).educt      = [species_9];
System.reaction(23).product      = [species_9, species_27];
System.reaction(23).propensity      = (1*parameter_38*species_9)/(8*parameter_37);

% (R24)
System.reaction(24).educt      = [species_21];
System.reaction(24).product      = [species_21, species_28];
System.reaction(24).propensity      = (1*parameter_9*species_21)/parameter_24;

% (R25)
System.reaction(25).educt      = [species_20];
System.reaction(25).product      = [species_20, species_29];
System.reaction(25).propensity      = (1*parameter_9*species_20)/parameter_24;

% (R26)
System.reaction(26).educt      = [species_22];
System.reaction(26).product      = [species_22, species_30];
System.reaction(26).propensity      = (1*parameter_9*species_22)/parameter_24;

% (R27)
System.reaction(27).educt      = [species_28, species_29, species_30];
System.reaction(27).product      = [species_11];
System.reaction(27).propensity      = 1*parameter_17*species_28*species_29*species_30;

% (R28)
System.reaction(28).educt      = [species_24];
System.reaction(28).product      = [species_24, species_13];
System.reaction(28).propensity      = (1*parameter_9*species_24)/parameter_24;

% (R29)
System.reaction(29).educt      = [species_26];
System.reaction(29).product      = [species_26, species_14];
System.reaction(29).propensity      = -(1*parameter_9*species_26*(parameter_25 - 1))/parameter_24;

% (R30)
System.reaction(30).educt      = [species_27];
System.reaction(30).product      = [species_27, species_31];
System.reaction(30).propensity      = (1*parameter_9*parameter_26*species_27)/parameter_24;

% (R31)
System.reaction(31).educt      = [species_23];
System.reaction(31).product      = [species_23, species_32];
System.reaction(31).propensity      = (1*parameter_9*species_23)/parameter_24;

% (R32)
System.reaction(32).educt      = [species_25];
System.reaction(32).product      = [species_25, species_33];
System.reaction(32).propensity      = (1*parameter_9*species_25)/parameter_24;

% (R33)
System.reaction(33).educt      = [species_26];
System.reaction(33).product      = [species_26, species_34];
System.reaction(33).propensity      = (1*parameter_9*parameter_25*species_26)/parameter_24;

% (R34)
System.reaction(34).educt      = [species_16, species_11, species_13, species_14, species_31, species_32, species_33, species_34];
System.reaction(34).product      = [species_35];
System.reaction(34).propensity      = (1*parameter_18*species_11*species_13*species_14*species_16*species_31*species_32*species_33*species_34)/((KmB + species_11)*(KmC + species_32)*(KmD + species_13)*(KmE + species_33)*(KmF + species_14)*(KmG + species_34)*(KmH + species_31));

% (R35)
System.reaction(35).educt      = [species_9];
System.reaction(35).product      = [];
System.reaction(35).propensity      = 1*parameter_20*species_9;

% (R36)
System.reaction(36).educt      = [species_10];
System.reaction(36).product      = [];
System.reaction(36).propensity      = 1*parameter_19*species_10;

% (R37)
System.reaction(37).educt      = [species_18];
System.reaction(37).product      = [];
System.reaction(37).propensity      = 1*parameter_19*species_18;

% (R38)
System.reaction(38).educt      = [species_12];
System.reaction(38).product      = [];
System.reaction(38).propensity      = 1*parameter_22*species_12;

% (R39)
System.reaction(39).educt      = [species_19];
System.reaction(39).product      = [];
System.reaction(39).propensity      = 1*parameter_22*species_19;

% (R40)
System.reaction(40).educt      = [species_17];
System.reaction(40).product      = [];
System.reaction(40).propensity      = 1*parameter_20*species_17;

% (R41)
System.reaction(41).educt      = [species_15];
System.reaction(41).product      = [];
System.reaction(41).propensity      = 1*parameter_20*species_15;

% (R42)
System.reaction(42).educt      = [species_16];
System.reaction(42).product      = [];
System.reaction(42).propensity      = 1*parameter_20*species_16;

% (R43)
System.reaction(43).educt      = [species_20];
System.reaction(43).product      = [];
System.reaction(43).propensity      = 1*parameter_21*species_20;

% (R44)
System.reaction(44).educt      = [species_21];
System.reaction(44).product      = [];
System.reaction(44).propensity      = 1*parameter_21*species_21;

% (R45)
System.reaction(45).educt      = [species_22];
System.reaction(45).product      = [];
System.reaction(45).propensity      = 1*parameter_21*species_22;

% (R46)
System.reaction(46).educt      = [species_23];
System.reaction(46).product      = [];
System.reaction(46).propensity      = 1*parameter_21*species_23;

% (R47)
System.reaction(47).educt      = [species_24];
System.reaction(47).product      = [];
System.reaction(47).propensity      = 1*parameter_21*species_24;

% (R48)
System.reaction(48).educt      = [species_25];
System.reaction(48).product      = [];
System.reaction(48).propensity      = 1*parameter_21*species_25;

% (R49)
System.reaction(49).educt      = [species_26];
System.reaction(49).product      = [];
System.reaction(49).propensity      = 1*parameter_21*species_26;

% (R50)
System.reaction(50).educt      = [species_27];
System.reaction(50).product      = [];
System.reaction(50).propensity      = 1*parameter_21*species_27;

System.output.variable = [species_1, species_2, species_3, species_4, species_5, species_6, species_7, species_8, species_9, species_10, species_11, species_12, species_13, species_14, species_15, species_16, species_17, species_18, species_19, species_20, species_21, species_22, species_23, species_24, species_25, species_26, species_27, species_28, species_29, species_30, species_31, species_32, species_33, species_34, species_35, species_36, species_37, species_38, species_39];
System.output.function = [species_1, species_2, species_3, species_4, species_5, species_6, species_7, species_8, species_9, species_10, species_11, species_12, species_13, species_14, species_15, species_16, species_17, species_18, species_19, species_20, species_21, species_22, species_23, species_24, species_25, species_26, species_27, species_28, species_29, species_30, species_31, species_32, species_33, species_34, species_35, species_36, species_37, species_38, species_39];
System.output.number   = length(System.output.variable);
System.output.name     = {'Bhi'; 'VattHi'; 'Vex'; 'Blo'; 'VattLo'; 'Ven'; 'Vfus'; 'VpCyt'; 'VpNuc'; 'Rc'; 'P_Rdrp'; 'RcRdrp'; 'P_Np'; 'P_M1'; 'VpNucM1'; 'VpCytM1'; 'Cp'; 'Rv'; 'RvRdrp'; 'Rm1'; 'Rm2'; 'Rm3'; 'Rm4'; 'Rm5'; 'Rm6'; 'Rm7'; 'Rm8'; 'P_Pb1'; 'P_Pb2'; 'P_Pa'; 'P_Nep'; 'P_Ha'; 'P_Na'; 'P_M2'; 'Vrel'; 'total cRNA'; 'total cRNA of a segment'; 'total vRNA'; 'total vRNA of a segment'};

System.input.function = [];
System.input.variable = [];
System.input.number = 0;
System.input.name     = {};

