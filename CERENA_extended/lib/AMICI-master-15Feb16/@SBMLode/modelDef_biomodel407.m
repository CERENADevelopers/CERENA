% Definition of symbolic variables
syms   TNFR_E TNF_E TNF_TNFR_E TNFR RIP TRADD TRAF2 FADD TNF_TNFR_TRADD TNFRC1 TNFRCint1 TNFRCint2 TNFRCint3 TNFRC2 FLIP TNFRC2_FLIP TNFRC2_pCasp8 TNFRC2_FLIP_FLIP TNFRC2_pCasp8_pCasp8 TNFRC2_FLIP_pCasp8 TNFRC2_FLIP_pCasp8_RIP_TRAF2 IKK IKKa A20 NFkB IkBa IkBa_NFkB PIkBa NFkB_N IkBa_N IkBa_NFkB_N A20_mRNA IkBa_mRNA XIAP_mRNA FLIP_mRNA BAR XIAP pCasp8 pCasp3 pCasp6 Casp8 Casp3 Casp6 BAR_Casp8 XIAP_Casp3 PARP cPARP
syms ka_1_J1  ka_2_J2  ka_3_J3  ka_4_J4  kd_4_J4  ka_5_J5  kd_5_J5  ka_6_J6  kd_6_J6  ka_7_J7  kd_7_J7  ka_8_J8  ka_9_J9  ka_10_J10  ka_11_J11  ka_12_J12  ka_13_J13  ka_14_J14  ka_15_J15  ka_16_J16  ka_17_J17  ka_18_J18  kd_18_J18  ka_19_J19  ka_20_J20  ka_21_J21  ka_22_J22  ka_23_J23  ka_24_J24  ka_25_J25  ka_26_J26  ka_27_J27  ka_28_J28  ka_29_J29  ka_30_J30  ka_31_J31  ka_32_J32  ka_33_J33  ka_34_J34  ka_35_J35  kd_35_J35  ka_36_J36  kd_36_J36  ka_37_J37  kd_37_J37  ka_38_J38  kd_38_J38  ka_39_J39  kd_39_J39  ka_40_J40  ka_41_J41  ka_42_J42  ka_43_J43  ka_44_J44  ka_45_J45  ka_46_J46  ka_47_J47  ka_48_J48  ka_49_J49  ka_50_J50  ka_51_J51  ka_52_J52  ka_53_J53  ka_54_J54  ka_55_J55  ka_56_J56  ka_57_J57  ka_58_J58  ka_59_J59  kd_59_J59  ka_60_J60  ka_61_J61  ka_62_J62  ka_63_J63  ka_64_J64  ka_65_J65  ka_66_J66  ka_67_J67  ka_68_J68  kd_68_J68  ka_69_J69  kd_69_J69  ka_70_J70  kd_70_J70  ka_71_J71  ka_72_J72  ka_73_J73  ka_74_J74  ka_75_J75  kd_75_J75  ka_76_J76  ka_77_J77  kd_77_J77  ka_78_J78  ka_79_J79  ka_80_J80  ka_81_J81  ka_82_J82  kd_82_J82  ka_83_J83  ka_84_J84  ka_85_J85  ka_86_J86  ka_87_J87  ka_88_J88 kd_88_J88
syms   time
% Define state vector
System.time = time;
System.compartments = {'cytoplasm', 'extracellular', 'nucleus'};
System.volumes      = [3.2           1344          1.056];
System.state.variable    = [TNFR_E; TNF_E; TNF_TNFR_E; TNFR; RIP; TRADD; TRAF2; FADD; TNF_TNFR_TRADD; TNFRC1; TNFRCint1; TNFRCint2; TNFRCint3; TNFRC2; FLIP; TNFRC2_FLIP; TNFRC2_pCasp8; TNFRC2_FLIP_FLIP; TNFRC2_pCasp8_pCasp8; TNFRC2_FLIP_pCasp8; TNFRC2_FLIP_pCasp8_RIP_TRAF2; IKK; IKKa; A20; NFkB; IkBa; IkBa_NFkB; PIkBa; NFkB_N; IkBa_N; IkBa_NFkB_N; A20_mRNA; IkBa_mRNA; XIAP_mRNA; FLIP_mRNA; BAR; XIAP; pCasp8; pCasp3; pCasp6; Casp8; Casp3; Casp6; BAR_Casp8; XIAP_Casp3; PARP; cPARP];
System.state.compartment = {'extracellular'; 'extracellular'; 'extracellular'; 'cytoplasm'; 'cytoplasm'; 'cytoplasm'; 'cytoplasm'; 'cytoplasm'; 'cytoplasm'; 'cytoplasm'; 'cytoplasm'; 'cytoplasm'; 'cytoplasm'; 'cytoplasm'; 'cytoplasm'; 'cytoplasm'; 'cytoplasm'; 'cytoplasm'; 'cytoplasm'; 'cytoplasm'; 'cytoplasm'; 'cytoplasm'; 'cytoplasm'; 'cytoplasm'; 'cytoplasm'; 'cytoplasm'; 'cytoplasm'; 'cytoplasm'; 'nucleus'; 'nucleus'; 'nucleus'; 'nucleus'; 'nucleus'; 'nucleus'; 'nucleus'; 'cytoplasm'; 'cytoplasm'; 'cytoplasm'; 'cytoplasm'; 'cytoplasm'; 'cytoplasm'; 'cytoplasm'; 'cytoplasm'; 'cytoplasm'; 'cytoplasm'; 'cytoplasm'; 'cytoplasm'};
System.state.number      = length(System.state.variable);
System.state.type        = {'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; };
System.state.name        = {'TNFR_E'; 'TNF_E'; 'TNF:TNFR_E'; 'TNFR'; 'RIP'; 'TRADD'; 'TRAF2'; 'FADD'; 'TNF:TNFR:TRADD'; 'TNFRC1'; 'TNFRCint1'; 'TNFRCint2'; 'TNFRCint3'; 'TNFRC2'; 'FLIP'; 'TNFRC2:FLIP'; 'TNFRC2:pCasp8'; 'TNFRC2:FLIP:FLIP'; 'TNFRC2:pCasp8:pCasp8'; 'TNFRC2:FLIP:pCasp8'; 'TNFRC2:FLIP:pCasp8:RIP:TRAF2'; 'IKK'; 'IKKa'; 'A20'; 'NFkB'; 'IkBa'; 'IkBa:NFkB'; 'PIkBa'; 'NFkB_N'; 'IkBa_N'; 'IkBa:NFkB_N'; 'A20_mRNA'; 'IkBa_mRNA'; 'XIAP_mRNA'; 'FLIP_mRNA'; 'BAR'; 'XIAP'; 'pCasp8'; 'pCasp3'; 'pCasp6'; 'Casp8'; 'Casp3'; 'Casp6'; 'BAR:Casp8'; 'XIAP:Casp3'; 'PARP'; 'cPARP'};
System.state.xmin        = transpose([0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0]);
System.state.xmax        = transpose([0.05        2.688            0       0.0028       2.0256       2.9344       3.3056       3.0944            0            0            0            0            0            0     0.320472            0            0            0            0            0            0          6.4            0      1.04434   0.00115365    0.0101518     0.151032            0   0.00691431     0.013839  0.000900189  0.000556657  0.000531517   0.00219646   0.00139056       2.8789      78.3371           32            8         0.64            0            0            0            0            0      16.6667            0]);
System.state.mu0         = transpose([0.005      0.2688           0     0.00028     0.20256     0.29344     0.33056     0.30944           0           0           0           0           0           0    0.032047           0           0           0           0           0           0        0.64           0     0.10443  0.00011537   0.0010152    0.015103           0  0.00069143   0.0013839  9.0019e-05  5.5666e-05  5.3152e-05  0.00021965  0.00013906     0.28789      7.8337         3.2         0.8       0.064           0           0           0           0           0      1.6667           0]);
System.state.C0          = zeros(System.state.number*(System.state.number+1)/2,1);

System.parameter.variable = [ ka_1_J1; ka_2_J2; ka_3_J3; ka_4_J4; kd_4_J4; ka_5_J5; kd_5_J5; ka_6_J6; kd_6_J6; ka_7_J7; kd_7_J7; ka_8_J8; ka_9_J9; ka_10_J10; ka_11_J11; ka_12_J12; ka_13_J13; ka_14_J14; ka_15_J15; ka_16_J16; ka_17_J17; ka_18_J18; kd_18_J18; ka_19_J19; ka_20_J20; ka_21_J21; ka_22_J22; ka_23_J23; ka_24_J24; ka_25_J25; ka_26_J26; ka_27_J27; ka_28_J28; ka_29_J29; ka_30_J30; ka_31_J31; ka_32_J32; ka_33_J33; ka_34_J34; ka_35_J35; kd_35_J35; ka_36_J36; kd_36_J36; ka_37_J37; kd_37_J37; ka_38_J38; kd_38_J38; ka_39_J39; kd_39_J39; ka_40_J40; ka_41_J41; ka_42_J42; ka_43_J43; ka_44_J44; ka_45_J45; ka_46_J46; ka_47_J47; ka_48_J48; ka_49_J49; ka_50_J50; ka_51_J51; ka_52_J52; ka_53_J53; ka_54_J54; ka_55_J55; ka_56_J56; ka_57_J57; ka_58_J58; ka_59_J59; kd_59_J59; ka_60_J60; ka_61_J61; ka_62_J62; ka_63_J63; ka_64_J64; ka_65_J65; ka_66_J66; ka_67_J67; ka_68_J68; kd_68_J68; ka_69_J69; kd_69_J69; ka_70_J70; kd_70_J70; ka_71_J71; ka_72_J72; ka_73_J73; ka_74_J74; ka_75_J75; kd_75_J75; ka_76_J76; ka_77_J77; kd_77_J77; ka_78_J78; ka_79_J79; ka_80_J80; ka_81_J81; ka_82_J82; kd_82_J82; ka_83_J83; ka_84_J84; ka_85_J85; ka_86_J86; ka_87_J87; ka_88_J88; kd_88_J88];

% Define reactions
System.scaleIndicator = 'macroscopic';
% (R1)
System.reaction(1).educt      = [TNFR];
System.reaction(1).product      = [TNFR_E];
System.reaction(1).propensity      = J1_ka_1*TNFR;

% (R2)
System.reaction(2).educt      = [];
System.reaction(2).product      = [TNFR];
System.reaction(2).propensity      = J2_ka_2;

% (R3)
System.reaction(3).educt      = [TNFR_E];
System.reaction(3).product      = [];
System.reaction(3).propensity      = J3_ka_3*TNFR_E;

% (R4)
System.reaction(4).educt      = [];
System.reaction(4).product      = [RIP];
System.reaction(4).propensity      = J4_ka_4 - J4_kd_4*RIP;

% (R5)
System.reaction(5).educt      = [];
System.reaction(5).product      = [TRADD];
System.reaction(5).propensity      = J5_ka_5 - J5_kd_5*TRADD;

% (R6)
System.reaction(6).educt      = [];
System.reaction(6).product      = [TRAF2];
System.reaction(6).propensity      = J6_ka_6 - J6_kd_6*TRAF2;

% (R7)
System.reaction(7).educt      = [];
System.reaction(7).product      = [FADD];
System.reaction(7).propensity      = J7_ka_7 - FADD*J7_kd_7;

% (R8)
System.reaction(8).educt      = [TNF_TNFR_E];
System.reaction(8).product      = [];
System.reaction(8).propensity      = J8_ka_8*TNF_TNFR_E;

% (R9)
System.reaction(9).educt      = [TNF_TNFR_TRADD];
System.reaction(9).product      = [];
System.reaction(9).propensity      = J9_ka_9*TNF_TNFR_TRADD;

% (R10)
System.reaction(10).educt      = [TNFRC1];
System.reaction(10).product      = [];
System.reaction(10).propensity      = J10_ka_10*TNFRC1;

% (R11)
System.reaction(11).educt      = [TNFRC2];
System.reaction(11).product      = [];
System.reaction(11).propensity      = J11_ka_11*TNFRC2;

% (R12)
System.reaction(12).educt      = [TNFRC2_FLIP];
System.reaction(12).product      = [];
System.reaction(12).propensity      = J12_ka_12*TNFRC2_FLIP;

% (R13)
System.reaction(13).educt      = [TNFRC2_FLIP_FLIP];
System.reaction(13).product      = [];
System.reaction(13).propensity      = J13_ka_13*TNFRC2_FLIP_FLIP;

% (R14)
System.reaction(14).educt      = [TNFRC2_pCasp8];
System.reaction(14).product      = [];
System.reaction(14).propensity      = J14_ka_14*TNFRC2_pCasp8;

% (R15)
System.reaction(15).educt      = [TNFRC2_pCasp8_pCasp8];
System.reaction(15).product      = [];
System.reaction(15).propensity      = J15_ka_15*TNFRC2_pCasp8_pCasp8;

% (R16)
System.reaction(16).educt      = [TNFRC2_FLIP_pCasp8];
System.reaction(16).product      = [];
System.reaction(16).propensity      = J16_ka_16*TNFRC2_FLIP_pCasp8;

% (R17)
System.reaction(17).educt      = [TNFRC2_FLIP_pCasp8_RIP_TRAF2];
System.reaction(17).product      = [];
System.reaction(17).propensity      = J17_ka_17*TNFRC2_FLIP_pCasp8_RIP_TRAF2;

% (R18)
System.reaction(18).educt      = [TNFR_E, TNF_E];
System.reaction(18).product      = [TNF_TNFR_E];
System.reaction(18).propensity      = J18_ka_18*TNF_E*TNFR_E - J18_kd_18*TNF_TNFR_E;

% (R19)
System.reaction(19).educt      = [TNF_TNFR_E, TRADD];
System.reaction(19).product      = [TNF_TNFR_TRADD];
System.reaction(19).propensity      = J19_ka_19*TNF_TNFR_E*TRADD;

% (R20)
System.reaction(20).educt      = [RIP, TRAF2, TNF_TNFR_TRADD];
System.reaction(20).product      = [TNFRC1];
System.reaction(20).propensity      = J20_ka_20*RIP*TNF_TNFR_TRADD*TRAF2;

% (R21)
System.reaction(21).educt      = [TNFRC1];
System.reaction(21).product      = [TNFRCint1];
System.reaction(21).propensity      = J21_ka_21*TNFRC1;

% (R22)
System.reaction(22).educt      = [TNFRCint1];
System.reaction(22).product      = [RIP, TRAF2, TNFRCint2];
System.reaction(22).propensity      = J22_ka_22*TNFRCint1;

% (R23)
System.reaction(23).educt      = [FADD, TNFRCint2];
System.reaction(23).product      = [TNFRCint3];
System.reaction(23).propensity      = FADD^2*J23_ka_23*TNFRCint2;

% (R24)
System.reaction(24).educt      = [TNFRCint3];
System.reaction(24).product      = [TNFRC2];
System.reaction(24).propensity      = J24_ka_24*TNFRCint3;

% (R25)
System.reaction(25).educt      = [TNFRC2, FLIP];
System.reaction(25).product      = [TNFRC2_FLIP];
System.reaction(25).propensity      = FLIP*J25_ka_25*TNFRC2;

% (R26)
System.reaction(26).educt      = [FLIP, TNFRC2_FLIP];
System.reaction(26).product      = [TNFRC2_FLIP_FLIP];
System.reaction(26).propensity      = FLIP*J26_ka_26*TNFRC2_FLIP;

% (R27)
System.reaction(27).educt      = [TNFRC2, pCasp8];
System.reaction(27).product      = [TNFRC2_pCasp8];
System.reaction(27).propensity      = J27_ka_27*TNFRC2*pCasp8;

% (R28)
System.reaction(28).educt      = [TNFRC2_pCasp8, pCasp8];
System.reaction(28).product      = [TNFRC2_pCasp8_pCasp8];
System.reaction(28).propensity      = J28_ka_28*TNFRC2_pCasp8*pCasp8;

% (R29)
System.reaction(29).educt      = [TNFRC2_pCasp8_pCasp8];
System.reaction(29).product      = [TNFRC2, Casp8];
System.reaction(29).propensity      = J29_ka_29*TNFRC2_pCasp8_pCasp8;

% (R30)
System.reaction(30).educt      = [FLIP, TNFRC2_pCasp8];
System.reaction(30).product      = [TNFRC2_FLIP_pCasp8];
System.reaction(30).propensity      = FLIP*J30_ka_30*TNFRC2_pCasp8;

% (R31)
System.reaction(31).educt      = [TNFRC2_FLIP, pCasp8];
System.reaction(31).product      = [TNFRC2_FLIP_pCasp8];
System.reaction(31).propensity      = J31_ka_31*TNFRC2_FLIP*pCasp8;

% (R32)
System.reaction(32).educt      = [TNFRC2_FLIP_pCasp8];
System.reaction(32).product      = [TNFRC2, Casp8];
System.reaction(32).propensity      = J32_ka_32*TNFRC2_FLIP_pCasp8;

% (R33)
System.reaction(33).educt      = [RIP, TRAF2, TNFRC2_FLIP_pCasp8];
System.reaction(33).product      = [TNFRC2_FLIP_pCasp8_RIP_TRAF2];
System.reaction(33).propensity      = J33_ka_33*RIP*TNFRC2_FLIP_pCasp8*TRAF2;

% (R34)
System.reaction(34).educt      = [IKK];
System.reaction(34).product      = [IKKa];
System.reaction(34).propensity      = IKK*J34_ka_34*TNFRC2_FLIP_pCasp8_RIP_TRAF2;

% (R35)
System.reaction(35).educt      = [];
System.reaction(35).product      = [IKK];
System.reaction(35).propensity      = J35_ka_35 - IKK*J35_kd_35;

% (R36)
System.reaction(36).educt      = [];
System.reaction(36).product      = [NFkB];
System.reaction(36).propensity      = J36_ka_36 - J36_kd_36*NFkB;

% (R37)
System.reaction(37).educt      = [];
System.reaction(37).product      = [FLIP];
System.reaction(37).propensity      = J37_ka_37 - FLIP*J37_kd_37;

% (R38)
System.reaction(38).educt      = [];
System.reaction(38).product      = [XIAP];
System.reaction(38).propensity      = J38_ka_38 - J38_kd_38*XIAP;

% (R39)
System.reaction(39).educt      = [];
System.reaction(39).product      = [A20];
System.reaction(39).propensity      = J39_ka_39 - A20*J39_kd_39;

% (R40)
System.reaction(40).educt      = [IKKa];
System.reaction(40).product      = [];
System.reaction(40).propensity      = IKKa*J40_ka_40;

% (R41)
System.reaction(41).educt      = [IkBa_NFkB];
System.reaction(41).product      = [];
System.reaction(41).propensity      = IkBa_NFkB*J41_ka_41;

% (R42)
System.reaction(42).educt      = [NFkB_N];
System.reaction(42).product      = [];
System.reaction(42).propensity      = J42_ka_42*NFkB_N;

% (R43)
System.reaction(43).educt      = [IkBa_mRNA];
System.reaction(43).product      = [];
System.reaction(43).propensity      = IkBa_mRNA*J43_ka_43;

% (R44)
System.reaction(44).educt      = [IkBa];
System.reaction(44).product      = [];
System.reaction(44).propensity      = IkBa*J44_ka_44;

% (R45)
System.reaction(45).educt      = [IkBa_N];
System.reaction(45).product      = [];
System.reaction(45).propensity      = IkBa_N*J45_ka_45;

% (R46)
System.reaction(46).educt      = [IkBa_NFkB_N];
System.reaction(46).product      = [];
System.reaction(46).propensity      = IkBa_NFkB_N*J46_ka_46;

% (R47)
System.reaction(47).educt      = [PIkBa];
System.reaction(47).product      = [];
System.reaction(47).propensity      = J47_ka_47*PIkBa;

% (R48)
System.reaction(48).educt      = [A20_mRNA];
System.reaction(48).product      = [];
System.reaction(48).propensity      = A20_mRNA*J48_ka_48;

% (R49)
System.reaction(49).educt      = [XIAP_mRNA];
System.reaction(49).product      = [];
System.reaction(49).propensity      = J49_ka_49*XIAP_mRNA;

% (R50)
System.reaction(50).educt      = [FLIP_mRNA];
System.reaction(50).product      = [];
System.reaction(50).propensity      = FLIP_mRNA*J50_ka_50;

% (R51)
System.reaction(51).educt      = [IKK];
System.reaction(51).product      = [IKKa];
System.reaction(51).propensity      = IKK*J51_ka_51*TNFRC1;

% (R52)
System.reaction(52).educt      = [IKKa];
System.reaction(52).product      = [IKK];
System.reaction(52).propensity      = IKKa*J52_ka_52;

% (R53)
System.reaction(53).educt      = [TNFRC1];
System.reaction(53).product      = [TRAF2, TNF_TNFR_TRADD];
System.reaction(53).propensity      = A20*J53_ka_53*TNFRC1;

% (R54)
System.reaction(54).educt      = [NFkB, IkBa];
System.reaction(54).product      = [IkBa_NFkB];
System.reaction(54).propensity      = IkBa*J54_ka_54*NFkB;

% (R55)
System.reaction(55).educt      = [IkBa_NFkB];
System.reaction(55).product      = [NFkB, PIkBa];
System.reaction(55).propensity      = IKKa*IkBa_NFkB*J55_ka_55;

% (R56)
System.reaction(56).educt      = [NFkB];
System.reaction(56).product      = [NFkB_N];
System.reaction(56).propensity      = J56_ka_56*NFkB;

% (R57)
System.reaction(57).educt      = [];
System.reaction(57).product      = [IkBa_mRNA];
System.reaction(57).propensity      = J57_ka_57*NFkB_N;

% (R58)
System.reaction(58).educt      = [];
System.reaction(58).product      = [IkBa];
System.reaction(58).propensity      = IkBa_mRNA*J58_ka_58;

% (R59)
System.reaction(59).educt      = [IkBa];
System.reaction(59).product      = [IkBa_N];
System.reaction(59).propensity      = IkBa*J59_ka_59 - IkBa_N*J59_kd_59;

% (R60)
System.reaction(60).educt      = [NFkB_N, IkBa_N];
System.reaction(60).product      = [IkBa_NFkB_N];
System.reaction(60).propensity      = IkBa_N*J60_ka_60*NFkB_N;

% (R61)
System.reaction(61).educt      = [IkBa_NFkB_N];
System.reaction(61).product      = [IkBa_NFkB];
System.reaction(61).propensity      = IkBa_NFkB_N*J61_ka_61;

% (R62)
System.reaction(62).educt      = [];
System.reaction(62).product      = [A20_mRNA];
System.reaction(62).propensity      = J62_ka_62*NFkB_N;

% (R63)
System.reaction(63).educt      = [];
System.reaction(63).product      = [A20];
System.reaction(63).propensity      = A20_mRNA*J63_ka_63;

% (R64)
System.reaction(64).educt      = [];
System.reaction(64).product      = [XIAP_mRNA];
System.reaction(64).propensity      = J64_ka_64*NFkB_N;

% (R65)
System.reaction(65).educt      = [];
System.reaction(65).product      = [XIAP];
System.reaction(65).propensity      = J65_ka_65*XIAP_mRNA;

% (R66)
System.reaction(66).educt      = [];
System.reaction(66).product      = [FLIP_mRNA];
System.reaction(66).propensity      = J66_ka_66*NFkB_N;

% (R67)
System.reaction(67).educt      = [];
System.reaction(67).product      = [FLIP];
System.reaction(67).propensity      = FLIP_mRNA*J67_ka_67;

% (R68)
System.reaction(68).educt      = [];
System.reaction(68).product      = [pCasp8];
System.reaction(68).propensity      = J68_ka_68 - J68_kd_68*pCasp8;

% (R69)
System.reaction(69).educt      = [];
System.reaction(69).product      = [pCasp3];
System.reaction(69).propensity      = J69_ka_69 - J69_kd_69*pCasp3;

% (R70)
System.reaction(70).educt      = [];
System.reaction(70).product      = [pCasp6];
System.reaction(70).propensity      = J70_ka_70 - J70_kd_70*pCasp6;

% (R71)
System.reaction(71).educt      = [Casp8];
System.reaction(71).product      = [];
System.reaction(71).propensity      = Casp8*J71_ka_71;

% (R72)
System.reaction(72).educt      = [Casp3];
System.reaction(72).product      = [];
System.reaction(72).propensity      = Casp3*J72_ka_72;

% (R73)
System.reaction(73).educt      = [Casp6];
System.reaction(73).product      = [];
System.reaction(73).propensity      = Casp6*J73_ka_73;

% (R74)
System.reaction(74).educt      = [XIAP_Casp3];
System.reaction(74).product      = [];
System.reaction(74).propensity      = J74_ka_74*XIAP_Casp3;

% (R75)
System.reaction(75).educt      = [];
System.reaction(75).product      = [BAR];
System.reaction(75).propensity      = J75_ka_75 - BAR*J75_kd_75;

% (R76)
System.reaction(76).educt      = [BAR_Casp8];
System.reaction(76).product      = [];
System.reaction(76).propensity      = BAR_Casp8*J76_ka_76;

% (R77)
System.reaction(77).educt      = [PARP];
System.reaction(77).product      = [];
System.reaction(77).propensity      = J77_ka_77*PARP - J77_kd_77;

% (R78)
System.reaction(78).educt      = [cPARP];
System.reaction(78).product      = [];
System.reaction(78).propensity      = J78_ka_78*cPARP;

% (R79)
System.reaction(79).educt      = [pCasp3];
System.reaction(79).product      = [Casp3];
System.reaction(79).propensity      = Casp8*J79_ka_79*pCasp3;

% (R80)
System.reaction(80).educt      = [pCasp6];
System.reaction(80).product      = [Casp6];
System.reaction(80).propensity      = Casp3*J80_ka_80*pCasp6;

% (R81)
System.reaction(81).educt      = [pCasp8];
System.reaction(81).product      = [Casp8];
System.reaction(81).propensity      = Casp6*J81_ka_81*pCasp8;

% (R82)
System.reaction(82).educt      = [XIAP, Casp3];
System.reaction(82).product      = [XIAP_Casp3];
System.reaction(82).propensity      = Casp3*J82_ka_82*XIAP - J82_kd_82*XIAP_Casp3;

% (R83)
System.reaction(83).educt      = [XIAP];
System.reaction(83).product      = [];
System.reaction(83).propensity      = Casp3*J83_ka_83*XIAP;

% (R84)
System.reaction(84).educt      = [XIAP_Casp3];
System.reaction(84).product      = [XIAP];
System.reaction(84).propensity      = J84_ka_84*XIAP_Casp3;

% (R85)
System.reaction(85).educt      = [RIP];
System.reaction(85).product      = [];
System.reaction(85).propensity      = Casp3*J85_ka_85*RIP;

% (R86)
System.reaction(86).educt      = [FLIP];
System.reaction(86).product      = [];
System.reaction(86).propensity      = Casp3*FLIP*J86_ka_86;

% (R87)
System.reaction(87).educt      = [PARP];
System.reaction(87).product      = [cPARP];
System.reaction(87).propensity      = Casp3*J87_ka_87*PARP;

% (R88)
System.reaction(88).educt      = [BAR, Casp8];
System.reaction(88).product      = [BAR_Casp8];
System.reaction(88).propensity      = BAR*Casp8*J88_ka_88 - BAR_Casp8*J88_kd_88;

System.output.variable = [Casp3];
System.output.function = [Casp3];
System.output.number   = length(System.output.variable);
System.output.name     = {'Casp3'};

System.input.function = [];
System.input.variable = [];
System.input.number = 0;
System.input.name     = {};

