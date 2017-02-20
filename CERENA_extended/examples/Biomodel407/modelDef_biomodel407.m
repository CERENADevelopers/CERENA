syms   TNFR_E TNF_E TNF_TNFR_E TNFR RIP TRADD TRAF2 FADD TNF_TNFR_TRADD TNFRC1 TNFRCint1 TNFRCint2 TNFRCint3 TNFRC2 FLIP TNFRC2_FLIP TNFRC2_pCasp8 TNFRC2_FLIP_FLIP TNFRC2_pCasp8_pCasp8 TNFRC2_FLIP_pCasp8 TNFRC2_FLIP_pCasp8_RIP_TRAF2 IKK IKKa A20 NFkB IkBa IkBa_NFkB PIkBa NFkB_N IkBa_N IkBa_NFkB_N A20_mRNA IkBa_mRNA XIAP_mRNA FLIP_mRNA BAR XIAP pCasp8 pCasp3 pCasp6 Casp8 Casp3 Casp6 BAR_Casp8 XIAP_Casp3 PARP cPARP
syms  ka_1_J1  ka_2_J2  ka_3_J3  ka_4_J4  kd_4_J4  ka_5_J5  kd_5_J5  ka_6_J6  kd_6_J6  ka_7_J7  kd_7_J7  ka_8_J8  ka_9_J9  ka_10_J10  ka_11_J11  ka_12_J12  ka_13_J13  ka_14_J14  ka_15_J15  ka_16_J16  ka_17_J17  ka_18_J18  kd_18_J18  ka_19_J19  ka_20_J20  ka_21_J21  ka_22_J22  ka_23_J23  ka_24_J24  ka_25_J25  ka_26_J26  ka_27_J27  ka_28_J28  ka_29_J29  ka_30_J30  ka_31_J31  ka_32_J32  ka_33_J33  ka_34_J34  ka_35_J35  kd_35_J35  ka_36_J36  kd_36_J36  ka_37_J37  kd_37_J37  ka_38_J38  kd_38_J38  ka_39_J39  kd_39_J39  ka_40_J40  ka_41_J41  ka_42_J42  ka_43_J43  ka_44_J44  ka_45_J45  ka_46_J46  ka_47_J47  ka_48_J48  ka_49_J49  ka_50_J50  ka_51_J51  ka_52_J52  ka_53_J53  ka_54_J54  ka_55_J55  ka_56_J56  ka_57_J57  ka_58_J58  ka_59_J59  kd_59_J59  ka_60_J60  ka_61_J61  ka_62_J62  ka_63_J63  ka_64_J64  ka_65_J65  ka_66_J66  ka_67_J67  ka_68_J68  kd_68_J68  ka_69_J69  kd_69_J69  ka_70_J70  kd_70_J70  ka_71_J71  ka_72_J72  ka_73_J73  ka_74_J74  ka_75_J75  kd_75_J75  ka_76_J76  ka_77_J77  kd_77_J77  ka_78_J78  ka_79_J79  ka_80_J80  ka_81_J81  ka_82_J82  kd_82_J82  ka_83_J83  ka_84_J84  ka_85_J85  ka_86_J86  ka_87_J87  ka_88_J88  kd_88_J88
syms totCasp3
syms   time
System.time = time;
System.compartments = {'cytoplasm', 'extracellular', 'nucleus'};
System.volumes      = [3.2           1344          1.056];
System.state.variable    = [TNFR_E; TNF_E; TNF_TNFR_E; TNFR; RIP; TRADD; TRAF2; FADD; TNF_TNFR_TRADD; TNFRC1; TNFRCint1; TNFRCint2; TNFRCint3; TNFRC2; FLIP; TNFRC2_FLIP; TNFRC2_pCasp8; TNFRC2_FLIP_FLIP; TNFRC2_pCasp8_pCasp8; TNFRC2_FLIP_pCasp8; TNFRC2_FLIP_pCasp8_RIP_TRAF2; IKK; IKKa; A20; NFkB; IkBa; IkBa_NFkB; PIkBa; NFkB_N; IkBa_N; IkBa_NFkB_N; A20_mRNA; IkBa_mRNA; XIAP_mRNA; FLIP_mRNA; BAR; XIAP; pCasp8; pCasp3; pCasp6; Casp8; Casp3; Casp6; BAR_Casp8; XIAP_Casp3; PARP; cPARP];
System.state.compartment = {'extracellular'; 'extracellular'; 'extracellular'; 'cytoplasm'; 'cytoplasm'; 'cytoplasm'; 'cytoplasm'; 'cytoplasm'; 'cytoplasm'; 'cytoplasm'; 'cytoplasm'; 'cytoplasm'; 'cytoplasm'; 'cytoplasm'; 'cytoplasm'; 'cytoplasm'; 'cytoplasm'; 'cytoplasm'; 'cytoplasm'; 'cytoplasm'; 'cytoplasm'; 'cytoplasm'; 'cytoplasm'; 'cytoplasm'; 'cytoplasm'; 'cytoplasm'; 'cytoplasm'; 'cytoplasm'; 'nucleus'; 'nucleus'; 'nucleus'; 'nucleus'; 'nucleus'; 'nucleus'; 'nucleus'; 'cytoplasm'; 'cytoplasm'; 'cytoplasm'; 'cytoplasm'; 'cytoplasm'; 'cytoplasm'; 'cytoplasm'; 'cytoplasm'; 'cytoplasm'; 'cytoplasm'; 'cytoplasm'; 'cytoplasm'};
System.state.number      = length(System.state.variable);
System.state.type        = {'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; 'moment'; };
System.state.name        = {'TNFR_E'; 'TNF_E'; 'TNF:TNFR_E'; 'TNFR'; 'RIP'; 'TRADD'; 'TRAF2'; 'FADD'; 'TNF:TNFR:TRADD'; 'TNFRC1'; 'TNFRCint1'; 'TNFRCint2'; 'TNFRCint3'; 'TNFRC2'; 'FLIP'; 'TNFRC2:FLIP'; 'TNFRC2:pCasp8'; 'TNFRC2:FLIP:FLIP'; 'TNFRC2:pCasp8:pCasp8'; 'TNFRC2:FLIP:pCasp8'; 'TNFRC2:FLIP:pCasp8:RIP:TRAF2'; 'IKK'; 'IKKa'; 'A20'; 'NFkB'; 'IkBa'; 'IkBa:NFkB'; 'PIkBa'; 'NFkB_N'; 'IkBa_N'; 'IkBa:NFkB_N'; 'A20_mRNA'; 'IkBa_mRNA'; 'XIAP_mRNA'; 'FLIP_mRNA'; 'BAR'; 'XIAP'; 'pCasp8'; 'pCasp3'; 'pCasp6'; 'Casp8'; 'Casp3'; 'Casp6'; 'BAR:Casp8'; 'XIAP:Casp3'; 'PARP'; 'cPARP'};
System.state.xmin        = transpose([0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  021965  0.00013906     0.28789      7.8337         3.2         0.8       0.064           0           0           0           0           0      1.6667           0]);
% System.state.C0          = zeros(System.state.number*(System.state.number+1)/2,1);

mu0_par =  sym('mu0',size(System.state.variable));
C0_par = sym('C0',[System.state.number*(System.state.number+1)/2,1]);
System.state.mu0         = mu0_par;
System.state.C0          = C0_par;

System.parameter.variable = [ ka_1_J1; ka_2_J2; ka_3_J3; ka_4_J4; kd_4_J4; ka_5_J5; kd_5_J5; ka_6_J6; kd_6_J6; ka_7_J7; kd_7_J7; ka_8_J8; ka_9_J9; ka_10_J10; ka_11_J11; ka_12_J12; ka_13_J13; ka_14_J14; ka_15_J15; ka_16_J16; ka_17_J17; ka_18_J18; kd_18_J18; ka_19_J19; ka_20_J20; ka_21_J21; ka_22_J22; ka_23_J23; ka_24_J24; ka_25_J25; ka_26_J26; ka_27_J27; ka_28_J28; ka_29_J29; ka_30_J30; ka_31_J31; ka_32_J32; ka_33_J33; ka_34_J34; ka_35_J35; kd_35_J35; ka_36_J36; kd_36_J36; ka_37_J37; kd_37_J37; ka_38_J38; kd_38_J38; ka_39_J39; kd_39_J39; ka_40_J40; ka_41_J41; ka_42_J42; ka_43_J43; ka_44_J44; ka_45_J45; ka_46_J46; ka_47_J47; ka_48_J48; ka_49_J49; ka_50_J50; ka_51_J51; ka_52_J52; ka_53_J53; ka_54_J54; ka_55_J55; ka_56_J56; ka_57_J57; ka_58_J58; ka_59_J59; kd_59_J59; ka_60_J60; ka_61_J61; ka_62_J62; ka_63_J63; ka_64_J64; ka_65_J65; ka_66_J66; ka_67_J67; ka_68_J68; kd_68_J68; ka_69_J69; kd_69_J69; ka_70_J70; kd_70_J70; ka_71_J71; ka_72_J72; ka_73_J73; ka_74_J74; ka_75_J75; kd_75_J75; ka_76_J76; ka_77_J77; kd_77_J77; ka_78_J78; ka_79_J79; ka_80_J80; ka_81_J81; ka_82_J82; kd_82_J82; ka_83_J83; ka_84_J84; ka_85_J85; ka_86_J86; ka_87_J87; ka_88_J88; kd_88_J88];
System.kappa.variable = [mu0_par;C0_par];

conservedSpecies = [];

System.scaleIndicator = 'microscopic';

System.reaction(1).educt      = [TNFR];
System.reaction(1).product      = [TNFR_E];
System.reaction(1).propensity      = ka_1_J1*TNFR;

System.reaction(2).educt      = [];
System.reaction(2).product      = [TNFR];
System.reaction(2).propensity      = ka_2_J2;

System.reaction(3).educt      = [TNFR_E];
System.reaction(3).product      = [];
System.reaction(3).propensity      = ka_3_J3*TNFR_E;

System.reaction(4).educt      = [];
System.reaction(4).product      = [RIP];
System.reaction(4).propensity      = ka_4_J4;

System.reaction(5).educt      = [RIP];
System.reaction(5).product      = [];
System.reaction(5).propensity      = kd_4_J4*RIP;

System.reaction(6).educt      = [];
System.reaction(6).product      = [TRADD];
System.reaction(6).propensity      = ka_5_J5;

System.reaction(7).educt      = [TRADD];
System.reaction(7).product      = [];
System.reaction(7).propensity      = kd_5_J5*TRADD;

System.reaction(8).educt      = [];
System.reaction(8).product      = [TRAF2];
System.reaction(8).propensity      = ka_6_J6;

System.reaction(9).educt      = [TRAF2];
System.reaction(9).product      = [];
System.reaction(9).propensity      = kd_6_J6*TRAF2;

System.reaction(10).educt      = [];
System.reaction(10).product      = [FADD];
System.reaction(10).propensity      = ka_7_J7;

System.reaction(11).educt      = [FADD];
System.reaction(11).product      = [];
System.reaction(11).propensity      = FADD*kd_7_J7;

System.reaction(12).educt      = [TNF_TNFR_E];
System.reaction(12).product      = [];
System.reaction(12).propensity      = ka_8_J8*TNF_TNFR_E;

System.reaction(13).educt      = [TNF_TNFR_TRADD];
System.reaction(13).product      = [];
System.reaction(13).propensity      = ka_9_J9*TNF_TNFR_TRADD;

System.reaction(14).educt      = [TNFRC1];
System.reaction(14).product      = [];
System.reaction(14).propensity      = ka_10_J10*TNFRC1;

System.reaction(15).educt      = [TNFRC2];
System.reaction(15).product      = [];
System.reaction(15).propensity      = ka_11_J11*TNFRC2;

System.reaction(16).educt      = [TNFRC2_FLIP];
System.reaction(16).product      = [];
System.reaction(16).propensity      = ka_12_J12*TNFRC2_FLIP;

System.reaction(17).educt      = [TNFRC2_FLIP_FLIP];
System.reaction(17).product      = [];
System.reaction(17).propensity      = ka_13_J13*TNFRC2_FLIP_FLIP;

System.reaction(18).educt      = [TNFRC2_pCasp8];
System.reaction(18).product      = [];
System.reaction(18).propensity      = ka_14_J14*TNFRC2_pCasp8;

System.reaction(19).educt      = [TNFRC2_pCasp8_pCasp8];
System.reaction(19).product      = [];
System.reaction(19).propensity      = ka_15_J15*TNFRC2_pCasp8_pCasp8;

System.reaction(20).educt      = [TNFRC2_FLIP_pCasp8];
System.reaction(20).product      = [];
System.reaction(20).propensity      = ka_16_J16*TNFRC2_FLIP_pCasp8;

System.reaction(21).educt      = [TNFRC2_FLIP_pCasp8_RIP_TRAF2];
System.reaction(21).product      = [];
System.reaction(21).propensity      = ka_17_J17*TNFRC2_FLIP_pCasp8_RIP_TRAF2;

System.reaction(22).educt      = [TNFR_E, TNF_E];
System.reaction(22).product      = [TNF_TNFR_E];
System.reaction(22).propensity      = ka_18_J18*TNF_E*TNFR_E;

System.reaction(23).educt      = [TNF_TNFR_E];
System.reaction(23).product      = [TNFR_E, TNF_E];
System.reaction(23).propensity      = kd_18_J18*TNF_TNFR_E;

System.reaction(24).educt      = [TNF_TNFR_E, TRADD];
System.reaction(24).product      = [TNF_TNFR_TRADD];
System.reaction(24).propensity      = ka_19_J19*TNF_TNFR_E*TRADD;

System.reaction(25).educt      = [RIP, TRAF2, TNF_TNFR_TRADD];
System.reaction(25).product      = [TNFRC1];
System.reaction(25).propensity      = ka_20_J20*RIP*TNF_TNFR_TRADD*TRAF2;

System.reaction(26).educt      = [TNFRC1];
System.reaction(26).product      = [TNFRCint1];
System.reaction(26).propensity      = ka_21_J21*TNFRC1;

System.reaction(27).educt      = [TNFRCint1];
System.reaction(27).product      = [RIP, TRAF2, TNFRCint2];
System.reaction(27).propensity      = ka_22_J22*TNFRCint1;

System.reaction(28).educt      = [FADD, FADD, TNFRCint2];
System.reaction(28).product      = [TNFRCint3];
System.reaction(28).propensity      = FADD^2*ka_23_J23*TNFRCint2;

System.reaction(29).educt      = [TNFRCint3];
System.reaction(29).product      = [TNFRC2];
System.reaction(29).propensity      = ka_24_J24*TNFRCint3;

System.reaction(30).educt      = [TNFRC2, FLIP];
System.reaction(30).product      = [TNFRC2_FLIP];
System.reaction(30).propensity      = FLIP*ka_25_J25*TNFRC2;

System.reaction(31).educt      = [FLIP, TNFRC2_FLIP];
System.reaction(31).product      = [TNFRC2_FLIP_FLIP];
System.reaction(31).propensity      = FLIP*ka_26_J26*TNFRC2_FLIP;

System.reaction(32).educt      = [TNFRC2, pCasp8];
System.reaction(32).product      = [TNFRC2_pCasp8];
System.reaction(32).propensity      = ka_27_J27*TNFRC2*pCasp8;

System.reaction(33).educt      = [TNFRC2_pCasp8, pCasp8];
System.reaction(33).product      = [TNFRC2_pCasp8_pCasp8];
System.reaction(33).propensity      = ka_28_J28*TNFRC2_pCasp8*pCasp8;

System.reaction(34).educt      = [TNFRC2_pCasp8_pCasp8];
System.reaction(34).product      = [TNFRC2, Casp8];
System.reaction(34).propensity      = ka_29_J29*TNFRC2_pCasp8_pCasp8;

System.reaction(35).educt      = [FLIP, TNFRC2_pCasp8];
System.reaction(35).product      = [TNFRC2_FLIP_pCasp8];
System.reaction(35).propensity      = FLIP*ka_30_J30*TNFRC2_pCasp8;

System.reaction(36).educt      = [TNFRC2_FLIP, pCasp8];
System.reaction(36).product      = [TNFRC2_FLIP_pCasp8];
System.reaction(36).propensity      = ka_31_J31*TNFRC2_FLIP*pCasp8;

System.reaction(37).educt      = [TNFRC2_FLIP_pCasp8];
System.reaction(37).product      = [TNFRC2, Casp8];
System.reaction(37).propensity      = ka_32_J32*TNFRC2_FLIP_pCasp8;

System.reaction(38).educt      = [RIP, TRAF2, TNFRC2_FLIP_pCasp8];
System.reaction(38).product      = [TNFRC2_FLIP_pCasp8_RIP_TRAF2];
System.reaction(38).propensity      = ka_33_J33*RIP*TNFRC2_FLIP_pCasp8*TRAF2;

System.reaction(39).educt      = [IKK];
System.reaction(39).product      = [IKKa];
System.reaction(39).propensity      = IKK*ka_34_J34*TNFRC2_FLIP_pCasp8_RIP_TRAF2;

System.reaction(40).educt      = [];
System.reaction(40).product      = [IKK];
System.reaction(40).propensity      = ka_35_J35;

System.reaction(41).educt      = [IKK];
System.reaction(41).product      = [];
System.reaction(41).propensity      = IKK*kd_35_J35;

System.reaction(42).educt      = [];
System.reaction(42).product      = [NFkB];
System.reaction(42).propensity      = ka_36_J36;

System.reaction(43).educt      = [NFkB];
System.reaction(43).product      = [];
System.reaction(43).propensity      = kd_36_J36*NFkB;

System.reaction(44).educt      = [];
System.reaction(44).product      = [FLIP];
System.reaction(44).propensity      = ka_37_J37;

System.reaction(45).educt      = [FLIP];
System.reaction(45).product      = [];
System.reaction(45).propensity      = FLIP*kd_37_J37;

System.reaction(46).educt      = [];
System.reaction(46).product      = [XIAP];
System.reaction(46).propensity      = ka_38_J38;

System.reaction(47).educt      = [XIAP];
System.reaction(47).product      = [];
System.reaction(47).propensity      = kd_38_J38*XIAP;

System.reaction(48).educt      = [];
System.reaction(48).product      = [A20];
System.reaction(48).propensity      = ka_39_J39;

System.reaction(49).educt      = [A20];
System.reaction(49).product      = [];
System.reaction(49).propensity      = A20*kd_39_J39;

System.reaction(50).educt      = [IKKa];
System.reaction(50).product      = [];
System.reaction(50).propensity      = IKKa*ka_40_J40;

System.reaction(51).educt      = [IkBa_NFkB];
System.reaction(51).product      = [];
System.reaction(51).propensity      = IkBa_NFkB*ka_41_J41;

System.reaction(52).educt      = [NFkB_N];
System.reaction(52).product      = [];
System.reaction(52).propensity      = ka_42_J42*NFkB_N;

System.reaction(53).educt      = [IkBa_mRNA];
System.reaction(53).product      = [];
System.reaction(53).propensity      = IkBa_mRNA*ka_43_J43;

System.reaction(54).educt      = [IkBa];
System.reaction(54).product      = [];
System.reaction(54).propensity      = IkBa*ka_44_J44;

System.reaction(55).educt      = [IkBa_N];
System.reaction(55).product      = [];
System.reaction(55).propensity      = IkBa_N*ka_45_J45;

System.reaction(56).educt      = [IkBa_NFkB_N];
System.reaction(56).product      = [];
System.reaction(56).propensity      = IkBa_NFkB_N*ka_46_J46;

System.reaction(57).educt      = [PIkBa];
System.reaction(57).product      = [];
System.reaction(57).propensity      = ka_47_J47*PIkBa;

System.reaction(58).educt      = [A20_mRNA];
System.reaction(58).product      = [];
System.reaction(58).propensity      = A20_mRNA*ka_48_J48;

System.reaction(59).educt      = [XIAP_mRNA];
System.reaction(59).product      = [];
System.reaction(59).propensity      = ka_49_J49*XIAP_mRNA;

System.reaction(60).educt      = [FLIP_mRNA];
System.reaction(60).product      = [];
System.reaction(60).propensity      = FLIP_mRNA*ka_50_J50;

System.reaction(61).educt      = [IKK];
System.reaction(61).product      = [IKKa];
System.reaction(61).propensity      = IKK*ka_51_J51*TNFRC1;

System.reaction(62).educt      = [IKKa];
System.reaction(62).product      = [IKK];
System.reaction(62).propensity      = IKKa*ka_52_J52;

System.reaction(63).educt      = [TNFRC1];
System.reaction(63).product      = [TRAF2, TNF_TNFR_TRADD];
System.reaction(63).propensity      = A20*ka_53_J53*TNFRC1;

System.reaction(64).educt      = [NFkB, IkBa];
System.reaction(64).product      = [IkBa_NFkB];
System.reaction(64).propensity      = IkBa*ka_54_J54*NFkB;

System.reaction(65).educt      = [IkBa_NFkB];
System.reaction(65).product      = [NFkB, PIkBa];
System.reaction(65).propensity      = IKKa*IkBa_NFkB*ka_55_J55;

System.reaction(66).educt      = [NFkB];
System.reaction(66).product      = [NFkB_N];
System.reaction(66).propensity      = ka_56_J56*NFkB;

System.reaction(67).educt      = [];
System.reaction(67).product      = [IkBa_mRNA];
System.reaction(67).propensity      = ka_57_J57*NFkB_N;

System.reaction(68).educt      = [];
System.reaction(68).product      = [IkBa];
System.reaction(68).propensity      = IkBa_mRNA*ka_58_J58;

System.reaction(69).educt      = [IkBa];
System.reaction(69).product      = [IkBa_N];
System.reaction(69).propensity      = IkBa*ka_59_J59;

System.reaction(70).educt      = [IkBa_N];
System.reaction(70).product      = [IkBa];
System.reaction(70).propensity      = IkBa_N*kd_59_J59;

System.reaction(71).educt      = [NFkB_N, IkBa_N];
System.reaction(71).product      = [IkBa_NFkB_N];
System.reaction(71).propensity      = IkBa_N*ka_60_J60*NFkB_N;

System.reaction(72).educt      = [IkBa_NFkB_N];
System.reaction(72).product      = [IkBa_NFkB];
System.reaction(72).propensity      = IkBa_NFkB_N*ka_61_J61;

System.reaction(73).educt      = [];
System.reaction(73).product      = [A20_mRNA];
System.reaction(73).propensity      = ka_62_J62*NFkB_N;

System.reaction(74).educt      = [];
System.reaction(74).product      = [A20];
System.reaction(74).propensity      = A20_mRNA*ka_63_J63;

System.reaction(75).educt      = [];
System.reaction(75).product      = [XIAP_mRNA];
System.reaction(75).propensity      = ka_64_J64*NFkB_N;

System.reaction(76).educt      = [];
System.reaction(76).product      = [XIAP];
System.reaction(76).propensity      = ka_65_J65*XIAP_mRNA;

System.reaction(77).educt      = [];
System.reaction(77).product      = [FLIP_mRNA];
System.reaction(77).propensity      = ka_66_J66*NFkB_N;

System.reaction(78).educt      = [];
System.reaction(78).product      = [FLIP];
System.reaction(78).propensity      = FLIP_mRNA*ka_67_J67;

System.reaction(79).educt      = [];
System.reaction(79).product      = [pCasp8];
System.reaction(79).propensity      = ka_68_J68;

System.reaction(80).educt      = [pCasp8];
System.reaction(80).product      = [];
System.reaction(80).propensity      = kd_68_J68*pCasp8;

System.reaction(81).educt      = [];
System.reaction(81).product      = [pCasp3];
System.reaction(81).propensity      = ka_69_J69;

System.reaction(82).educt      = [pCasp3];
System.reaction(82).product      = [];
System.reaction(82).propensity      = kd_69_J69*pCasp3;

System.reaction(83).educt      = [];
System.reaction(83).product      = [pCasp6];
System.reaction(83).propensity      = ka_70_J70;

System.reaction(84).educt      = [pCasp6];
System.reaction(84).product      = [];
System.reaction(84).propensity      = kd_70_J70*pCasp6;

System.reaction(85).educt      = [Casp8];
System.reaction(85).product      = [];
System.reaction(85).propensity      = Casp8*ka_71_J71;

System.reaction(86).educt      = [Casp3];
System.reaction(86).product      = [];
System.reaction(86).propensity      = Casp3*ka_72_J72;

System.reaction(87).educt      = [Casp6];
System.reaction(87).product      = [];
System.reaction(87).propensity      = Casp6*ka_73_J73;

System.reaction(88).educt      = [XIAP_Casp3];
System.reaction(88).product      = [];
System.reaction(88).propensity      = ka_74_J74*XIAP_Casp3;

System.reaction(89).educt      = [];
System.reaction(89).product      = [BAR];
System.reaction(89).propensity      = ka_75_J75;

System.reaction(90).educt      = [BAR];
System.reaction(90).product      = [];
System.reaction(90).propensity      = BAR*kd_75_J75;

System.reaction(91).educt      = [BAR_Casp8];
System.reaction(91).product      = [];
System.reaction(91).propensity      = BAR_Casp8*ka_76_J76;

System.reaction(92).educt      = [PARP];
System.reaction(92).product      = [];
System.reaction(92).propensity      = ka_77_J77*PARP;

System.reaction(93).educt      = [];
System.reaction(93).product      = [PARP];
System.reaction(93).propensity      = kd_77_J77;

System.reaction(94).educt      = [cPARP];
System.reaction(94).product      = [];
System.reaction(94).propensity      = ka_78_J78*cPARP;

System.reaction(95).educt      = [pCasp3];
System.reaction(95).product      = [Casp3];
System.reaction(95).propensity      = Casp8*ka_79_J79*pCasp3;

System.reaction(96).educt      = [pCasp6];
System.reaction(96).product      = [Casp6];
System.reaction(96).propensity      = Casp3*ka_80_J80*pCasp6;

System.reaction(97).educt      = [pCasp8];
System.reaction(97).product      = [Casp8];
System.reaction(97).propensity      = Casp6*ka_81_J81*pCasp8;

System.reaction(98).educt      = [XIAP, Casp3];
System.reaction(98).product      = [XIAP_Casp3];
System.reaction(98).propensity      = Casp3*ka_82_J82*XIAP;

System.reaction(99).educt      = [XIAP_Casp3];
System.reaction(99).product      = [XIAP, Casp3];
System.reaction(99).propensity      = kd_82_J82*XIAP_Casp3;

System.reaction(100).educt      = [XIAP];
System.reaction(100).product      = [];
System.reaction(100).propensity      = Casp3*ka_83_J83*XIAP;

System.reaction(101).educt      = [XIAP_Casp3];
System.reaction(101).product      = [XIAP];
System.reaction(101).propensity      = ka_84_J84*XIAP_Casp3;

System.reaction(102).educt      = [RIP];
System.reaction(102).product      = [];
System.reaction(102).propensity      = Casp3*ka_85_J85*RIP;

System.reaction(103).educt      = [FLIP];
System.reaction(103).product      = [];
System.reaction(103).propensity      = Casp3*FLIP*ka_86_J86;

System.reaction(104).educt      = [PARP];
System.reaction(104).product      = [cPARP];
System.reaction(104).propensity      = Casp3*ka_87_J87*PARP;

System.reaction(105).educt      = [BAR, Casp8];
System.reaction(105).product      = [BAR_Casp8];
System.reaction(105).propensity      = BAR*Casp8*ka_88_J88;

System.reaction(106).educt      = [BAR_Casp8];
System.reaction(106).product      = [BAR, Casp8];
System.reaction(106).propensity      = BAR_Casp8*kd_88_J88;

System.output.variable = [totCasp3];
System.output.function = [pCasp3+Casp3];
System.output.number   = length(System.output.variable);
System.output.name     = {'totCasp3'};

System.input.function = [];
System.input.variable = [];
System.input.number = 0;
System.input.name     = {};

