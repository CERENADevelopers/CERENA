% function [model] = RRE_JakStat_RRE_syms(f0_user)
function [model] = RRE_JakStat_RRE_syms(varargin)

% CVODES OPTIONS

model.atol = 1e-8;
model.rtol = 1e-8;
model.maxsteps = 1e4;

% STATES

syms STAT pSTAT pSTAT_pSTAT npSTAT_npSTAT nSTAT1 nSTAT2 nSTAT3 nSTAT4 nSTAT5

x = [
STAT, pSTAT, pSTAT_pSTAT, npSTAT_npSTAT, nSTAT1, nSTAT2, nSTAT3, nSTAT4, nSTAT5 ...
];

% PARAMETERS

syms k_1 k_2 k_3 k_4 STAT_total_conc Omega_cyt Omega_nuc offset_pSTAT offset_tSTAT scale_pSTAT scale_tSTAT sp1 sp2 sp3 sp4 sp5 a b 

% KAPPA (constant parameters)

syms indmu1 indmu2 indmu3 indmu4 indmu5 indmu6 indmu7 indmu8 indmu9 indC1 indC2 indC3 indC4 indC5 indC6 indC7 indC8 indC9 indC10 indC11 indC12 indC13 indC14 indC15 indC16 indC17 indC18 indC19 indC20 indC21 indC22 indC23 indC24 indC25 indC26 indC27 indC28 indC29 indC30 indC31 indC32 indC33 indC34 indC35 indC36 indC37 indC38 indC39 indC40 indC41 indC42 indC43 indC44 indC45 kmu01 kmu02 kmu03 kmu04 kmu05 kmu06 kmu07 kmu08 kmu09 kC01 kC02 kC03 kC04 kC05 kC06 kC07 kC08 kC09 kC010 kC011 kC012 kC013 kC014 kC015 kC016 kC017 kC018 kC019 kC020 kC021 kC022 kC023 kC024 kC025 kC026 kC027 kC028 kC029 kC030 kC031 kC032 kC033 kC034 kC035 kC036 kC037 kC038 kC039 kC040 kC041 kC042 kC043 kC044 kC045 

syms t

p = [k_1,k_2,k_3,k_4,STAT_total_conc,Omega_cyt,Omega_nuc,offset_pSTAT,offset_tSTAT,scale_pSTAT,scale_tSTAT,sp1,sp2,sp3,sp4,sp5,a,b];

k = [indmu1,indmu2,indmu3,indmu4,indmu5,indmu6,indmu7,indmu8,indmu9,indC1,indC2,indC3,indC4,indC5,indC6,indC7,indC8,indC9,indC10,indC11,indC12,indC13,indC14,indC15,indC16,indC17,indC18,indC19,indC20,indC21,indC22,indC23,indC24,indC25,indC26,indC27,indC28,indC29,indC30,indC31,indC32,indC33,indC34,indC35,indC36,indC37,indC38,indC39,indC40,indC41,indC42,indC43,indC44,indC45,kmu01,kmu02,kmu03,kmu04,kmu05,kmu06,kmu07,kmu08,kmu09,kC01,kC02,kC03,kC04,kC05,kC06,kC07,kC08,kC09,kC010,kC011,kC012,kC013,kC014,kC015,kC016,kC017,kC018,kC019,kC020,kC021,kC022,kC023,kC024,kC025,kC026,kC027,kC028,kC029,kC030,kC031,kC032,kC033,kC034,kC035,kC036,kC037,kC038,kC039,kC040,kC041,kC042,kC043,kC044,kC045];

if nargin > 0
   f0_user = varargin{1};
   if ~isnumeric(f0_user)
      p_user = setdiff(symvar(f0_user),p);
      % ADDITIONAL PARAMETERS IN INITIAL CONDITIONS
      p = [p,p_user];
   end
	fmu01 = f0_user(1); 
	fmu02 = f0_user(2); 
	fmu03 = f0_user(3); 
	fmu04 = f0_user(4); 
	fmu05 = f0_user(5); 
	fmu06 = f0_user(6); 
	fmu07 = f0_user(7); 
	fmu08 = f0_user(8); 
	fmu09 = f0_user(9); 
	fC01 = f0_user(10); 
	fC02 = f0_user(11); 
	fC03 = f0_user(12); 
	fC04 = f0_user(13); 
	fC05 = f0_user(14); 
	fC06 = f0_user(15); 
	fC07 = f0_user(16); 
	fC08 = f0_user(17); 
	fC09 = f0_user(18); 
	fC010 = f0_user(19); 
	fC011 = f0_user(20); 
	fC012 = f0_user(21); 
	fC013 = f0_user(22); 
	fC014 = f0_user(23); 
	fC015 = f0_user(24); 
	fC016 = f0_user(25); 
	fC017 = f0_user(26); 
	fC018 = f0_user(27); 
	fC019 = f0_user(28); 
	fC020 = f0_user(29); 
	fC021 = f0_user(30); 
	fC022 = f0_user(31); 
	fC023 = f0_user(32); 
	fC024 = f0_user(33); 
	fC025 = f0_user(34); 
	fC026 = f0_user(35); 
	fC027 = f0_user(36); 
	fC028 = f0_user(37); 
	fC029 = f0_user(38); 
	fC030 = f0_user(39); 
	fC031 = f0_user(40); 
	fC032 = f0_user(41); 
	fC033 = f0_user(42); 
	fC034 = f0_user(43); 
	fC035 = f0_user(44); 
	fC036 = f0_user(45); 
	fC037 = f0_user(46); 
	fC038 = f0_user(47); 
	fC039 = f0_user(48); 
	fC040 = f0_user(49); 
	fC041 = f0_user(50); 
	fC042 = f0_user(51); 
	fC043 = f0_user(52); 
	fC044 = f0_user(53); 
	fC045 = f0_user(54); 
else
	fmu01 = Omega_cyt*STAT_total_conc; 
	fmu02 = 0; 
	fmu03 = 0; 
	fmu04 = 0; 
	fmu05 = 0; 
	fmu06 = 0; 
	fmu07 = 0; 
	fmu08 = 0; 
	fmu09 = 0; 
	fC01 = 0; 
	fC02 = 0; 
	fC03 = 0; 
	fC04 = 0; 
	fC05 = 0; 
	fC06 = 0; 
	fC07 = 0; 
	fC08 = 0; 
	fC09 = 0; 
	fC010 = 0; 
	fC011 = 0; 
	fC012 = 0; 
	fC013 = 0; 
	fC014 = 0; 
	fC015 = 0; 
	fC016 = 0; 
	fC017 = 0; 
	fC018 = 0; 
	fC019 = 0; 
	fC020 = 0; 
	fC021 = 0; 
	fC022 = 0; 
	fC023 = 0; 
	fC024 = 0; 
	fC025 = 0; 
	fC026 = 0; 
	fC027 = 0; 
	fC028 = 0; 
	fC029 = 0; 
	fC030 = 0; 
	fC031 = 0; 
	fC032 = 0; 
	fC033 = 0; 
	fC034 = 0; 
	fC035 = 0; 
	fC036 = 0; 
	fC037 = 0; 
	fC038 = 0; 
	fC039 = 0; 
	fC040 = 0; 
	fC041 = 0; 
	fC042 = 0; 
	fC043 = 0; 
	fC044 = 0; 
	fC045 = 0; 
end
% INPUT 

u = sym.empty(0,0);

% SYSTEM EQUATIONS

xdot = sym(zeros(size(x)));

xdot(1) = (Omega_nuc*k_4*nSTAT5 - Omega_cyt*STAT*b*k_1*t^2*exp(-a*t))/Omega_cyt;
xdot(2) = -(2*Omega_cyt*k_2*pSTAT^2 - STAT*STAT_total_conc*b*k_1*t^2*exp(-a*t))/STAT_total_conc;
xdot(3) = -(STAT_total_conc*k_3*pSTAT_pSTAT - Omega_cyt*k_2*pSTAT^2)/STAT_total_conc;
xdot(4) = -(Omega_nuc*k_4*npSTAT_npSTAT - Omega_cyt*k_3*pSTAT_pSTAT)/Omega_nuc;
xdot(5) = -k_4*(nSTAT1 - 2*npSTAT_npSTAT);
xdot(6) = k_4*(nSTAT1 - nSTAT2);
xdot(7) = k_4*(nSTAT2 - nSTAT3);
xdot(8) = k_4*(nSTAT3 - nSTAT4);
xdot(9) = k_4*(nSTAT4 - nSTAT5);
% INITIAL CONDITIONS

x0 = sym(zeros(size(x)));

x0(1) = (indmu1*kmu01 - fmu01*(indmu1 - 1))/Omega_cyt;
x0(2) = (indmu2*kmu02 - fmu02*(indmu2 - 1))/Omega_cyt;
x0(3) = (indmu3*kmu03 - fmu03*(indmu3 - 1))/Omega_cyt;
x0(4) = (indmu4*kmu04 - fmu04*(indmu4 - 1))/Omega_nuc;
x0(5) = (indmu5*kmu05 - fmu05*(indmu5 - 1))/Omega_nuc;
x0(6) = (indmu6*kmu06 - fmu06*(indmu6 - 1))/Omega_nuc;
x0(7) = (indmu7*kmu07 - fmu07*(indmu7 - 1))/Omega_nuc;
x0(8) = (indmu8*kmu08 - fmu08*(indmu8 - 1))/Omega_nuc;
x0(9) = (indmu9*kmu09 - fmu09*(indmu9 - 1))/Omega_nuc;

% OBSERVABLES

y = sym(zeros(11,1));

y(1) = STAT;
y(2) = pSTAT;
y(3) = pSTAT_pSTAT;
y(4) = npSTAT_npSTAT;
y(5) = nSTAT1;
y(6) = nSTAT2;
y(7) = nSTAT3;
y(8) = nSTAT4;
y(9) = nSTAT5;
y(10) = (STAT_total_conc*offset_pSTAT + pSTAT*scale_pSTAT + 2*pSTAT_pSTAT*scale_pSTAT)/STAT_total_conc;
y(11) = (STAT_total_conc*offset_tSTAT + STAT*scale_tSTAT + pSTAT*scale_tSTAT + 2*pSTAT_pSTAT*scale_tSTAT)/STAT_total_conc;

% SYSTEM STRUCT

model.sym.nmx = 9;
model.sym.x = x;
model.sym.u = u;
model.sym.xdot = xdot;
model.sym.p = p;
model.sym.k = k;
model.sym.x0 = x0;
model.sym.y = y;
% Additional fields for the prespecified length of kappa
model.sym.nk1 = 0;
end