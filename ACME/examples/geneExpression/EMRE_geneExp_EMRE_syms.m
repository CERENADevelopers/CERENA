% function [model] = EMRE_geneExp_EMRE_syms(f0_user)
function [model] = EMRE_geneExp_EMRE_syms(varargin)

% CVODES OPTIONS

model.atol = 1e-8;
model.rtol = 1e-8;
model.maxsteps = 1e4;

% STATES

syms DNA_off DNA_on mRNA Protein COV11_LNA COV21_LNA COV22_LNA COV31_LNA COV32_LNA COV33_LNA COV41_LNA COV42_LNA COV43_LNA COV44_LNA DNA_off_CORR_EMRE DNA_on_CORR_EMRE mRNA_CORR_EMRE Protein_CORR_EMRE

x = [
DNA_off, DNA_on, mRNA, Protein, COV11_LNA, COV21_LNA, COV22_LNA, COV31_LNA, COV32_LNA, COV33_LNA, COV41_LNA, COV42_LNA, COV43_LNA, COV44_LNA, DNA_off_CORR_EMRE, DNA_on_CORR_EMRE, mRNA_CORR_EMRE, Protein_CORR_EMRE ...
];

% PARAMETERS

syms tau_on tau_off k_m gamma_m k_p gamma_p tau_on_p scaleP offsetP r0 

% KAPPA (constant parameters)

syms Omega indmu1 indmu2 indmu3 indmu4 indC1 indC2 indC3 indC4 indC5 indC6 indC7 indC8 indC9 indC10 kmu01 kmu02 kmu03 kmu04 kC01 kC02 kC03 kC04 kC05 kC06 kC07 kC08 kC09 kC010 

syms t

p = [tau_on,tau_off,k_m,gamma_m,k_p,gamma_p,tau_on_p,scaleP,offsetP,r0];

k = [Omega,indmu1,indmu2,indmu3,indmu4,indC1,indC2,indC3,indC4,indC5,indC6,indC7,indC8,indC9,indC10,kmu01,kmu02,kmu03,kmu04,kC01,kC02,kC03,kC04,kC05,kC06,kC07,kC08,kC09,kC010];

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
	fC01 = f0_user(5); 
	fC02 = f0_user(6); 
	fC03 = f0_user(7); 
	fC04 = f0_user(8); 
	fC05 = f0_user(9); 
	fC06 = f0_user(10); 
	fC07 = f0_user(11); 
	fC08 = f0_user(12); 
	fC09 = f0_user(13); 
	fC010 = f0_user(14); 
else
	fmu01 = 1; 
	fmu02 = 0; 
	fmu03 = r0; 
	fmu04 = 0; 
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
end
% INPUT 

u = sym.empty(0,0);

% SYSTEM EQUATIONS

xdot = sym(zeros(size(x)));

xdot(1) = DNA_on*tau_off - DNA_off*tau_on - DNA_off*Omega*Protein*tau_on_p;
xdot(2) = DNA_off*tau_on - DNA_on*tau_off + DNA_off*Omega*Protein*tau_on_p;
xdot(3) = DNA_on*k_m - gamma_m*mRNA;
xdot(4) = k_p*mRNA - Protein*gamma_p;
xdot(5) = (DNA_on*tau_off + DNA_off*tau_on - 2*COV11_LNA*Omega*tau_on + 2*COV21_LNA*Omega*tau_off - 2*COV41_LNA*DNA_off*Omega^2*tau_on_p - 2*COV11_LNA*Omega^2*Protein*tau_on_p + DNA_off*Omega*Protein*tau_on_p)/Omega;
xdot(6) = -(DNA_on*Omega*tau_off + DNA_off*Omega*tau_on - COV11_LNA*Omega^2*(tau_on + Omega*Protein*tau_on_p) + COV21_LNA*Omega^2*(tau_on + Omega*Protein*tau_on_p) + COV21_LNA*Omega^2*tau_off - COV22_LNA*Omega^2*tau_off - COV41_LNA*DNA_off*Omega^3*tau_on_p + COV42_LNA*DNA_off*Omega^3*tau_on_p + DNA_off*Omega^2*Protein*tau_on_p)/Omega^2;
xdot(7) = (DNA_on*tau_off + DNA_off*tau_on + 2*COV21_LNA*Omega*tau_on - 2*COV22_LNA*Omega*tau_off + 2*COV42_LNA*DNA_off*Omega^2*tau_on_p + 2*COV21_LNA*Omega^2*Protein*tau_on_p + DNA_off*Omega*Protein*tau_on_p)/Omega;
xdot(8) = COV21_LNA*k_m - COV31_LNA*gamma_m - COV31_LNA*tau_on + COV32_LNA*tau_off - COV43_LNA*DNA_off*Omega*tau_on_p - COV31_LNA*Omega*Protein*tau_on_p;
xdot(9) = COV22_LNA*k_m - COV32_LNA*gamma_m + COV31_LNA*tau_on - COV32_LNA*tau_off + COV43_LNA*DNA_off*Omega*tau_on_p + COV31_LNA*Omega*Protein*tau_on_p;
xdot(10) = (DNA_on*k_m + gamma_m*mRNA - 2*COV33_LNA*Omega*gamma_m + 2*COV32_LNA*Omega*k_m)/Omega;
xdot(11) = COV31_LNA*k_p - COV41_LNA*gamma_p - COV41_LNA*tau_on + COV42_LNA*tau_off - COV44_LNA*DNA_off*Omega*tau_on_p - COV41_LNA*Omega*Protein*tau_on_p;
xdot(12) = COV32_LNA*k_p - COV42_LNA*gamma_p + COV41_LNA*tau_on - COV42_LNA*tau_off + COV44_LNA*DNA_off*Omega*tau_on_p + COV41_LNA*Omega*Protein*tau_on_p;
xdot(13) = COV42_LNA*k_m - COV43_LNA*gamma_p - COV43_LNA*gamma_m + COV33_LNA*k_p;
xdot(14) = (Protein*gamma_p + k_p*mRNA - 2*COV44_LNA*Omega*gamma_p + 2*COV43_LNA*Omega*k_p)/Omega;
xdot(15) = DNA_on_CORR_EMRE*tau_off - DNA_off_CORR_EMRE*tau_on - COV41_LNA*Omega*tau_on_p - DNA_off*Omega*Protein_CORR_EMRE*tau_on_p - DNA_off_CORR_EMRE*Omega*Protein*tau_on_p;
xdot(16) = DNA_off_CORR_EMRE*tau_on - DNA_on_CORR_EMRE*tau_off + COV41_LNA*Omega*tau_on_p + DNA_off*Omega*Protein_CORR_EMRE*tau_on_p + DNA_off_CORR_EMRE*Omega*Protein*tau_on_p;
xdot(17) = DNA_on_CORR_EMRE*k_m - gamma_m*mRNA_CORR_EMRE;
xdot(18) = k_p*mRNA_CORR_EMRE - Protein_CORR_EMRE*gamma_p;
% INITIAL CONDITIONS

x0 = sym(zeros(size(x)));

x0(1) = (indmu1*kmu01 - fmu01*(indmu1 - 1))/Omega;
x0(2) = (indmu2*kmu02 - fmu02*(indmu2 - 1))/Omega;
x0(3) = (indmu3*kmu03 - fmu03*(indmu3 - 1))/Omega;
x0(4) = (indmu4*kmu04 - fmu04*(indmu4 - 1))/Omega;

% OBSERVABLES

y = sym(zeros(16,1));

y(1) = DNA_off + DNA_off_CORR_EMRE;
y(2) = DNA_on + DNA_on_CORR_EMRE;
y(3) = mRNA + mRNA_CORR_EMRE;
y(4) = Protein + Protein_CORR_EMRE;
y(5) = COV11_LNA;
y(6) = COV21_LNA;
y(7) = COV22_LNA;
y(8) = COV31_LNA;
y(9) = COV32_LNA;
y(10) = COV33_LNA;
y(11) = COV41_LNA;
y(12) = COV42_LNA;
y(13) = COV43_LNA;
y(14) = COV44_LNA;
y(15) = offsetP + scaleP*(Protein + Protein_CORR_EMRE);
y(16) = COV44_LNA*scaleP^2;

% SYSTEM STRUCT

model.sym.nmx = 14;
model.sym.x = x;
model.sym.u = u;
model.sym.xdot = xdot;
model.sym.p = p;
model.sym.k = k;
model.sym.x0 = x0;
model.sym.y = y;
% Additional fields for the prespecified length of kappa
model.sym.nk1 = 1;
end