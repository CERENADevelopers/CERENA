% function [model] = IOS_geneExp_IOS_syms(f0_user)
function [model] = IOS_geneExp_IOS_syms(varargin)

% CVODES OPTIONS

model.atol = 1e-8;
model.rtol = 1e-8;
model.maxsteps = 1e4;

% STATES

syms DNA_off DNA_on mRNA Protein COV11_LNA COV21_LNA COV22_LNA COV31_LNA COV32_LNA COV33_LNA COV41_LNA COV42_LNA COV43_LNA COV44_LNA DNA_off_CORR_EMRE DNA_on_CORR_EMRE mRNA_CORR_EMRE Protein_CORR_EMRE COV11_CORR_IOS COV21_CORR_IOS COV22_CORR_IOS COV31_CORR_IOS COV32_CORR_IOS COV33_CORR_IOS COV41_CORR_IOS COV42_CORR_IOS COV43_CORR_IOS COV44_CORR_IOS SKEW111_IOS SKEW211_IOS SKEW221_IOS SKEW222_IOS SKEW311_IOS SKEW321_IOS SKEW322_IOS SKEW331_IOS SKEW332_IOS SKEW333_IOS SKEW411_IOS SKEW421_IOS SKEW422_IOS SKEW431_IOS SKEW432_IOS SKEW433_IOS SKEW441_IOS SKEW442_IOS SKEW443_IOS SKEW444_IOS

x = [
DNA_off, DNA_on, mRNA, Protein, COV11_LNA, COV21_LNA, COV22_LNA, COV31_LNA, COV32_LNA, COV33_LNA, COV41_LNA, COV42_LNA, COV43_LNA, COV44_LNA, DNA_off_CORR_EMRE, DNA_on_CORR_EMRE, mRNA_CORR_EMRE, Protein_CORR_EMRE, COV11_CORR_IOS, COV21_CORR_IOS, COV22_CORR_IOS, COV31_CORR_IOS, COV32_CORR_IOS, COV33_CORR_IOS, COV41_CORR_IOS, COV42_CORR_IOS, COV43_CORR_IOS, COV44_CORR_IOS, SKEW111_IOS, SKEW211_IOS, SKEW221_IOS, SKEW222_IOS, SKEW311_IOS, SKEW321_IOS, SKEW322_IOS, SKEW331_IOS, SKEW332_IOS, SKEW333_IOS, SKEW411_IOS, SKEW421_IOS, SKEW422_IOS, SKEW431_IOS, SKEW432_IOS, SKEW433_IOS, SKEW441_IOS, SKEW442_IOS, SKEW443_IOS, SKEW444_IOS ...
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
xdot(19) = (DNA_on_CORR_EMRE*Omega*tau_off - 2*COV11_CORR_IOS*Omega^2*(tau_on + Omega*Protein*tau_on_p) + 2*DNA_off_CORR_EMRE*Omega*(tau_on/2 + (Omega*Protein*tau_on_p)/2) + COV41_LNA*Omega^2*tau_on_p + 2*COV21_CORR_IOS*Omega^2*tau_off - 2*Omega^3*SKEW411_IOS*tau_on_p - 2*COV41_CORR_IOS*DNA_off*Omega^3*tau_on_p + DNA_off*Omega^2*Protein_CORR_EMRE*tau_on_p)/Omega^2;
xdot(20) = -(DNA_on_CORR_EMRE*Omega*tau_off - COV11_CORR_IOS*Omega^2*(tau_on + Omega*Protein*tau_on_p) + COV21_CORR_IOS*Omega^2*(tau_on + Omega*Protein*tau_on_p) + 2*DNA_off_CORR_EMRE*Omega*(tau_on/2 + (Omega*Protein*tau_on_p)/2) + COV41_LNA*Omega^2*tau_on_p + COV21_CORR_IOS*Omega^2*tau_off - COV22_CORR_IOS*Omega^2*tau_off - Omega^3*SKEW411_IOS*tau_on_p + Omega^3*SKEW421_IOS*tau_on_p - COV41_CORR_IOS*DNA_off*Omega^3*tau_on_p + COV42_CORR_IOS*DNA_off*Omega^3*tau_on_p + DNA_off*Omega^2*Protein_CORR_EMRE*tau_on_p)/Omega^2;
xdot(21) = (DNA_on_CORR_EMRE*Omega*tau_off + 2*COV21_CORR_IOS*Omega^2*(tau_on + Omega*Protein*tau_on_p) + 2*DNA_off_CORR_EMRE*Omega*(tau_on/2 + (Omega*Protein*tau_on_p)/2) + COV41_LNA*Omega^2*tau_on_p - 2*COV22_CORR_IOS*Omega^2*tau_off + 2*Omega^3*SKEW421_IOS*tau_on_p + 2*COV42_CORR_IOS*DNA_off*Omega^3*tau_on_p + DNA_off*Omega^2*Protein_CORR_EMRE*tau_on_p)/Omega^2;
xdot(22) = COV21_CORR_IOS*k_m - COV31_CORR_IOS*gamma_m - COV31_CORR_IOS*tau_on + COV32_CORR_IOS*tau_off - Omega*SKEW431_IOS*tau_on_p - COV43_CORR_IOS*DNA_off*Omega*tau_on_p - COV31_CORR_IOS*Omega*Protein*tau_on_p;
xdot(23) = COV22_CORR_IOS*k_m - COV32_CORR_IOS*gamma_m + COV31_CORR_IOS*tau_on - COV32_CORR_IOS*tau_off + Omega*SKEW431_IOS*tau_on_p + COV43_CORR_IOS*DNA_off*Omega*tau_on_p + COV31_CORR_IOS*Omega*Protein*tau_on_p;
xdot(24) = (DNA_on_CORR_EMRE*k_m + gamma_m*mRNA_CORR_EMRE - 2*COV33_CORR_IOS*Omega*gamma_m + 2*COV32_CORR_IOS*Omega*k_m)/Omega;
xdot(25) = COV31_CORR_IOS*k_p - COV41_CORR_IOS*gamma_p - COV41_CORR_IOS*tau_on + COV42_CORR_IOS*tau_off - Omega*SKEW441_IOS*tau_on_p - COV44_CORR_IOS*DNA_off*Omega*tau_on_p - COV41_CORR_IOS*Omega*Protein*tau_on_p;
xdot(26) = COV32_CORR_IOS*k_p - COV42_CORR_IOS*gamma_p + COV41_CORR_IOS*tau_on - COV42_CORR_IOS*tau_off + Omega*SKEW441_IOS*tau_on_p + COV44_CORR_IOS*DNA_off*Omega*tau_on_p + COV41_CORR_IOS*Omega*Protein*tau_on_p;
xdot(27) = COV42_CORR_IOS*k_m - COV43_CORR_IOS*gamma_p - COV43_CORR_IOS*gamma_m + COV33_CORR_IOS*k_p;
xdot(28) = (Protein_CORR_EMRE*gamma_p + k_p*mRNA_CORR_EMRE - 2*COV44_CORR_IOS*Omega*gamma_p + 2*COV43_CORR_IOS*Omega*k_p)/Omega;
xdot(29) = (3*DNA_off_CORR_EMRE*Omega*(DNA_on*Omega*tau_off + DNA_off*Omega*tau_on + DNA_off*Omega^2*Protein*tau_on_p) + DNA_on*Omega*tau_off - DNA_off*Omega*tau_on + 3*COV11_LNA*Omega^2*(tau_on + Omega*Protein*tau_on_p) - 3*Omega^3*SKEW111_IOS*(tau_on + Omega*Protein*tau_on_p) + 3*COV21_LNA*Omega^2*tau_off + 3*Omega^3*SKEW211_IOS*tau_off - 9*COV11_LNA*COV41_LNA*Omega^4*tau_on_p + 3*COV41_LNA*DNA_off*Omega^3*tau_on_p - DNA_off*Omega^2*Protein*tau_on_p - 3*DNA_off*Omega^4*SKEW411_IOS*tau_on_p)/Omega^3;
xdot(30) = (DNA_off*Omega*tau_on - DNA_on*Omega*tau_off - 2*COV11_LNA*Omega^2*(tau_on + Omega*Protein*tau_on_p) + COV21_LNA*Omega^2*(tau_on + Omega*Protein*tau_on_p) + DNA_on_CORR_EMRE*Omega^2*(DNA_on*tau_off + DNA_off*tau_on + DNA_off*Omega*Protein*tau_on_p) - 2*DNA_off_CORR_EMRE*Omega^2*(DNA_on*tau_off + DNA_off*tau_on + DNA_off*Omega*Protein*tau_on_p) + Omega^3*SKEW111_IOS*(tau_on + Omega*Protein*tau_on_p) - 2*Omega^3*SKEW211_IOS*(tau_on + Omega*Protein*tau_on_p) - 2*Omega^4*tau_on_p*(COV11_LNA*COV42_LNA + 2*COV21_LNA*COV41_LNA) - 2*COV21_LNA*Omega^2*tau_off + COV22_LNA*Omega^2*tau_off - Omega^3*SKEW211_IOS*tau_off + 2*Omega^3*SKEW221_IOS*tau_off + 3*COV11_LNA*COV41_LNA*Omega^4*tau_on_p - 2*COV41_LNA*DNA_off*Omega^3*tau_on_p + COV42_LNA*DNA_off*Omega^3*tau_on_p + DNA_off*Omega^2*Protein*tau_on_p + DNA_off*Omega^4*SKEW411_IOS*tau_on_p - 2*DNA_off*Omega^4*SKEW421_IOS*tau_on_p)/Omega^3;
xdot(31) = -(DNA_off*Omega*tau_on - DNA_on*Omega*tau_off - COV11_LNA*Omega^2*(tau_on + Omega*Protein*tau_on_p) + 2*COV21_LNA*Omega^2*(tau_on + Omega*Protein*tau_on_p) + 2*DNA_on_CORR_EMRE*Omega^2*(DNA_on*tau_off + DNA_off*tau_on + DNA_off*Omega*Protein*tau_on_p) - DNA_off_CORR_EMRE*Omega^2*(DNA_on*tau_off + DNA_off*tau_on + DNA_off*Omega*Protein*tau_on_p) - 2*Omega^3*SKEW211_IOS*(tau_on + Omega*Protein*tau_on_p) + Omega^3*SKEW221_IOS*(tau_on + Omega*Protein*tau_on_p) - 2*Omega^4*tau_on_p*(COV11_LNA*COV42_LNA + 2*COV21_LNA*COV41_LNA) + Omega^4*tau_on_p*(2*COV21_LNA*COV42_LNA + COV22_LNA*COV41_LNA) - COV21_LNA*Omega^2*tau_off + 2*COV22_LNA*Omega^2*tau_off + 2*Omega^3*SKEW221_IOS*tau_off - Omega^3*SKEW222_IOS*tau_off - COV41_LNA*DNA_off*Omega^3*tau_on_p + 2*COV42_LNA*DNA_off*Omega^3*tau_on_p + DNA_off*Omega^2*Protein*tau_on_p - 2*DNA_off*Omega^4*SKEW421_IOS*tau_on_p + DNA_off*Omega^4*SKEW422_IOS*tau_on_p)/Omega^3;
xdot(32) = (3*tau_on_p*(2*COV21_LNA*COV42_LNA*Omega^4 + COV22_LNA*COV41_LNA*Omega^4) + 3*DNA_on_CORR_EMRE*Omega*(DNA_on*Omega*tau_off + DNA_off*Omega*tau_on + DNA_off*Omega^2*Protein*tau_on_p) - DNA_on*Omega*tau_off + DNA_off*Omega*tau_on + 3*COV21_LNA*Omega^2*(tau_on + Omega*Protein*tau_on_p) + 3*Omega^3*SKEW221_IOS*(tau_on + Omega*Protein*tau_on_p) + 3*COV22_LNA*Omega^2*tau_off - 3*Omega^3*SKEW222_IOS*tau_off + 3*COV42_LNA*DNA_off*Omega^3*tau_on_p + DNA_off*Omega^2*Protein*tau_on_p + 3*DNA_off*Omega^4*SKEW422_IOS*tau_on_p)/Omega^3;
xdot(33) = (COV31_LNA*tau_on + COV32_LNA*tau_off - Omega*SKEW311_IOS*gamma_m + Omega*SKEW211_IOS*k_m - 2*Omega*SKEW311_IOS*tau_on + 2*Omega*SKEW321_IOS*tau_off + DNA_on*mRNA_CORR_EMRE*tau_off + DNA_off*mRNA_CORR_EMRE*tau_on - 2*COV11_LNA*COV43_LNA*Omega^2*tau_on_p - 4*COV31_LNA*COV41_LNA*Omega^2*tau_on_p - 2*DNA_off*Omega^2*SKEW431_IOS*tau_on_p - 2*Omega^2*Protein*SKEW311_IOS*tau_on_p + COV43_LNA*DNA_off*Omega*tau_on_p + COV31_LNA*Omega*Protein*tau_on_p + DNA_off*Omega*Protein*mRNA_CORR_EMRE*tau_on_p)/Omega;
xdot(34) = -(COV31_LNA*tau_on + COV32_LNA*tau_off + Omega*SKEW321_IOS*gamma_m - Omega*SKEW221_IOS*k_m - Omega*SKEW311_IOS*tau_on + Omega*SKEW321_IOS*tau_on + Omega*SKEW321_IOS*tau_off - Omega*SKEW322_IOS*tau_off + DNA_on*mRNA_CORR_EMRE*tau_off + DNA_off*mRNA_CORR_EMRE*tau_on - COV11_LNA*COV43_LNA*Omega^2*tau_on_p + COV21_LNA*COV43_LNA*Omega^2*tau_on_p - 2*COV31_LNA*COV41_LNA*Omega^2*tau_on_p + COV31_LNA*COV42_LNA*Omega^2*tau_on_p + COV32_LNA*COV41_LNA*Omega^2*tau_on_p - DNA_off*Omega^2*SKEW431_IOS*tau_on_p + DNA_off*Omega^2*SKEW432_IOS*tau_on_p - Omega^2*Protein*SKEW311_IOS*tau_on_p + Omega^2*Protein*SKEW321_IOS*tau_on_p + COV43_LNA*DNA_off*Omega*tau_on_p + COV31_LNA*Omega*Protein*tau_on_p + DNA_off*Omega*Protein*mRNA_CORR_EMRE*tau_on_p)/Omega;
xdot(35) = (COV31_LNA*tau_on + COV32_LNA*tau_off - Omega*SKEW322_IOS*gamma_m + Omega*SKEW222_IOS*k_m + 2*Omega*SKEW321_IOS*tau_on - 2*Omega*SKEW322_IOS*tau_off + DNA_on*mRNA_CORR_EMRE*tau_off + DNA_off*mRNA_CORR_EMRE*tau_on + 2*COV21_LNA*COV43_LNA*Omega^2*tau_on_p + 2*COV31_LNA*COV42_LNA*Omega^2*tau_on_p + 2*COV32_LNA*COV41_LNA*Omega^2*tau_on_p + 2*DNA_off*Omega^2*SKEW432_IOS*tau_on_p + 2*Omega^2*Protein*SKEW321_IOS*tau_on_p + COV43_LNA*DNA_off*Omega*tau_on_p + COV31_LNA*Omega*Protein*tau_on_p + DNA_off*Omega*Protein*mRNA_CORR_EMRE*tau_on_p)/Omega;
xdot(36) = (COV31_LNA*gamma_m + COV21_LNA*k_m + DNA_on*DNA_off_CORR_EMRE*k_m - 2*Omega*SKEW331_IOS*gamma_m + 2*Omega*SKEW321_IOS*k_m - Omega*SKEW331_IOS*tau_on + Omega*SKEW332_IOS*tau_off + DNA_off_CORR_EMRE*gamma_m*mRNA - 2*COV31_LNA*COV43_LNA*Omega^2*tau_on_p - COV33_LNA*COV41_LNA*Omega^2*tau_on_p - DNA_off*Omega^2*SKEW433_IOS*tau_on_p - Omega^2*Protein*SKEW331_IOS*tau_on_p)/Omega;
xdot(37) = (COV32_LNA*gamma_m + COV22_LNA*k_m + DNA_on*DNA_on_CORR_EMRE*k_m - 2*Omega*SKEW332_IOS*gamma_m + 2*Omega*SKEW322_IOS*k_m + Omega*SKEW331_IOS*tau_on - Omega*SKEW332_IOS*tau_off + DNA_on_CORR_EMRE*gamma_m*mRNA + 2*COV31_LNA*COV43_LNA*Omega^2*tau_on_p + COV33_LNA*COV41_LNA*Omega^2*tau_on_p + DNA_off*Omega^2*SKEW433_IOS*tau_on_p + Omega^2*Protein*SKEW331_IOS*tau_on_p)/Omega;
xdot(38) = (DNA_on*k_m - gamma_m*mRNA + 3*COV33_LNA*Omega*gamma_m + 3*COV32_LNA*Omega*k_m - 3*Omega^2*SKEW333_IOS*gamma_m + 3*Omega^2*SKEW332_IOS*k_m + 3*DNA_on*Omega*k_m*mRNA_CORR_EMRE + 3*Omega*gamma_m*mRNA*mRNA_CORR_EMRE)/Omega^2;
xdot(39) = (COV41_LNA*tau_on + COV42_LNA*tau_off - 4*COV41_LNA^2*Omega^2*tau_on_p + DNA_on*Protein_CORR_EMRE*tau_off + DNA_off*Protein_CORR_EMRE*tau_on - Omega*SKEW411_IOS*gamma_p + Omega*SKEW311_IOS*k_p - 2*Omega*SKEW411_IOS*tau_on + 2*Omega*SKEW421_IOS*tau_off - 2*COV11_LNA*COV44_LNA*Omega^2*tau_on_p - 2*DNA_off*Omega^2*SKEW441_IOS*tau_on_p - 2*Omega^2*Protein*SKEW411_IOS*tau_on_p + COV44_LNA*DNA_off*Omega*tau_on_p + COV41_LNA*Omega*Protein*tau_on_p + DNA_off*Omega*Protein*Protein_CORR_EMRE*tau_on_p)/Omega;
xdot(40) = -(tau_on_p*(COV21_LNA*COV44_LNA*Omega^4 + 2*COV41_LNA*COV42_LNA*Omega^4) - tau_on_p*(2*COV41_LNA^2*Omega^4 + COV11_LNA*COV44_LNA*Omega^4) + Omega*Protein_CORR_EMRE*(DNA_on*Omega*tau_off + DNA_off*Omega*tau_on + DNA_off*Omega^2*Protein*tau_on_p) + COV41_LNA*Omega^2*(tau_on + Omega*Protein*tau_on_p) - Omega^3*SKEW411_IOS*(tau_on + Omega*Protein*tau_on_p) + Omega^3*SKEW421_IOS*(tau_on + Omega*Protein*tau_on_p) + COV42_LNA*Omega^2*tau_off + Omega^3*SKEW421_IOS*gamma_p - Omega^3*SKEW321_IOS*k_p + Omega^3*SKEW421_IOS*tau_off - Omega^3*SKEW422_IOS*tau_off + COV44_LNA*DNA_off*Omega^3*tau_on_p - DNA_off*Omega^4*SKEW441_IOS*tau_on_p + DNA_off*Omega^4*SKEW442_IOS*tau_on_p)/Omega^3;
xdot(41) = (COV41_LNA*tau_on + COV42_LNA*tau_off + DNA_on*Protein_CORR_EMRE*tau_off + DNA_off*Protein_CORR_EMRE*tau_on - Omega*SKEW422_IOS*gamma_p + Omega*SKEW322_IOS*k_p + 2*Omega*SKEW421_IOS*tau_on - 2*Omega*SKEW422_IOS*tau_off + 2*COV21_LNA*COV44_LNA*Omega^2*tau_on_p + 4*COV41_LNA*COV42_LNA*Omega^2*tau_on_p + 2*DNA_off*Omega^2*SKEW442_IOS*tau_on_p + 2*Omega^2*Protein*SKEW421_IOS*tau_on_p + COV44_LNA*DNA_off*Omega*tau_on_p + COV41_LNA*Omega*Protein*tau_on_p + DNA_off*Omega*Protein*Protein_CORR_EMRE*tau_on_p)/Omega;
xdot(42) = SKEW421_IOS*k_m - SKEW431_IOS*gamma_p - SKEW431_IOS*gamma_m + SKEW331_IOS*k_p - SKEW431_IOS*tau_on + SKEW432_IOS*tau_off - COV31_LNA*COV44_LNA*Omega*tau_on_p - 2*COV41_LNA*COV43_LNA*Omega*tau_on_p - DNA_off*Omega*SKEW443_IOS*tau_on_p - Omega*Protein*SKEW431_IOS*tau_on_p;
xdot(43) = SKEW422_IOS*k_m - SKEW432_IOS*gamma_p - SKEW432_IOS*gamma_m + SKEW332_IOS*k_p + SKEW431_IOS*tau_on - SKEW432_IOS*tau_off + COV31_LNA*COV44_LNA*Omega*tau_on_p + 2*COV41_LNA*COV43_LNA*Omega*tau_on_p + DNA_off*Omega*SKEW443_IOS*tau_on_p + Omega*Protein*SKEW431_IOS*tau_on_p;
xdot(44) = (COV43_LNA*gamma_m + COV42_LNA*k_m + DNA_on*Protein_CORR_EMRE*k_m - 2*Omega*SKEW433_IOS*gamma_m - Omega*SKEW433_IOS*gamma_p + 2*Omega*SKEW432_IOS*k_m + Omega*SKEW333_IOS*k_p + Protein_CORR_EMRE*gamma_m*mRNA)/Omega;
xdot(45) = (COV41_LNA*gamma_p + COV31_LNA*k_p + DNA_off_CORR_EMRE*Protein*gamma_p - 2*Omega*SKEW441_IOS*gamma_p + 2*Omega*SKEW431_IOS*k_p - Omega*SKEW441_IOS*tau_on + Omega*SKEW442_IOS*tau_off + DNA_off_CORR_EMRE*k_p*mRNA - 3*COV41_LNA*COV44_LNA*Omega^2*tau_on_p - DNA_off*Omega^2*SKEW444_IOS*tau_on_p - Omega^2*Protein*SKEW441_IOS*tau_on_p)/Omega;
xdot(46) = (COV42_LNA*gamma_p + COV32_LNA*k_p + DNA_on_CORR_EMRE*Protein*gamma_p - 2*Omega*SKEW442_IOS*gamma_p + 2*Omega*SKEW432_IOS*k_p + Omega*SKEW441_IOS*tau_on - Omega*SKEW442_IOS*tau_off + DNA_on_CORR_EMRE*k_p*mRNA + 3*COV41_LNA*COV44_LNA*Omega^2*tau_on_p + DNA_off*Omega^2*SKEW444_IOS*tau_on_p + Omega^2*Protein*SKEW441_IOS*tau_on_p)/Omega;
xdot(47) = (COV43_LNA*gamma_p + COV33_LNA*k_p - Omega*SKEW443_IOS*gamma_m - 2*Omega*SKEW443_IOS*gamma_p + Omega*SKEW442_IOS*k_m + 2*Omega*SKEW433_IOS*k_p + Protein*gamma_p*mRNA_CORR_EMRE + k_p*mRNA*mRNA_CORR_EMRE)/Omega;
xdot(48) = (k_p*mRNA - Protein*gamma_p + 3*COV44_LNA*Omega*gamma_p + 3*COV43_LNA*Omega*k_p - 3*Omega^2*SKEW444_IOS*gamma_p + 3*Omega^2*SKEW443_IOS*k_p + 3*Omega*Protein*Protein_CORR_EMRE*gamma_p + 3*Omega*Protein_CORR_EMRE*k_p*mRNA)/Omega^2;
% INITIAL CONDITIONS

x0 = sym(zeros(size(x)));

x0(1) = (indmu1*kmu01 - fmu01*(indmu1 - 1))/Omega;
x0(2) = (indmu2*kmu02 - fmu02*(indmu2 - 1))/Omega;
x0(3) = (indmu3*kmu03 - fmu03*(indmu3 - 1))/Omega;
x0(4) = (indmu4*kmu04 - fmu04*(indmu4 - 1))/Omega;

% OBSERVABLES

y = sym(zeros(37,1));

y(1) = DNA_off + DNA_off_CORR_EMRE;
y(2) = DNA_on + DNA_on_CORR_EMRE;
y(3) = mRNA + mRNA_CORR_EMRE;
y(4) = Protein + Protein_CORR_EMRE;
y(5) = COV11_LNA + COV11_CORR_IOS - DNA_off_CORR_EMRE^2;
y(6) = COV21_LNA + COV21_CORR_IOS - DNA_on_CORR_EMRE*DNA_off_CORR_EMRE;
y(7) = COV22_LNA + COV22_CORR_IOS - DNA_on_CORR_EMRE^2;
y(8) = COV31_LNA + COV31_CORR_IOS - DNA_off_CORR_EMRE*mRNA_CORR_EMRE;
y(9) = COV32_LNA + COV32_CORR_IOS - DNA_on_CORR_EMRE*mRNA_CORR_EMRE;
y(10) = COV33_LNA + COV33_CORR_IOS - mRNA_CORR_EMRE^2;
y(11) = COV41_LNA + COV41_CORR_IOS - DNA_off_CORR_EMRE*Protein_CORR_EMRE;
y(12) = COV42_LNA + COV42_CORR_IOS - DNA_on_CORR_EMRE*Protein_CORR_EMRE;
y(13) = COV43_LNA + COV43_CORR_IOS - Protein_CORR_EMRE*mRNA_CORR_EMRE;
y(14) = COV44_LNA + COV44_CORR_IOS - Protein_CORR_EMRE^2;
y(15) = SKEW111_IOS;
y(16) = SKEW211_IOS;
y(17) = SKEW221_IOS;
y(18) = SKEW222_IOS;
y(19) = SKEW311_IOS;
y(20) = SKEW321_IOS;
y(21) = SKEW322_IOS;
y(22) = SKEW331_IOS;
y(23) = SKEW332_IOS;
y(24) = SKEW333_IOS;
y(25) = SKEW411_IOS;
y(26) = SKEW421_IOS;
y(27) = SKEW422_IOS;
y(28) = SKEW431_IOS;
y(29) = SKEW432_IOS;
y(30) = SKEW433_IOS;
y(31) = SKEW441_IOS;
y(32) = SKEW442_IOS;
y(33) = SKEW443_IOS;
y(34) = SKEW444_IOS;
y(35) = offsetP + scaleP*(Protein + Protein_CORR_EMRE);
y(36) = scaleP^2*(COV44_LNA + COV44_CORR_IOS - Protein_CORR_EMRE^2);
y(37) = SKEW444_IOS*scaleP^3;

% SYSTEM STRUCT

model.sym.nmx = 34;
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