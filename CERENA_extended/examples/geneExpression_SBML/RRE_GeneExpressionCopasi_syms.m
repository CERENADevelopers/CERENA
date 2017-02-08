% function [model] = RRE_GeneExpressionCopasi_syms(f0_user)
function [model] = RRE_GeneExpressionCopasi_syms(varargin)

% CVODES OPTIONS

model.atol = 1e-8;
model.rtol = 1e-8;
model.maxsteps = 1e4;

% STATES

syms DNA_off DNA_on mRNA Protein

x = [
DNA_off, DNA_on, mRNA, Protein ...
];

% PARAMETERS

syms tau_off tau_on k_m gamma_m k_p gamma_p tau_on_p 

% KAPPA (constant parameters)

syms indmu1 indmu2 indmu3 indmu4 indC1 indC2 indC3 indC4 indC5 indC6 indC7 indC8 indC9 indC10 kmu01 kmu02 kmu03 kmu04 kC01 kC02 kC03 kC04 kC05 kC06 kC07 kC08 kC09 kC010 

syms t

p = [tau_off,tau_on,k_m,gamma_m,k_p,gamma_p,tau_on_p];

k = [indmu1,indmu2,indmu3,indmu4,indC1,indC2,indC3,indC4,indC5,indC6,indC7,indC8,indC9,indC10,kmu01,kmu02,kmu03,kmu04,kC01,kC02,kC03,kC04,kC05,kC06,kC07,kC08,kC09,kC010];

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
	fmu03 = 4; 
	fmu04 = 10; 
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

xdot(1) = DNA_on*tau_off - DNA_off*tau_on_p - DNA_off*Protein*tau_on_p;
xdot(2) = DNA_off*tau_on_p - DNA_on*tau_off + DNA_off*Protein*tau_on_p;
xdot(3) = DNA_on*k_m - gamma_m*mRNA;
xdot(4) = k_p*mRNA - Protein*gamma_p;
% INITIAL CONDITIONS

x0 = sym(zeros(size(x)));

x0(1) = indmu1*kmu01 - fmu01*(indmu1 - 1);
x0(2) = indmu2*kmu02 - fmu02*(indmu2 - 1);
x0(3) = indmu3*kmu03 - fmu03*(indmu3 - 1);
x0(4) = indmu4*kmu04 - fmu04*(indmu4 - 1);

% OBSERVABLES

y = sym(zeros(8,1));

y(1) = DNA_off;
y(2) = DNA_on;
y(3) = mRNA;
y(4) = Protein;
y(5) = DNA_off;
y(6) = DNA_on;
y(7) = mRNA;
y(8) = Protein;

% SYSTEM STRUCT

model.sym.nmx = 4;
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