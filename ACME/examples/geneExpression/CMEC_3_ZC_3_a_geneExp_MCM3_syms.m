% function [model] = CMEC_3_ZC_3_a_geneExp_MCM3_syms(f0_user)
function [model] = CMEC_3_ZC_3_a_geneExp_MCM3_syms(varargin)

% CVODES OPTIONS

model.atol = 1e-8;
model.rtol = 1e-8;
model.maxsteps = 1e4;

% STATES

syms p_y_1_0 p_y_0_1 mu_1_y_1_0 mu_2_y_1_0 C_1_1_y_1_0 C_1_2_y_1_0 C_2_2_y_1_0 C_1_1_1_y_1_0 C_1_1_2_y_1_0 C_1_2_2_y_1_0 C_2_2_2_y_1_0 mu_1_y_0_1 mu_2_y_0_1 C_1_1_y_0_1 C_1_2_y_0_1 C_2_2_y_0_1 C_1_1_1_y_0_1 C_1_1_2_y_0_1 C_1_2_2_y_0_1 C_2_2_2_y_0_1

x = [
p_y_1_0, p_y_0_1, mu_1_y_1_0, mu_2_y_1_0, C_1_1_y_1_0, C_1_2_y_1_0, C_2_2_y_1_0, C_1_1_1_y_1_0, C_1_1_2_y_1_0, C_1_2_2_y_1_0, C_2_2_2_y_1_0, mu_1_y_0_1, mu_2_y_0_1, C_1_1_y_0_1, C_1_2_y_0_1, C_2_2_y_0_1, C_1_1_1_y_0_1, C_1_1_2_y_0_1, C_1_2_2_y_0_1, C_2_2_2_y_0_1 ...
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

f = sym(zeros(size(x)));

f(1) = p_y_0_1*tau_off - p_y_1_0*tau_on - mu_2_y_1_0*p_y_1_0*tau_on_p;
f(2) = p_y_1_0*tau_on - p_y_0_1*tau_off + mu_2_y_1_0*p_y_1_0*tau_on_p;
f(3) = p_y_0_1*(mu_1_y_0_1*tau_off - mu_1_y_1_0*tau_off) - C_1_2_y_1_0*p_y_1_0*tau_on_p - gamma_m*mu_1_y_1_0*p_y_1_0;
f(4) = p_y_0_1*(mu_2_y_0_1*tau_off - mu_2_y_1_0*tau_off) - C_2_2_y_1_0*p_y_1_0*tau_on_p - gamma_p*mu_2_y_1_0*p_y_1_0 + k_p*mu_1_y_1_0*p_y_1_0;
f(5) = p_y_0_1*(C_1_1_y_0_1*tau_off - C_1_1_y_1_0*tau_off + mu_1_y_0_1^2*tau_off + mu_1_y_1_0^2*tau_off - 2*mu_1_y_0_1*mu_1_y_1_0*tau_off) - 2*C_1_1_y_1_0*gamma_m*p_y_1_0 - C_1_1_2_y_1_0*p_y_1_0*tau_on_p + gamma_m*mu_1_y_1_0*p_y_1_0;
f(6) = p_y_0_1*(C_1_2_y_0_1*tau_off - C_1_2_y_1_0*tau_off + mu_1_y_0_1*mu_2_y_0_1*tau_off - mu_1_y_0_1*mu_2_y_1_0*tau_off - mu_1_y_1_0*mu_2_y_0_1*tau_off + mu_1_y_1_0*mu_2_y_1_0*tau_off) - C_1_2_y_1_0*gamma_m*p_y_1_0 - C_1_2_y_1_0*gamma_p*p_y_1_0 + C_1_1_y_1_0*k_p*p_y_1_0 - C_1_2_2_y_1_0*p_y_1_0*tau_on_p;
f(7) = p_y_0_1*(C_2_2_y_0_1*tau_off - C_2_2_y_1_0*tau_off + mu_2_y_0_1^2*tau_off + mu_2_y_1_0^2*tau_off - 2*mu_2_y_0_1*mu_2_y_1_0*tau_off) - 2*C_2_2_y_1_0*gamma_p*p_y_1_0 + 2*C_1_2_y_1_0*k_p*p_y_1_0 - C_2_2_2_y_1_0*p_y_1_0*tau_on_p + gamma_p*mu_2_y_1_0*p_y_1_0 + k_p*mu_1_y_1_0*p_y_1_0;
f(8) = p_y_0_1*(C_1_1_1_y_0_1*tau_off - C_1_1_1_y_1_0*tau_off + mu_1_y_0_1^3*tau_off - mu_1_y_1_0^3*tau_off + 3*mu_1_y_0_1*mu_1_y_1_0^2*tau_off - 3*mu_1_y_0_1^2*mu_1_y_1_0*tau_off + 3*C_1_1_y_0_1*mu_1_y_0_1*tau_off - 3*C_1_1_y_0_1*mu_1_y_1_0*tau_off - 3*C_1_1_y_1_0*mu_1_y_0_1*p_y_1_0*tau_off + 3*C_1_1_y_1_0*mu_1_y_1_0*p_y_1_0*tau_off) + 3*C_1_1_y_1_0*gamma_m*p_y_1_0 - 3*C_1_1_1_y_1_0*gamma_m*p_y_1_0 - gamma_m*mu_1_y_1_0*p_y_1_0 + 3*C_1_1_y_1_0*C_1_2_y_1_0*p_y_1_0^2*tau_on_p + 3*C_1_1_y_1_0*gamma_m*mu_1_y_1_0*p_y_1_0^2 - 3*C_1_1_y_1_0*C_1_2_y_1_0*p_y_1_0*tau_on_p - 3*C_1_1_y_1_0*gamma_m*mu_1_y_1_0*p_y_1_0;
f(9) = p_y_0_1*(C_1_1_2_y_0_1*tau_off - C_1_1_2_y_1_0*tau_off + mu_1_y_0_1^2*mu_2_y_0_1*tau_off - mu_1_y_0_1^2*mu_2_y_1_0*tau_off + mu_1_y_1_0^2*mu_2_y_0_1*tau_off - mu_1_y_1_0^2*mu_2_y_1_0*tau_off + C_1_1_y_0_1*mu_2_y_0_1*tau_off + 2*C_1_2_y_0_1*mu_1_y_0_1*tau_off - C_1_1_y_0_1*mu_2_y_1_0*tau_off - 2*C_1_2_y_0_1*mu_1_y_1_0*tau_off - 2*mu_1_y_0_1*mu_1_y_1_0*mu_2_y_0_1*tau_off + 2*mu_1_y_0_1*mu_1_y_1_0*mu_2_y_1_0*tau_off - C_1_1_y_1_0*mu_2_y_0_1*p_y_1_0*tau_off - 2*C_1_2_y_1_0*mu_1_y_0_1*p_y_1_0*tau_off + C_1_1_y_1_0*mu_2_y_1_0*p_y_1_0*tau_off + 2*C_1_2_y_1_0*mu_1_y_1_0*p_y_1_0*tau_off) + 2*C_1_2_y_1_0^2*p_y_1_0^2*tau_on_p + C_1_2_y_1_0*gamma_m*p_y_1_0 - 2*C_1_1_2_y_1_0*gamma_m*p_y_1_0 - C_1_1_2_y_1_0*gamma_p*p_y_1_0 + C_1_1_1_y_1_0*k_p*p_y_1_0 - 2*C_1_2_y_1_0^2*p_y_1_0*tau_on_p + C_1_1_y_1_0*C_2_2_y_1_0*p_y_1_0^2*tau_on_p + 2*C_1_2_y_1_0*gamma_m*mu_1_y_1_0*p_y_1_0^2 + C_1_1_y_1_0*gamma_p*mu_2_y_1_0*p_y_1_0^2 - C_1_1_y_1_0*k_p*mu_1_y_1_0*p_y_1_0^2 - C_1_1_y_1_0*C_2_2_y_1_0*p_y_1_0*tau_on_p - 2*C_1_2_y_1_0*gamma_m*mu_1_y_1_0*p_y_1_0 - C_1_1_y_1_0*gamma_p*mu_2_y_1_0*p_y_1_0 + C_1_1_y_1_0*k_p*mu_1_y_1_0*p_y_1_0;
f(10) = p_y_0_1*(C_1_2_2_y_0_1*tau_off - C_1_2_2_y_1_0*tau_off + mu_1_y_0_1*mu_2_y_0_1^2*tau_off + mu_1_y_0_1*mu_2_y_1_0^2*tau_off - mu_1_y_1_0*mu_2_y_0_1^2*tau_off - mu_1_y_1_0*mu_2_y_1_0^2*tau_off + 2*C_1_2_y_0_1*mu_2_y_0_1*tau_off - 2*C_1_2_y_0_1*mu_2_y_1_0*tau_off + C_2_2_y_0_1*mu_1_y_0_1*tau_off - C_2_2_y_0_1*mu_1_y_1_0*tau_off - 2*mu_1_y_0_1*mu_2_y_0_1*mu_2_y_1_0*tau_off + 2*mu_1_y_1_0*mu_2_y_0_1*mu_2_y_1_0*tau_off - 2*C_1_2_y_1_0*mu_2_y_0_1*p_y_1_0*tau_off + 2*C_1_2_y_1_0*mu_2_y_1_0*p_y_1_0*tau_off - C_2_2_y_1_0*mu_1_y_0_1*p_y_1_0*tau_off + C_2_2_y_1_0*mu_1_y_1_0*p_y_1_0*tau_off) + C_1_2_y_1_0*gamma_p*p_y_1_0 - C_1_2_2_y_1_0*gamma_m*p_y_1_0 - 2*C_1_2_2_y_1_0*gamma_p*p_y_1_0 + C_1_1_y_1_0*k_p*p_y_1_0 + 2*C_1_1_2_y_1_0*k_p*p_y_1_0 + 3*C_1_2_y_1_0*C_2_2_y_1_0*p_y_1_0^2*tau_on_p + 2*C_1_2_y_1_0*gamma_p*mu_2_y_1_0*p_y_1_0^2 + C_2_2_y_1_0*gamma_m*mu_1_y_1_0*p_y_1_0^2 - 2*C_1_2_y_1_0*k_p*mu_1_y_1_0*p_y_1_0^2 - 3*C_1_2_y_1_0*C_2_2_y_1_0*p_y_1_0*tau_on_p - 2*C_1_2_y_1_0*gamma_p*mu_2_y_1_0*p_y_1_0 - C_2_2_y_1_0*gamma_m*mu_1_y_1_0*p_y_1_0 + 2*C_1_2_y_1_0*k_p*mu_1_y_1_0*p_y_1_0;
f(11) = p_y_0_1*(C_2_2_2_y_0_1*tau_off - C_2_2_2_y_1_0*tau_off + mu_2_y_0_1^3*tau_off - mu_2_y_1_0^3*tau_off + 3*mu_2_y_0_1*mu_2_y_1_0^2*tau_off - 3*mu_2_y_0_1^2*mu_2_y_1_0*tau_off + 3*C_2_2_y_0_1*mu_2_y_0_1*tau_off - 3*C_2_2_y_0_1*mu_2_y_1_0*tau_off - 3*C_2_2_y_1_0*mu_2_y_0_1*p_y_1_0*tau_off + 3*C_2_2_y_1_0*mu_2_y_1_0*p_y_1_0*tau_off) + 3*C_2_2_y_1_0^2*p_y_1_0^2*tau_on_p + 3*C_2_2_y_1_0*gamma_p*p_y_1_0 - 3*C_2_2_2_y_1_0*gamma_p*p_y_1_0 + 3*C_1_2_y_1_0*k_p*p_y_1_0 + 3*C_1_2_2_y_1_0*k_p*p_y_1_0 - gamma_p*mu_2_y_1_0*p_y_1_0 + k_p*mu_1_y_1_0*p_y_1_0 - 3*C_2_2_y_1_0^2*p_y_1_0*tau_on_p + 3*C_2_2_y_1_0*gamma_p*mu_2_y_1_0*p_y_1_0^2 - 3*C_2_2_y_1_0*k_p*mu_1_y_1_0*p_y_1_0^2 - 3*C_2_2_y_1_0*gamma_p*mu_2_y_1_0*p_y_1_0 + 3*C_2_2_y_1_0*k_p*mu_1_y_1_0*p_y_1_0;
f(12) = p_y_0_1*(k_m - gamma_m*mu_1_y_0_1) + C_1_2_y_1_0*p_y_1_0*tau_on_p - mu_1_y_0_1*p_y_1_0*tau_on + mu_1_y_1_0*p_y_1_0*tau_on - mu_1_y_0_1*mu_2_y_1_0*p_y_1_0*tau_on_p + mu_1_y_1_0*mu_2_y_1_0*p_y_1_0*tau_on_p;
f(13) = mu_2_y_1_0^2*p_y_1_0*tau_on_p - p_y_0_1*(gamma_p*mu_2_y_0_1 - k_p*mu_1_y_0_1) + C_2_2_y_1_0*p_y_1_0*tau_on_p - mu_2_y_0_1*p_y_1_0*tau_on + mu_2_y_1_0*p_y_1_0*tau_on - mu_2_y_0_1*mu_2_y_1_0*p_y_1_0*tau_on_p;
f(14) = p_y_0_1*(k_m - 2*C_1_1_y_0_1*gamma_m + gamma_m*mu_1_y_0_1) - C_1_1_y_0_1*(p_y_1_0*tau_on + mu_2_y_1_0*p_y_1_0*tau_on_p) + p_y_1_0*tau_on*(mu_1_y_0_1 - mu_1_y_1_0)^2 + C_1_1_y_1_0*p_y_1_0*tau_on + p_y_1_0*tau_on_p*(C_1_1_2_y_1_0 + C_1_1_y_1_0*mu_2_y_1_0) - C_1_2_y_1_0*p_y_1_0*tau_on_p*(2*mu_1_y_0_1 - 2*mu_1_y_1_0) + mu_2_y_1_0*p_y_1_0*tau_on_p*(mu_1_y_0_1 - mu_1_y_1_0)^2;
f(15) = C_1_2_y_1_0*p_y_1_0*tau_on - C_1_2_y_0_1*(p_y_1_0*tau_on + mu_2_y_1_0*p_y_1_0*tau_on_p) - p_y_0_1*(C_1_2_y_0_1*gamma_m + C_1_2_y_0_1*gamma_p - C_1_1_y_0_1*k_p) + p_y_1_0*tau_on_p*(C_1_2_2_y_1_0 + C_1_2_y_1_0*mu_2_y_1_0) - C_1_2_y_1_0*p_y_1_0*tau_on_p*(mu_2_y_0_1 - mu_2_y_1_0) - C_2_2_y_1_0*p_y_1_0*tau_on_p*(mu_1_y_0_1 - mu_1_y_1_0) + p_y_1_0*tau_on*(mu_1_y_0_1 - mu_1_y_1_0)*(mu_2_y_0_1 - mu_2_y_1_0) + mu_2_y_1_0*p_y_1_0*tau_on_p*(mu_1_y_0_1 - mu_1_y_1_0)*(mu_2_y_0_1 - mu_2_y_1_0);
f(16) = p_y_0_1*(2*C_1_2_y_0_1*k_p - 2*C_2_2_y_0_1*gamma_p + gamma_p*mu_2_y_0_1 + k_p*mu_1_y_0_1) - C_2_2_y_0_1*(p_y_1_0*tau_on + mu_2_y_1_0*p_y_1_0*tau_on_p) + p_y_1_0*tau_on*(mu_2_y_0_1 - mu_2_y_1_0)^2 + C_2_2_y_1_0*p_y_1_0*tau_on + p_y_1_0*tau_on_p*(C_2_2_2_y_1_0 + C_2_2_y_1_0*mu_2_y_1_0) - C_2_2_y_1_0*p_y_1_0*tau_on_p*(2*mu_2_y_0_1 - 2*mu_2_y_1_0) + mu_2_y_1_0*p_y_1_0*tau_on_p*(mu_2_y_0_1 - mu_2_y_1_0)^2;
f(17) = mu_1_y_1_0^3*p_y_1_0*tau_on - p_y_0_1*(3*C_1_1_1_y_0_1*gamma_m - 3*C_1_1_y_0_1*gamma_m - k_m - 3*C_1_1_y_0_1*k_m + gamma_m*mu_1_y_0_1 + 3*C_1_1_y_0_1*gamma_m*mu_1_y_0_1 + 3*C_1_1_y_0_1*C_1_2_y_1_0*p_y_1_0*tau_on_p - 3*C_1_1_y_0_1*mu_1_y_0_1*p_y_1_0*tau_on + 3*C_1_1_y_0_1*mu_1_y_1_0*p_y_1_0*tau_on - 3*C_1_1_y_0_1*mu_1_y_0_1*mu_2_y_1_0*p_y_1_0*tau_on_p + 3*C_1_1_y_0_1*mu_1_y_1_0*mu_2_y_1_0*p_y_1_0*tau_on_p) - mu_1_y_0_1^3*p_y_1_0*tau_on - p_y_0_1^2*(3*C_1_1_y_0_1*k_m - 3*C_1_1_y_0_1*gamma_m*mu_1_y_0_1) - C_1_1_1_y_0_1*p_y_1_0*tau_on + C_1_1_1_y_1_0*p_y_1_0*tau_on + 3*C_1_2_y_1_0*mu_1_y_0_1^2*p_y_1_0*tau_on_p + 3*C_1_2_y_1_0*mu_1_y_1_0^2*p_y_1_0*tau_on_p - 3*mu_1_y_0_1*mu_1_y_1_0^2*p_y_1_0*tau_on + 3*mu_1_y_0_1^2*mu_1_y_1_0*p_y_1_0*tau_on - mu_1_y_0_1^3*mu_2_y_1_0*p_y_1_0*tau_on_p + mu_1_y_1_0^3*mu_2_y_1_0*p_y_1_0*tau_on_p + 3*C_1_1_y_1_0*C_1_2_y_1_0*p_y_1_0*tau_on_p - 3*C_1_1_y_1_0*mu_1_y_0_1*p_y_1_0*tau_on + 3*C_1_1_y_1_0*mu_1_y_1_0*p_y_1_0*tau_on - C_1_1_1_y_0_1*mu_2_y_1_0*p_y_1_0*tau_on_p - 3*C_1_1_2_y_1_0*mu_1_y_0_1*p_y_1_0*tau_on_p + C_1_1_1_y_1_0*mu_2_y_1_0*p_y_1_0*tau_on_p + 3*C_1_1_2_y_1_0*mu_1_y_1_0*p_y_1_0*tau_on_p - 3*C_1_1_y_1_0*mu_1_y_0_1*mu_2_y_1_0*p_y_1_0*tau_on_p - 6*C_1_2_y_1_0*mu_1_y_0_1*mu_1_y_1_0*p_y_1_0*tau_on_p + 3*C_1_1_y_1_0*mu_1_y_1_0*mu_2_y_1_0*p_y_1_0*tau_on_p - 3*mu_1_y_0_1*mu_1_y_1_0^2*mu_2_y_1_0*p_y_1_0*tau_on_p + 3*mu_1_y_0_1^2*mu_1_y_1_0*mu_2_y_1_0*p_y_1_0*tau_on_p;
f(18) = C_1_1_2_y_1_0*p_y_1_0*tau_on - p_y_0_1*(2*C_1_1_2_y_0_1*gamma_m - C_1_2_y_0_1*gamma_m + C_1_1_2_y_0_1*gamma_p - 2*C_1_2_y_0_1*k_m - C_1_1_1_y_0_1*k_p + 2*C_1_2_y_0_1*gamma_m*mu_1_y_0_1 + C_1_1_y_0_1*gamma_p*mu_2_y_0_1 - C_1_1_y_0_1*k_p*mu_1_y_0_1 + C_1_1_y_0_1*mu_2_y_1_0^2*p_y_1_0*tau_on_p + 2*C_1_2_y_0_1*C_1_2_y_1_0*p_y_1_0*tau_on_p + C_1_1_y_0_1*C_2_2_y_1_0*p_y_1_0*tau_on_p - C_1_1_y_0_1*mu_2_y_0_1*p_y_1_0*tau_on - 2*C_1_2_y_0_1*mu_1_y_0_1*p_y_1_0*tau_on + C_1_1_y_0_1*mu_2_y_1_0*p_y_1_0*tau_on + 2*C_1_2_y_0_1*mu_1_y_1_0*p_y_1_0*tau_on - C_1_1_y_0_1*mu_2_y_0_1*mu_2_y_1_0*p_y_1_0*tau_on_p - 2*C_1_2_y_0_1*mu_1_y_0_1*mu_2_y_1_0*p_y_1_0*tau_on_p + 2*C_1_2_y_0_1*mu_1_y_1_0*mu_2_y_1_0*p_y_1_0*tau_on_p) - C_1_1_2_y_0_1*p_y_1_0*tau_on - p_y_0_1^2*(2*C_1_2_y_0_1*k_m - 2*C_1_2_y_0_1*gamma_m*mu_1_y_0_1 - C_1_1_y_0_1*gamma_p*mu_2_y_0_1 + C_1_1_y_0_1*k_p*mu_1_y_0_1) + 2*C_1_2_y_1_0^2*p_y_1_0*tau_on_p + C_1_1_y_1_0*mu_2_y_1_0^2*p_y_1_0*tau_on_p + C_2_2_y_1_0*mu_1_y_0_1^2*p_y_1_0*tau_on_p + C_2_2_y_1_0*mu_1_y_1_0^2*p_y_1_0*tau_on_p - mu_1_y_0_1^2*mu_2_y_0_1*p_y_1_0*tau_on + mu_1_y_0_1^2*mu_2_y_1_0*p_y_1_0*tau_on - mu_1_y_1_0^2*mu_2_y_0_1*p_y_1_0*tau_on + mu_1_y_1_0^2*mu_2_y_1_0*p_y_1_0*tau_on + C_1_1_y_1_0*C_2_2_y_1_0*p_y_1_0*tau_on_p + mu_1_y_0_1^2*mu_2_y_1_0^2*p_y_1_0*tau_on_p + mu_1_y_1_0^2*mu_2_y_1_0^2*p_y_1_0*tau_on_p - C_1_1_y_1_0*mu_2_y_0_1*p_y_1_0*tau_on - 2*C_1_2_y_1_0*mu_1_y_0_1*p_y_1_0*tau_on + C_1_1_y_1_0*mu_2_y_1_0*p_y_1_0*tau_on + 2*C_1_2_y_1_0*mu_1_y_1_0*p_y_1_0*tau_on - C_1_1_2_y_0_1*mu_2_y_1_0*p_y_1_0*tau_on_p - C_1_1_2_y_1_0*mu_2_y_0_1*p_y_1_0*tau_on_p + 2*C_1_1_2_y_1_0*mu_2_y_1_0*p_y_1_0*tau_on_p - 2*C_1_2_2_y_1_0*mu_1_y_0_1*p_y_1_0*tau_on_p + 2*C_1_2_2_y_1_0*mu_1_y_1_0*p_y_1_0*tau_on_p + 2*C_1_2_y_1_0*mu_1_y_0_1*mu_2_y_0_1*p_y_1_0*tau_on_p - C_1_1_y_1_0*mu_2_y_0_1*mu_2_y_1_0*p_y_1_0*tau_on_p - 4*C_1_2_y_1_0*mu_1_y_0_1*mu_2_y_1_0*p_y_1_0*tau_on_p - 2*C_1_2_y_1_0*mu_1_y_1_0*mu_2_y_0_1*p_y_1_0*tau_on_p + 4*C_1_2_y_1_0*mu_1_y_1_0*mu_2_y_1_0*p_y_1_0*tau_on_p - 2*C_2_2_y_1_0*mu_1_y_0_1*mu_1_y_1_0*p_y_1_0*tau_on_p + 2*mu_1_y_0_1*mu_1_y_1_0*mu_2_y_0_1*p_y_1_0*tau_on - 2*mu_1_y_0_1*mu_1_y_1_0*mu_2_y_1_0*p_y_1_0*tau_on - 2*mu_1_y_0_1*mu_1_y_1_0*mu_2_y_1_0^2*p_y_1_0*tau_on_p - mu_1_y_0_1^2*mu_2_y_0_1*mu_2_y_1_0*p_y_1_0*tau_on_p - mu_1_y_1_0^2*mu_2_y_0_1*mu_2_y_1_0*p_y_1_0*tau_on_p + 2*mu_1_y_0_1*mu_1_y_1_0*mu_2_y_0_1*mu_2_y_1_0*p_y_1_0*tau_on_p;
f(19) = C_1_2_2_y_1_0*p_y_1_0*tau_on - p_y_0_1*(C_1_2_2_y_0_1*gamma_m - C_1_2_y_0_1*gamma_p + 2*C_1_2_2_y_0_1*gamma_p - C_1_1_y_0_1*k_p - C_2_2_y_0_1*k_m - 2*C_1_1_2_y_0_1*k_p + 2*C_1_2_y_0_1*gamma_p*mu_2_y_0_1 + C_2_2_y_0_1*gamma_m*mu_1_y_0_1 - 2*C_1_2_y_0_1*k_p*mu_1_y_0_1 + 2*C_1_2_y_0_1*mu_2_y_1_0^2*p_y_1_0*tau_on_p + 2*C_1_2_y_0_1*C_2_2_y_1_0*p_y_1_0*tau_on_p + C_1_2_y_1_0*C_2_2_y_0_1*p_y_1_0*tau_on_p - 2*C_1_2_y_0_1*mu_2_y_0_1*p_y_1_0*tau_on + 2*C_1_2_y_0_1*mu_2_y_1_0*p_y_1_0*tau_on - C_2_2_y_0_1*mu_1_y_0_1*p_y_1_0*tau_on + C_2_2_y_0_1*mu_1_y_1_0*p_y_1_0*tau_on - 2*C_1_2_y_0_1*mu_2_y_0_1*mu_2_y_1_0*p_y_1_0*tau_on_p - C_2_2_y_0_1*mu_1_y_0_1*mu_2_y_1_0*p_y_1_0*tau_on_p + C_2_2_y_0_1*mu_1_y_1_0*mu_2_y_1_0*p_y_1_0*tau_on_p) - C_1_2_2_y_0_1*p_y_1_0*tau_on - p_y_0_1^2*(C_2_2_y_0_1*k_m - 2*C_1_2_y_0_1*gamma_p*mu_2_y_0_1 - C_2_2_y_0_1*gamma_m*mu_1_y_0_1 + 2*C_1_2_y_0_1*k_p*mu_1_y_0_1) + C_1_2_y_1_0*mu_2_y_0_1^2*p_y_1_0*tau_on_p + 3*C_1_2_y_1_0*mu_2_y_1_0^2*p_y_1_0*tau_on_p - mu_1_y_0_1*mu_2_y_0_1^2*p_y_1_0*tau_on - mu_1_y_0_1*mu_2_y_1_0^2*p_y_1_0*tau_on + mu_1_y_1_0*mu_2_y_0_1^2*p_y_1_0*tau_on + mu_1_y_1_0*mu_2_y_1_0^2*p_y_1_0*tau_on - mu_1_y_0_1*mu_2_y_1_0^3*p_y_1_0*tau_on_p + mu_1_y_1_0*mu_2_y_1_0^3*p_y_1_0*tau_on_p + 3*C_1_2_y_1_0*C_2_2_y_1_0*p_y_1_0*tau_on_p - 2*C_1_2_y_1_0*mu_2_y_0_1*p_y_1_0*tau_on + 2*C_1_2_y_1_0*mu_2_y_1_0*p_y_1_0*tau_on - C_2_2_y_1_0*mu_1_y_0_1*p_y_1_0*tau_on + C_2_2_y_1_0*mu_1_y_1_0*p_y_1_0*tau_on - C_1_2_2_y_0_1*mu_2_y_1_0*p_y_1_0*tau_on_p - 2*C_1_2_2_y_1_0*mu_2_y_0_1*p_y_1_0*tau_on_p + 3*C_1_2_2_y_1_0*mu_2_y_1_0*p_y_1_0*tau_on_p - C_2_2_2_y_1_0*mu_1_y_0_1*p_y_1_0*tau_on_p + C_2_2_2_y_1_0*mu_1_y_1_0*p_y_1_0*tau_on_p - 4*C_1_2_y_1_0*mu_2_y_0_1*mu_2_y_1_0*p_y_1_0*tau_on_p + 2*C_2_2_y_1_0*mu_1_y_0_1*mu_2_y_0_1*p_y_1_0*tau_on_p - 3*C_2_2_y_1_0*mu_1_y_0_1*mu_2_y_1_0*p_y_1_0*tau_on_p - 2*C_2_2_y_1_0*mu_1_y_1_0*mu_2_y_0_1*p_y_1_0*tau_on_p + 3*C_2_2_y_1_0*mu_1_y_1_0*mu_2_y_1_0*p_y_1_0*tau_on_p + 2*mu_1_y_0_1*mu_2_y_0_1*mu_2_y_1_0*p_y_1_0*tau_on - 2*mu_1_y_1_0*mu_2_y_0_1*mu_2_y_1_0*p_y_1_0*tau_on + 2*mu_1_y_0_1*mu_2_y_0_1*mu_2_y_1_0^2*p_y_1_0*tau_on_p - mu_1_y_0_1*mu_2_y_0_1^2*mu_2_y_1_0*p_y_1_0*tau_on_p - 2*mu_1_y_1_0*mu_2_y_0_1*mu_2_y_1_0^2*p_y_1_0*tau_on_p + mu_1_y_1_0*mu_2_y_0_1^2*mu_2_y_1_0*p_y_1_0*tau_on_p;
f(20) = p_y_0_1*(3*C_2_2_y_0_1*gamma_p - 3*C_2_2_2_y_0_1*gamma_p + 3*C_1_2_y_0_1*k_p + 3*C_1_2_2_y_0_1*k_p - gamma_p*mu_2_y_0_1 + k_p*mu_1_y_0_1 - 3*C_2_2_y_0_1*gamma_p*mu_2_y_0_1 + 3*C_2_2_y_0_1*k_p*mu_1_y_0_1 - 3*C_2_2_y_0_1*mu_2_y_1_0^2*p_y_1_0*tau_on_p - 3*C_2_2_y_0_1*C_2_2_y_1_0*p_y_1_0*tau_on_p + 3*C_2_2_y_0_1*mu_2_y_0_1*p_y_1_0*tau_on - 3*C_2_2_y_0_1*mu_2_y_1_0*p_y_1_0*tau_on + 3*C_2_2_y_0_1*mu_2_y_0_1*mu_2_y_1_0*p_y_1_0*tau_on_p) + p_y_0_1^2*(3*C_2_2_y_0_1*gamma_p*mu_2_y_0_1 - 3*C_2_2_y_0_1*k_p*mu_1_y_0_1) - mu_2_y_0_1^3*p_y_1_0*tau_on + mu_2_y_1_0^3*p_y_1_0*tau_on + mu_2_y_1_0^4*p_y_1_0*tau_on_p - C_2_2_2_y_0_1*p_y_1_0*tau_on + C_2_2_2_y_1_0*p_y_1_0*tau_on + 3*C_2_2_y_1_0^2*p_y_1_0*tau_on_p + 3*C_2_2_y_1_0*mu_2_y_0_1^2*p_y_1_0*tau_on_p + 6*C_2_2_y_1_0*mu_2_y_1_0^2*p_y_1_0*tau_on_p - 3*mu_2_y_0_1*mu_2_y_1_0^2*p_y_1_0*tau_on + 3*mu_2_y_0_1^2*mu_2_y_1_0*p_y_1_0*tau_on - 3*mu_2_y_0_1*mu_2_y_1_0^3*p_y_1_0*tau_on_p - mu_2_y_0_1^3*mu_2_y_1_0*p_y_1_0*tau_on_p + 3*mu_2_y_0_1^2*mu_2_y_1_0^2*p_y_1_0*tau_on_p - 3*C_2_2_y_1_0*mu_2_y_0_1*p_y_1_0*tau_on + 3*C_2_2_y_1_0*mu_2_y_1_0*p_y_1_0*tau_on - C_2_2_2_y_0_1*mu_2_y_1_0*p_y_1_0*tau_on_p - 3*C_2_2_2_y_1_0*mu_2_y_0_1*p_y_1_0*tau_on_p + 4*C_2_2_2_y_1_0*mu_2_y_1_0*p_y_1_0*tau_on_p - 9*C_2_2_y_1_0*mu_2_y_0_1*mu_2_y_1_0*p_y_1_0*tau_on_p;
M = diag(sym([[1], [1], [p_y_1_0], [p_y_1_0], [p_y_1_0], [p_y_1_0], [p_y_1_0], [p_y_1_0], [p_y_1_0], [p_y_1_0], [p_y_1_0], [p_y_0_1], [p_y_0_1], [p_y_0_1], [p_y_0_1], [p_y_0_1], [p_y_0_1], [p_y_0_1], [p_y_0_1], [p_y_0_1]]));
% INITIAL CONDITIONS

x0 = sym(zeros(size(x)));

x0(1) = 562949953365017/562949953421312;
x0(2) = 1/10000000000;
x0(3) = indmu3*kmu03 - fmu03*(indmu3 - 1);
x0(4) = indmu4*kmu04 - fmu04*(indmu4 - 1);
x0(5) = indC8*kC08 - fC08*(indC8 - 1);
x0(6) = indC9*kC09 - fC09*(indC9 - 1);
x0(7) = indC10*kC010 - fC010*(indC10 - 1);
x0(8) = 1/10000000000;
x0(9) = 1/10000000000;
x0(10) = 1/10000000000;
x0(11) = 1/10000000000;
x0(12) = indmu3*kmu03 - fmu03*(indmu3 - 1);
x0(13) = indmu4*kmu04 - fmu04*(indmu4 - 1);
x0(14) = indC8*kC08 - fC08*(indC8 - 1);
x0(15) = indC9*kC09 - fC09*(indC9 - 1);
x0(16) = indC10*kC010 - fC010*(indC10 - 1);
x0(17) = 1/10000000000;
x0(18) = 1/10000000000;
x0(19) = 1/10000000000;
x0(20) = 1/10000000000;
dx0 = sym(zeros(size(x)));
dx0(1) = tau_off/10000000000 - (562949953365017*tau_on)/562949953421312 - (562949953365017*tau_on_p*(indmu4*kmu04 - fmu04*(indmu4 - 1)))/562949953421312;
dx0(2) = (562949953365017*tau_on)/562949953421312 - tau_off/10000000000 + (562949953365017*tau_on_p*(indmu4*kmu04 - fmu04*(indmu4 - 1)))/562949953421312;
dx0(3) = - gamma_m*(indmu3*kmu03 - fmu03*(indmu3 - 1)) - tau_on_p*(indC9*kC09 - fC09*(indC9 - 1));
dx0(4) = k_p*(indmu3*kmu03 - fmu03*(indmu3 - 1)) - gamma_p*(indmu4*kmu04 - fmu04*(indmu4 - 1)) - tau_on_p*(indC10*kC010 - fC010*(indC10 - 1));
dx0(5) = gamma_m*(indmu3*kmu03 - fmu03*(indmu3 - 1)) - 2*gamma_m*(indC8*kC08 - fC08*(indC8 - 1)) - tau_on_p/10000000000;
dx0(6) = k_p*(indC8*kC08 - fC08*(indC8 - 1)) - gamma_m*(indC9*kC09 - fC09*(indC9 - 1)) - gamma_p*(indC9*kC09 - fC09*(indC9 - 1)) - tau_on_p/10000000000;
dx0(7) = gamma_p*(indmu4*kmu04 - fmu04*(indmu4 - 1)) - 2*gamma_p*(indC10*kC010 - fC010*(indC10 - 1)) - tau_on_p/10000000000 + 2*k_p*(indC9*kC09 - fC09*(indC9 - 1)) + k_p*(indmu3*kmu03 - fmu03*(indmu3 - 1));
dx0(8) = 3*gamma_m*(indC8*kC08 - fC08*(indC8 - 1)) - (3*gamma_m)/10000000000 - gamma_m*(indmu3*kmu03 - fmu03*(indmu3 - 1)) - (168885*gamma_m*(indC8*kC08 - fC08*(indC8 - 1))*(indmu3*kmu03 - fmu03*(indmu3 - 1)))/562949953421312 - (168885*tau_on_p*(indC8*kC08 - fC08*(indC8 - 1))*(indC9*kC09 - fC09*(indC9 - 1)))/562949953421312;
dx0(9) = k_p/10000000000 - gamma_p/10000000000 - gamma_m/5000000000 - (56295*tau_on_p*(indC9*kC09 - fC09*(indC9 - 1))^2)/281474976710656 + gamma_m*(indC9*kC09 - fC09*(indC9 - 1)) - (56295*gamma_m*(indC9*kC09 - fC09*(indC9 - 1))*(indmu3*kmu03 - fmu03*(indmu3 - 1)))/281474976710656 - (56295*gamma_p*(indC8*kC08 - fC08*(indC8 - 1))*(indmu4*kmu04 - fmu04*(indmu4 - 1)))/562949953421312 + (56295*k_p*(indC8*kC08 - fC08*(indC8 - 1))*(indmu3*kmu03 - fmu03*(indmu3 - 1)))/562949953421312 - (56295*tau_on_p*(indC8*kC08 - fC08*(indC8 - 1))*(indC10*kC010 - fC010*(indC10 - 1)))/562949953421312;
dx0(10) = k_p/5000000000 - gamma_p/5000000000 - gamma_m/10000000000 + gamma_p*(indC9*kC09 - fC09*(indC9 - 1)) + k_p*(indC8*kC08 - fC08*(indC8 - 1)) - (56295*gamma_m*(indC10*kC010 - fC010*(indC10 - 1))*(indmu3*kmu03 - fmu03*(indmu3 - 1)))/562949953421312 - (56295*gamma_p*(indC9*kC09 - fC09*(indC9 - 1))*(indmu4*kmu04 - fmu04*(indmu4 - 1)))/281474976710656 + (56295*k_p*(indC9*kC09 - fC09*(indC9 - 1))*(indmu3*kmu03 - fmu03*(indmu3 - 1)))/281474976710656 - (168885*tau_on_p*(indC9*kC09 - fC09*(indC9 - 1))*(indC10*kC010 - fC010*(indC10 - 1)))/562949953421312;
dx0(11) = (3*k_p)/10000000000 - (3*gamma_p)/10000000000 - (168885*tau_on_p*(indC10*kC010 - fC010*(indC10 - 1))^2)/562949953421312 + 3*gamma_p*(indC10*kC010 - fC010*(indC10 - 1)) - gamma_p*(indmu4*kmu04 - fmu04*(indmu4 - 1)) + 3*k_p*(indC9*kC09 - fC09*(indC9 - 1)) + k_p*(indmu3*kmu03 - fmu03*(indmu3 - 1)) - (168885*gamma_p*(indC10*kC010 - fC010*(indC10 - 1))*(indmu4*kmu04 - fmu04*(indmu4 - 1)))/562949953421312 + (168885*k_p*(indC10*kC010 - fC010*(indC10 - 1))*(indmu3*kmu03 - fmu03*(indmu3 - 1)))/562949953421312;
dx0(12) = k_m - gamma_m*(indmu3*kmu03 - fmu03*(indmu3 - 1)) + (5497558138330244140625*tau_on_p*(indC9*kC09 - fC09*(indC9 - 1)))/549755813888;
dx0(13) = k_p*(indmu3*kmu03 - fmu03*(indmu3 - 1)) - gamma_p*(indmu4*kmu04 - fmu04*(indmu4 - 1)) + (5497558138330244140625*tau_on_p*(indC10*kC010 - fC010*(indC10 - 1)))/549755813888;
dx0(14) = k_m + (5497558138330244140625*tau_on_p*((indC8*kC08 - fC08*(indC8 - 1))*(indmu4*kmu04 - fmu04*(indmu4 - 1)) + 1/10000000000))/549755813888 - 10000000000*(indC8*kC08 - fC08*(indC8 - 1))*((562949953365017*tau_on)/562949953421312 + (562949953365017*tau_on_p*(indmu4*kmu04 - fmu04*(indmu4 - 1)))/562949953421312) - 2*gamma_m*(indC8*kC08 - fC08*(indC8 - 1)) + gamma_m*(indmu3*kmu03 - fmu03*(indmu3 - 1)) + (5497558138330244140625*tau_on*(indC8*kC08 - fC08*(indC8 - 1)))/549755813888;
dx0(15) = (5497558138330244140625*tau_on_p*((indC9*kC09 - fC09*(indC9 - 1))*(indmu4*kmu04 - fmu04*(indmu4 - 1)) + 1/10000000000))/549755813888 - 10000000000*(indC9*kC09 - fC09*(indC9 - 1))*((562949953365017*tau_on)/562949953421312 + (562949953365017*tau_on_p*(indmu4*kmu04 - fmu04*(indmu4 - 1)))/562949953421312) - gamma_m*(indC9*kC09 - fC09*(indC9 - 1)) - gamma_p*(indC9*kC09 - fC09*(indC9 - 1)) + k_p*(indC8*kC08 - fC08*(indC8 - 1)) + (5497558138330244140625*tau_on*(indC9*kC09 - fC09*(indC9 - 1)))/549755813888;
dx0(16) = (5497558138330244140625*tau_on_p*((indC10*kC010 - fC010*(indC10 - 1))*(indmu4*kmu04 - fmu04*(indmu4 - 1)) + 1/10000000000))/549755813888 - 10000000000*(indC10*kC010 - fC010*(indC10 - 1))*((562949953365017*tau_on)/562949953421312 + (562949953365017*tau_on_p*(indmu4*kmu04 - fmu04*(indmu4 - 1)))/562949953421312) - 2*gamma_p*(indC10*kC010 - fC010*(indC10 - 1)) + gamma_p*(indmu4*kmu04 - fmu04*(indmu4 - 1)) + 2*k_p*(indC9*kC09 - fC09*(indC9 - 1)) + k_p*(indmu3*kmu03 - fmu03*(indmu3 - 1)) + (5497558138330244140625*tau_on*(indC10*kC010 - fC010*(indC10 - 1)))/549755813888;
dx0(17) = k_m - (3*gamma_m)/10000000000 + 3*gamma_m*(indC8*kC08 - fC08*(indC8 - 1)) - gamma_m*(indmu3*kmu03 - fmu03*(indmu3 - 1)) + (29999999997*k_m*(indC8*kC08 - fC08*(indC8 - 1)))/10000000000 - (29999999997*gamma_m*(indC8*kC08 - fC08*(indC8 - 1))*(indmu3*kmu03 - fmu03*(indmu3 - 1)))/10000000000 + (16888498599261660139904949*tau_on_p*(indC8*kC08 - fC08*(indC8 - 1))*(indC9*kC09 - fC09*(indC9 - 1)))/562949953421312;
dx0(18) = k_p/10000000000 - gamma_p/10000000000 - gamma_m/5000000000 + (5629499533087220046634983*tau_on_p*(indC9*kC09 - fC09*(indC9 - 1))^2)/281474976710656 + gamma_m*(indC9*kC09 - fC09*(indC9 - 1)) + (9999999999*k_m*(indC9*kC09 - fC09*(indC9 - 1)))/5000000000 - (9999999999*gamma_m*(indC9*kC09 - fC09*(indC9 - 1))*(indmu3*kmu03 - fmu03*(indmu3 - 1)))/5000000000 - (9999999999*gamma_p*(indC8*kC08 - fC08*(indC8 - 1))*(indmu4*kmu04 - fmu04*(indmu4 - 1)))/10000000000 + (9999999999*k_p*(indC8*kC08 - fC08*(indC8 - 1))*(indmu3*kmu03 - fmu03*(indmu3 - 1)))/10000000000 + (5629499533087220046634983*tau_on_p*(indC8*kC08 - fC08*(indC8 - 1))*(indC10*kC010 - fC010*(indC10 - 1)))/562949953421312;
dx0(19) = k_p/5000000000 - gamma_p/5000000000 - gamma_m/10000000000 + gamma_p*(indC9*kC09 - fC09*(indC9 - 1)) + (9999999999*k_m*(indC10*kC010 - fC010*(indC10 - 1)))/10000000000 + k_p*(indC8*kC08 - fC08*(indC8 - 1)) - (9999999999*gamma_m*(indC10*kC010 - fC010*(indC10 - 1))*(indmu3*kmu03 - fmu03*(indmu3 - 1)))/10000000000 - (9999999999*gamma_p*(indC9*kC09 - fC09*(indC9 - 1))*(indmu4*kmu04 - fmu04*(indmu4 - 1)))/5000000000 + (9999999999*k_p*(indC9*kC09 - fC09*(indC9 - 1))*(indmu3*kmu03 - fmu03*(indmu3 - 1)))/5000000000 + (16888498599261660139904949*tau_on_p*(indC9*kC09 - fC09*(indC9 - 1))*(indC10*kC010 - fC010*(indC10 - 1)))/562949953421312;
dx0(20) = (3*k_p)/10000000000 - (3*gamma_p)/10000000000 + (16888498599261660139904949*tau_on_p*(indC10*kC010 - fC010*(indC10 - 1))^2)/562949953421312 + 3*gamma_p*(indC10*kC010 - fC010*(indC10 - 1)) - gamma_p*(indmu4*kmu04 - fmu04*(indmu4 - 1)) + 3*k_p*(indC9*kC09 - fC09*(indC9 - 1)) + k_p*(indmu3*kmu03 - fmu03*(indmu3 - 1)) - (29999999997*gamma_p*(indC10*kC010 - fC010*(indC10 - 1))*(indmu4*kmu04 - fmu04*(indmu4 - 1)))/10000000000 + (29999999997*k_p*(indC10*kC010 - fC010*(indC10 - 1))*(indmu3*kmu03 - fmu03*(indmu3 - 1)))/10000000000;

% OBSERVABLES

y = sym(zeros(37,1));

y(1) = p_y_1_0;
y(2) = p_y_0_1;
y(3) = mu_1_y_0_1*p_y_0_1 + mu_1_y_1_0*p_y_1_0;
y(4) = mu_2_y_0_1*p_y_0_1 + mu_2_y_1_0*p_y_1_0;
y(5) = p_y_0_1*p_y_1_0^2 + p_y_1_0*(p_y_1_0 - 1)^2;
y(6) = p_y_0_1*p_y_1_0*(p_y_0_1 + p_y_1_0 - 2);
y(7) = p_y_1_0*(p_y_1_0 - 1)*(mu_1_y_0_1*p_y_0_1 - mu_1_y_1_0 + mu_1_y_1_0*p_y_1_0) + p_y_0_1*p_y_1_0*(mu_1_y_0_1*p_y_0_1 - mu_1_y_0_1 + mu_1_y_1_0*p_y_1_0);
y(8) = p_y_1_0*(p_y_1_0 - 1)*(mu_2_y_0_1*p_y_0_1 - mu_2_y_1_0 + mu_2_y_1_0*p_y_1_0) + p_y_0_1*p_y_1_0*(mu_2_y_0_1*p_y_0_1 - mu_2_y_0_1 + mu_2_y_1_0*p_y_1_0);
y(9) = p_y_0_1^2*p_y_1_0 + p_y_0_1*(p_y_0_1 - 1)^2;
y(10) = p_y_0_1*(p_y_0_1 - 1)*(mu_1_y_0_1*p_y_0_1 - mu_1_y_0_1 + mu_1_y_1_0*p_y_1_0) + p_y_0_1*p_y_1_0*(mu_1_y_0_1*p_y_0_1 - mu_1_y_1_0 + mu_1_y_1_0*p_y_1_0);
y(11) = p_y_0_1*(p_y_0_1 - 1)*(mu_2_y_0_1*p_y_0_1 - mu_2_y_0_1 + mu_2_y_1_0*p_y_1_0) + p_y_0_1*p_y_1_0*(mu_2_y_0_1*p_y_0_1 - mu_2_y_1_0 + mu_2_y_1_0*p_y_1_0);
y(12) = p_y_0_1*(C_1_1_y_0_1 + (mu_1_y_0_1*p_y_0_1 - mu_1_y_0_1 + mu_1_y_1_0*p_y_1_0)^2) + p_y_1_0*(C_1_1_y_1_0 + (mu_1_y_0_1*p_y_0_1 - mu_1_y_1_0 + mu_1_y_1_0*p_y_1_0)^2);
y(13) = p_y_0_1*(C_1_2_y_0_1 + (mu_1_y_0_1*p_y_0_1 - mu_1_y_0_1 + mu_1_y_1_0*p_y_1_0)*(mu_2_y_0_1*p_y_0_1 - mu_2_y_0_1 + mu_2_y_1_0*p_y_1_0)) + p_y_1_0*(C_1_2_y_1_0 + (mu_1_y_0_1*p_y_0_1 - mu_1_y_1_0 + mu_1_y_1_0*p_y_1_0)*(mu_2_y_0_1*p_y_0_1 - mu_2_y_1_0 + mu_2_y_1_0*p_y_1_0));
y(14) = p_y_0_1*(C_2_2_y_0_1 + (mu_2_y_0_1*p_y_0_1 - mu_2_y_0_1 + mu_2_y_1_0*p_y_1_0)^2) + p_y_1_0*(C_2_2_y_1_0 + (mu_2_y_0_1*p_y_0_1 - mu_2_y_1_0 + mu_2_y_1_0*p_y_1_0)^2);
y(15) = - p_y_0_1*p_y_1_0^3 - p_y_1_0*(p_y_1_0 - 1)^3;
y(16) = - p_y_0_1*p_y_1_0^2*(p_y_0_1 - 1) - p_y_0_1*p_y_1_0*(p_y_1_0 - 1)^2;
y(17) = - p_y_0_1*p_y_1_0^2*(mu_1_y_0_1*p_y_0_1 - mu_1_y_0_1 + mu_1_y_1_0*p_y_1_0) - p_y_1_0*(p_y_1_0 - 1)^2*(mu_1_y_0_1*p_y_0_1 - mu_1_y_1_0 + mu_1_y_1_0*p_y_1_0);
y(18) = - p_y_0_1*p_y_1_0^2*(mu_2_y_0_1*p_y_0_1 - mu_2_y_0_1 + mu_2_y_1_0*p_y_1_0) - p_y_1_0*(p_y_1_0 - 1)^2*(mu_2_y_0_1*p_y_0_1 - mu_2_y_1_0 + mu_2_y_1_0*p_y_1_0);
y(19) = - p_y_0_1*p_y_1_0*(p_y_0_1 - 1)^2 - p_y_0_1^2*p_y_1_0*(p_y_1_0 - 1);
y(20) = - p_y_0_1*p_y_1_0*(p_y_0_1 - 1)*(mu_1_y_0_1*p_y_0_1 - mu_1_y_0_1 + mu_1_y_1_0*p_y_1_0) - p_y_0_1*p_y_1_0*(p_y_1_0 - 1)*(mu_1_y_0_1*p_y_0_1 - mu_1_y_1_0 + mu_1_y_1_0*p_y_1_0);
y(21) = - p_y_0_1*p_y_1_0*(p_y_0_1 - 1)*(mu_2_y_0_1*p_y_0_1 - mu_2_y_0_1 + mu_2_y_1_0*p_y_1_0) - p_y_0_1*p_y_1_0*(p_y_1_0 - 1)*(mu_2_y_0_1*p_y_0_1 - mu_2_y_1_0 + mu_2_y_1_0*p_y_1_0);
y(22) = - p_y_1_0*(C_1_1_y_1_0 + (mu_1_y_0_1*p_y_0_1 - mu_1_y_1_0 + mu_1_y_1_0*p_y_1_0)^2)*(p_y_1_0 - 1) - p_y_0_1*p_y_1_0*(C_1_1_y_0_1 + (mu_1_y_0_1*p_y_0_1 - mu_1_y_0_1 + mu_1_y_1_0*p_y_1_0)^2);
y(23) = - p_y_0_1*p_y_1_0*(C_1_2_y_0_1 + (mu_1_y_0_1*p_y_0_1 - mu_1_y_0_1 + mu_1_y_1_0*p_y_1_0)*(mu_2_y_0_1*p_y_0_1 - mu_2_y_0_1 + mu_2_y_1_0*p_y_1_0)) - p_y_1_0*(C_1_2_y_1_0 + (mu_1_y_0_1*p_y_0_1 - mu_1_y_1_0 + mu_1_y_1_0*p_y_1_0)*(mu_2_y_0_1*p_y_0_1 - mu_2_y_1_0 + mu_2_y_1_0*p_y_1_0))*(p_y_1_0 - 1);
y(24) = - p_y_1_0*(C_2_2_y_1_0 + (mu_2_y_0_1*p_y_0_1 - mu_2_y_1_0 + mu_2_y_1_0*p_y_1_0)^2)*(p_y_1_0 - 1) - p_y_0_1*p_y_1_0*(C_2_2_y_0_1 + (mu_2_y_0_1*p_y_0_1 - mu_2_y_0_1 + mu_2_y_1_0*p_y_1_0)^2);
y(25) = - p_y_0_1^3*p_y_1_0 - p_y_0_1*(p_y_0_1 - 1)^3;
y(26) = - p_y_0_1^2*p_y_1_0*(mu_1_y_0_1*p_y_0_1 - mu_1_y_1_0 + mu_1_y_1_0*p_y_1_0) - p_y_0_1*(p_y_0_1 - 1)^2*(mu_1_y_0_1*p_y_0_1 - mu_1_y_0_1 + mu_1_y_1_0*p_y_1_0);
y(27) = - p_y_0_1^2*p_y_1_0*(mu_2_y_0_1*p_y_0_1 - mu_2_y_1_0 + mu_2_y_1_0*p_y_1_0) - p_y_0_1*(p_y_0_1 - 1)^2*(mu_2_y_0_1*p_y_0_1 - mu_2_y_0_1 + mu_2_y_1_0*p_y_1_0);
y(28) = - p_y_0_1*(C_1_1_y_0_1 + (mu_1_y_0_1*p_y_0_1 - mu_1_y_0_1 + mu_1_y_1_0*p_y_1_0)^2)*(p_y_0_1 - 1) - p_y_0_1*p_y_1_0*(C_1_1_y_1_0 + (mu_1_y_0_1*p_y_0_1 - mu_1_y_1_0 + mu_1_y_1_0*p_y_1_0)^2);
y(29) = - p_y_0_1*p_y_1_0*(C_1_2_y_1_0 + (mu_1_y_0_1*p_y_0_1 - mu_1_y_1_0 + mu_1_y_1_0*p_y_1_0)*(mu_2_y_0_1*p_y_0_1 - mu_2_y_1_0 + mu_2_y_1_0*p_y_1_0)) - p_y_0_1*(C_1_2_y_0_1 + (mu_1_y_0_1*p_y_0_1 - mu_1_y_0_1 + mu_1_y_1_0*p_y_1_0)*(mu_2_y_0_1*p_y_0_1 - mu_2_y_0_1 + mu_2_y_1_0*p_y_1_0))*(p_y_0_1 - 1);
y(30) = - p_y_0_1*(C_2_2_y_0_1 + (mu_2_y_0_1*p_y_0_1 - mu_2_y_0_1 + mu_2_y_1_0*p_y_1_0)^2)*(p_y_0_1 - 1) - p_y_0_1*p_y_1_0*(C_2_2_y_1_0 + (mu_2_y_0_1*p_y_0_1 - mu_2_y_1_0 + mu_2_y_1_0*p_y_1_0)^2);
y(31) = - p_y_0_1*((mu_1_y_0_1*p_y_0_1 - mu_1_y_0_1 + mu_1_y_1_0*p_y_1_0)^3 - C_1_1_1_y_0_1 + C_1_1_y_0_1*(3*mu_1_y_0_1*p_y_0_1 - 3*mu_1_y_0_1 + 3*mu_1_y_1_0*p_y_1_0)) - p_y_1_0*((mu_1_y_0_1*p_y_0_1 - mu_1_y_1_0 + mu_1_y_1_0*p_y_1_0)^3 - C_1_1_1_y_1_0 + C_1_1_y_1_0*(3*mu_1_y_0_1*p_y_0_1 - 3*mu_1_y_1_0 + 3*mu_1_y_1_0*p_y_1_0));
y(32) = - p_y_0_1*((mu_1_y_0_1*p_y_0_1 - mu_1_y_0_1 + mu_1_y_1_0*p_y_1_0)^2*(mu_2_y_0_1*p_y_0_1 - mu_2_y_0_1 + mu_2_y_1_0*p_y_1_0) - C_1_1_2_y_0_1 + C_1_2_y_0_1*(2*mu_1_y_0_1*p_y_0_1 - 2*mu_1_y_0_1 + 2*mu_1_y_1_0*p_y_1_0) + C_1_1_y_0_1*(mu_2_y_0_1*p_y_0_1 - mu_2_y_0_1 + mu_2_y_1_0*p_y_1_0)) - p_y_1_0*((mu_1_y_0_1*p_y_0_1 - mu_1_y_1_0 + mu_1_y_1_0*p_y_1_0)^2*(mu_2_y_0_1*p_y_0_1 - mu_2_y_1_0 + mu_2_y_1_0*p_y_1_0) - C_1_1_2_y_1_0 + C_1_2_y_1_0*(2*mu_1_y_0_1*p_y_0_1 - 2*mu_1_y_1_0 + 2*mu_1_y_1_0*p_y_1_0) + C_1_1_y_1_0*(mu_2_y_0_1*p_y_0_1 - mu_2_y_1_0 + mu_2_y_1_0*p_y_1_0));
y(33) = - p_y_0_1*((mu_1_y_0_1*p_y_0_1 - mu_1_y_0_1 + mu_1_y_1_0*p_y_1_0)*(mu_2_y_0_1*p_y_0_1 - mu_2_y_0_1 + mu_2_y_1_0*p_y_1_0)^2 - C_1_2_2_y_0_1 + C_1_2_y_0_1*(2*mu_2_y_0_1*p_y_0_1 - 2*mu_2_y_0_1 + 2*mu_2_y_1_0*p_y_1_0) + C_2_2_y_0_1*(mu_1_y_0_1*p_y_0_1 - mu_1_y_0_1 + mu_1_y_1_0*p_y_1_0)) - p_y_1_0*((mu_1_y_0_1*p_y_0_1 - mu_1_y_1_0 + mu_1_y_1_0*p_y_1_0)*(mu_2_y_0_1*p_y_0_1 - mu_2_y_1_0 + mu_2_y_1_0*p_y_1_0)^2 - C_1_2_2_y_1_0 + C_1_2_y_1_0*(2*mu_2_y_0_1*p_y_0_1 - 2*mu_2_y_1_0 + 2*mu_2_y_1_0*p_y_1_0) + C_2_2_y_1_0*(mu_1_y_0_1*p_y_0_1 - mu_1_y_1_0 + mu_1_y_1_0*p_y_1_0));
y(34) = - p_y_0_1*((mu_2_y_0_1*p_y_0_1 - mu_2_y_0_1 + mu_2_y_1_0*p_y_1_0)^3 - C_2_2_2_y_0_1 + C_2_2_y_0_1*(3*mu_2_y_0_1*p_y_0_1 - 3*mu_2_y_0_1 + 3*mu_2_y_1_0*p_y_1_0)) - p_y_1_0*((mu_2_y_0_1*p_y_0_1 - mu_2_y_1_0 + mu_2_y_1_0*p_y_1_0)^3 - C_2_2_2_y_1_0 + C_2_2_y_1_0*(3*mu_2_y_0_1*p_y_0_1 - 3*mu_2_y_1_0 + 3*mu_2_y_1_0*p_y_1_0));
y(35) = offsetP + scaleP*(mu_2_y_0_1*p_y_0_1 + mu_2_y_1_0*p_y_1_0);
y(36) = scaleP^2*(p_y_0_1*(C_2_2_y_0_1 + (mu_2_y_0_1*p_y_0_1 - mu_2_y_0_1 + mu_2_y_1_0*p_y_1_0)^2) + p_y_1_0*(C_2_2_y_1_0 + (mu_2_y_0_1*p_y_0_1 - mu_2_y_1_0 + mu_2_y_1_0*p_y_1_0)^2));
y(37) = -scaleP^3*(p_y_0_1*((mu_2_y_0_1*p_y_0_1 - mu_2_y_0_1 + mu_2_y_1_0*p_y_1_0)^3 - C_2_2_2_y_0_1 + C_2_2_y_0_1*(3*mu_2_y_0_1*p_y_0_1 - 3*mu_2_y_0_1 + 3*mu_2_y_1_0*p_y_1_0)) + p_y_1_0*((mu_2_y_0_1*p_y_0_1 - mu_2_y_1_0 + mu_2_y_1_0*p_y_1_0)^3 - C_2_2_2_y_1_0 + C_2_2_y_1_0*(3*mu_2_y_0_1*p_y_0_1 - 3*mu_2_y_1_0 + 3*mu_2_y_1_0*p_y_1_0)));

% SYSTEM STRUCT

model.sym.nmx = 34;
model.sym.x = x;
model.sym.u = u;
model.sym.f = f;
model.sym.M = M;
model.sym.p = p;
model.sym.k = k;
model.sym.x0 = x0;
model.sym.dx0 = dx0;
model.sym.y = y;
% Additional fields for the prespecified length of kappa
model.sym.nk1 = 1;
end