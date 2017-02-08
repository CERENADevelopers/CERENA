% function MM = getMM_centered(system,options)
function MM = getMM_centered_greedy(varargin)

%% CHECK / ASSIGN INPUTS
% Check required inputs
if nargin >= 1
    if ~isempty(varargin{1})
        % Assign inputs
        system = varargin{1};
    else
        % Error message:
        error('The first input must be non-empty!');
    end
else
    % Error message:
    error('At least one input is required!');
end

% Options
options.moment_order = 2;
options.moment_closure = 'none';
options.filename = 'simMM';
if nargin >= 2
    options = setdefault(varargin{2},options);
end
if ~isfield(options,'moment_order_max') || isempty(options.moment_order_max)
    options.moment_order_max = options.moment_order + 2;
end
%% INITIALIZATION (1)
n_c = length(system.parameter.variable);
n_kappa = length(system.kappa.variable);

%% INITIALIZATION (2)
% Generation of I-index
I = system.MM.sym.state.order;
n_I = size(I,1);
mu_ind = I(find(sum(I>=1,2) == 1),:);
% Mean
mu  = getMu(mu_ind);
% Moments with order >= 2
% All moments
M = system.MM.sym.state.moments;

% Parameters
c = sym(zeros(n_c,1));
for i = 1:n_c
    c(i) = sym(['c' num2str(i,'%d')]);
end

% Kappa
c_kappa = sym(zeros(n_kappa,1));
for i = 1:n_kappa
    c_kappa(i) = sym(['c_kappa' num2str(i,'%d')]);
end

%% REACTION STOICHIOMETRY
% S_e = system.eductStoichiometry;
% S_p = system.productStoichiometry;
% Overall reaction stoichiometry
M0 = system.MM.sym.state.M0;
H = system.MM.sym.output.function;
%% CONSTRUCTION OF MOMENT EQUATION
dMdt = system.MM.sym.state.dMdt_raw;
disp('dMdt reused!')

dMdt_raw = dMdt;
disp('construction of dMdt done.')
%% EXTRACTING IMPORTANT COVARIANCES
Iind = system.options.Iind;
[~,ind_M] =ismember(Iind,I,'rows');
ind_hoM = sort(setdiff(1:n_I,ind_M));
MAll = M;
IAll = I;
M = M(ind_M);
M0 = M0(ind_M);
dMdt = dMdt(ind_M);
I = I(ind_M,:);
hoM = MAll(ind_hoM);
hoI = IAll(ind_hoM,:);

%% MOMENT CLOSURE
% Determine moments of order > options.moment_order

hoM_used = hoM;
% hoM_closure2 = sym([]);
hoM_closure = sym(zeros(size(hoM_used)));

hoM_used = transpose(hoM_used);
hoM_closure = transpose(hoM_closure);
dMdt = mysubs(dMdt,hoM_used,hoM_closure);
H = mysubs(H,hoM_used,hoM_closure);

% dMdt = simplify(dMdt);
% dMdt = simplify(dMdt);
disp('moment closure done.')
%% EVALUATION OF JACOBIAN
f = dMdt;
if isfield(system,'input')
    f = mysubs(f,system.input.variable,system.input.function);
    f = mysubs(f,system.state.variable,mu);
    f = mysubs(f,system.parameter.variable,c);
    f = mysubs(f,system.kappa.variable,c_kappa);
end

%% ASSEMBLE OUTPUT

MM.sym.state.order = I;
MM.sym.state.moments = M;
MM.sym.state.ho_moments = hoM_used;
MM.sym.state.ho_moments_closure = hoM_closure;
MM.sym.state.dMdt_raw = dMdt_raw;
MM.sym.state.derivative = f;
MM.sym.state.M0 = M0;
if isfield(system,'output')
    MM.sym.output.function = H;
else
    MM.sym.output.order = [];
    MM.sym.output.moments = [];
    MM.sym.output.function = [];
end

end

% better subs
function out = mysubs(in, old, new)
if(~isnumeric(in) && ~isempty(old) && ~isempty(findsym(in)))
    matVer = ver('MATLAB');
    if(str2double(matVer.Version)>=8.1)
        out = subs(in, old(:), new(:));
    else
        out = subs(in, old(:), new(:), 0);
    end
else
    out = in;
end
end