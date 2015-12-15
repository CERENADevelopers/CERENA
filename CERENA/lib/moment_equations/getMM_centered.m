% function MM = getMM_centered(system,options)
function MM = getMM_centered(varargin)

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
n_s  = length(system.state.variable);
if isfield(system,'output')
    n_sy = length(system.output.variable);
end
n_r = length(system.reaction);
n_c = length(system.parameter.variable);
n_kappa = length(system.kappa.variable);

%% INITIALIZATION (2)
% Generation of I-index
[I ,~] = getMomentIndexSet(n_s ,options.moment_order_max);
hoI = I(find(sum(I~=0,2)> options.moment_order),:);
I   = I(find(sum(I~=0,2)<=options.moment_order),:);
alpha   = convertI2alpha(  I,n_s);
hoalpha = convertI2alpha(hoI,n_s);
mu_ind = I(find(sum(I>=1,2) == 1),:);
C_ind  = I(find((2 <= sum(I>=1,2)).*(sum(I>=1,2) <= options.moment_order)),options.moment_order_max-options.moment_order+1:end);
covar_ind  = I(2 == sum(I>=1,2),:);
% Mean
mu  = getMu(mu_ind);
% Moments with order >= 2
C  = getC(C_ind);
% All moments
M  = [mu ;C ];
n_I = length(M);
hoM = getC(hoI);
n_hoI = length(hoM);

% sym var for the states
X = sym('X',[n_s,1]);
XM   = prod(repmat(transpose(X),  n_I,1).^alpha  ,2);
hoXM = prod(repmat(transpose(X),n_hoI,1).^hoalpha,2);
uncentM = convertUncent2Cent(I,alpha,mu,n_s);
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
S_e = system.eductStoichiometry;
S_p = system.productStoichiometry;
% Overall reaction stoichiometry
S = system.stoichiometry;

% Initial conditions
if options.moment_order == 1
    M0 = [system.state.mu0];
    M0 = mysubs(M0,system.parameter.variable,c);
    M0 = mysubs(M0,system.kappa.variable,c_kappa);
    M0 = sym(M0);
elseif options.moment_order > 1
    M0 = sym(zeros(n_I,1));
    M0(1:length(system.state.mu0)+length(system.state.C0)) = [system.state.mu0;system.state.C0];
    M0 = mysubs(M0,system.parameter.variable,c);
    M0 = mysubs(M0,system.kappa.variable,c_kappa);
    M0 = sym(M0);
end
%% ASSIGN PROPENSITIES AND COMPUTE DERIVATIVES
for r = 1:n_r
    term{r} = sym([]);
    L{r} = [];
    % Propensities
    w{r} = system.reaction(r).propensity;
    w{r} = mysubs(w{r},system.state.variable,mu);
    %     w{r} = mysubs(w{r},system.input.variable,system.input.function);
    w{r} = mysubs(w{r},system.parameter.variable,c);
    w{r} = mysubs(w{r},system.kappa.variable,c_kappa);
    dwdx{r} = simplify(jacobian(w{r},mu));
    d2wdx2{r} = simplify(hessian(w{r},mu));
    %
    if w{r} ~= 0
        term{r}(1) = w{r};
        L{r}(1,1:2) = 0;
        for k = 1:n_s
            term{r}(end+1) = dwdx{r}(k);
            L{r}(end+1,1:2) = [0,k];
        end
        for k = 1:n_s
            for l = k:n_s
                if k == l
                    term{r}(end+1) = 0.5*d2wdx2{r}(k,l);
                else
                    term{r}(end+1) = d2wdx2{r}(k,l);
                end
                L{r}(end+1,1:2) = [k,l];
            end
        end
    end
    ind = isAlways(term{r}==0);
    term{r}(ind) = [];
    L{r}(ind,:) = [];
end

%% ASSIGN OUTPUT MAP AND COMPUTE DERIVATIVES
if isfield(system,'output')
    for o = 1:n_sy
        termy{o} = sym([]);
        Ly{o} = [];
        % Propensities
        h{o} = system.output.function(o);
        h{o} = mysubs(h{o},system.parameter.variable,c);
        h{o} = mysubs(h{o},system.kappa.variable,c_kappa);
        %         exph{o} = mysubs(h{o},system.state.variable,X);
        h{o} = mysubs(h{o},system.state.variable,mu);
        dhdx{o} = simplify(jacobian(h{o},mu));
        %
        if h{o} ~= 0
            termy{o}(1) = h{o};
            Ly{o}(1,1:2) = 0;
            %         for k = 1:n_sy
            for k = 1:n_s
                termy{o}(end+1) = dhdx{o}(k);
                Ly{o}(end+1,1:2) = [0,k];
            end
        end
        ind = isAlways(termy{o}==0);
        termy{o}(ind) = [];
        Ly{o}(ind,:) = [];
    end
    disp('output derivatives done.')
end
%% CONSTRUCTION OF MOMENT EQUATION
dMdt  = sym(zeros(n_I,1));
% Loop: Moments
for i = 1:n_I
    % Assignment of I-index and alpha-index
    Ii = I(i,find(I(i,:)~=0));
    for k = 1:n_s
        alphai(k) = sum(Ii == k);
    end
    %% MEAN
    if sum(alphai) == 1
        % Loop: reactions
        for r = 1:n_r
            dMdt(i) = dMdt(i) + S(Ii,r)*term{r}*getC(L{r});
        end
        %% MOMENTS WITH ORDER >= 2
    else
        for r = 1:n_r
            % Generation of combinations
            [L_alpha,L_I] = getLset(Ii,n_s);
            % Loop: combinations
            for l = 1:size(L_alpha,1)
                l_alpha = L_alpha(l,:);
                l_I = L_I(l,:);
                l_I = l_I(find(l_I~=0));
                % Prefactor
                factor = nchoosek_vec(alphai,l_alpha)*prod(S(:,r).^(alphai(:)-l_alpha(:)));
                if factor ~= 0
                    dMdt(i) = dMdt(i) + factor*term{r}*getC(sort([L{r},repmat(l_I,size(L{r},1),1)],2));
                end
            end
        end
        % Summand 4
        if length(Ii) >= 3
            for q = 1:n_s
                if alphai(q) >= 1
                    Iiq = Ii; Iiq(find(Iiq==q,1)) = [];
                    dMdt(i) = dMdt(i) - alphai(q)*dMdt(q)*getC(Iiq);
                end
            end
        end
    end
end
disp('construction of dMdt done.')
%% CONSTRUCT OUTPUT MAP
if isfield(system,'output')
    options.output.calculate = 0;
    [H,My,Iy,n_Iy] = getOutputMoments(system,I,M,X,XM,uncentM,c,c_kappa,options);
    disp('output map done.')
else
    H = [];
end

%% MOMENT CLOSURE
% Determine moments of order > options.moment_order
var_system = [c;c_kappa;M;system.time];
if isfield(system,'input')
    var_input = [system.input.variable;transpose(symvar(system.input.function))];
    var_system = [var_system;var_input];
end
if isfield(system,'output')
    var_output = [system.output.variable;transpose(symvar(system.output.function))];
    var_system = [var_system;var_output];
end
hoM_used = transpose(setdiff(symvar(dMdt),var_system));
ind_hoM_used = 1:length(hoM_used);
% hoM_closure2 = sym([]);
hoM_closure = sym(zeros(size(hoM_used)));
if ~isempty(hoM_used)
    disp('Moment-closure used!')
    %%%%   USER-DEFINED CLOSURE
    if strcmp(options.moment_closure,'user-defined')
        Ii = hoI(ismember(hoM,hoM_used),:);
        nargout_ud = nargout('ud_closure');
        if nargout_ud == 1
            hoM_closure = ud_closure(Ii,options.moment_order,system);
            ind_hoM_used = [];
        elseif nargout_ud == 2
            [hoM_closure_user,flag_closure] = ud_closure(Ii,options.moment_order,system);
            hoM_closure(flag_closure~=0)= hoM_closure_user(flag_closure~=0);
            ind_hoM_used = transpose(find(flag_closure==0));
            options.moment_closure = 'zero-cumulants';
        elseif nargout_ud >= 3
            [hoM_closure_user,flag_closure,moment_closure_ud] = ud_closure(Ii,options.moment_order,system);
            hoM_closure(flag_closure~=0)= hoM_closure_user(flag_closure~=0);
            ind_hoM_used = transpose(find(flag_closure==0));
            options.moment_closure = moment_closure_ud;
        end
        % CHECK FOR WHETHER THE SUBSTITUITION STILL
        % INCLUDES ANY HIGHER-ORDER MOMENTS. IS IT NECESSARY,
        % OR SHOULD THE USER BE REQUIRED TO MAKE SURE OF THAT?
        %                     symvar_hoM_closure = symvar(hoM_closure2(i));
        %                         while ~ismember(symvar_hoM_closure(isv),var_system)
        %
        %                         end
    end
    if ~isempty(ind_hoM_used)
        if strcmp(options.moment_closure,'low-dispersion')
            %         hoM_closure2 = sym(zeros(size(hoM_used2)));
        elseif strcmp(options.moment_closure,'log-normal')
            if options.moment_order >= 2
                var_ind  = covar_ind(covar_ind(:,end-1)==covar_ind(:,end),:);
                epsilon = 1e-6;
                munorm = log(epsilon + mu.^2./sqrt(mu.^2+getC(var_ind)+epsilon));
                covarnorm = log(1+getC(covar_ind)./(mu(covar_ind(:,end-1)).*mu(covar_ind(:,end))+epsilon));
                covarnormMat = sym(zeros(n_s));
                k=1;
                for i=1:n_s
                    covarnormMat(i,i:n_s) = covarnorm(k:k+n_s-i);
                    k = k+n_s-i+1;
                end
                covarnormMat = covarnormMat + tril(transpose(covarnormMat),-1);
            else
                warning('Log-Normal closure is not available for 1st-order moment equations! Closure scheme was changed to Derivative-Matching!')
                options.moment_closure = 'derivative-matching';
            end
        end
        if strcmp(options.moment_closure,'derivative-matching')
            I_alpha_C = convertI2alpha(C_ind,n_s);
            I_alpha_mu = convertI2alpha(mu_ind,n_s);
            I_alpha = [I_alpha_mu;I_alpha_C];
            uncentMom = convertUncent2Cent(C_ind,I_alpha_C,mu,n_s);
            uncentMom = [mu;uncentMom];
            
        end
        
        % Closure
        if ~strcmp(options.moment_closure,'low-dispersion')
            %         for i = 1:length(hoM_used2)
            for i = ind_hoM_used
                j = find(hoM==hoM_used(i));
                Ii = hoI(j,(find(hoI(j,:)~=0)));
                switch options.moment_closure
                    %             case 'low-dispersion'
                    %
                    case 'mean-field'
                        uI = unique(Ii);
                        hoM_parts = sym(zeros(length(uI),1));
                        for j = 1:length(uI)
                            Iij = Ii(Ii==uI(j));
                            if length(Iij)==1
                                hoM_parts(j) = getMu(Iij);
                            else
                                hoM_parts(j) = getC(Iij);
                            end
                        end
                        
                        for j = 1:length(uI)
                            nj = sum(Ii==uI(j));
                            hoM_parts(j) = getMu(uI(j))^nj;
                        end
                        
                        hoM_closure(i) = prod(hoM_parts);
                        
                    case 'zero-cumulants'
                        % generate all possible partitions in YI
                        B = partitions(Ii);
                        % generating cumulant corresponding to this moment
                        K = sym(0);
                        for b = 1:length(B)
                            parts = [];
                            for ib = 1:length(B{b})
                                alpha_ib = convertI2alpha(B{b}{ib},n_s);
                                uncentMom_ib = convertUncent2Cent(B{b}{ib},alpha_ib,mu,n_s);
                                parts = [parts,uncentMom_ib];
                            end
                            npart = length(parts);
                            K = K + factorial(npart-1) * (-1)^(npart-1) * prod(parts);
                        end
                        hoM_closure(i) = solve(K,hoM_used(i));
                    case 'derivative-matching'
                        Ii = hoI(j,:);
                        alpha_bar = convertI2alpha(Ii,n_s);
                        gamma_p = getDerMatchExponent(alpha_bar,I_alpha);
                        hoM_uncent = convertUncent2Cent(Ii,alpha_bar,mu,n_s);
                        %                 eq = hoM_uncent - prod(uncentMom.^gamma_p);
                        eq = hoM_uncent;
                        hoM_closure(i) = solve(eq,hoM_used(i));
                        hoM_closure(i) = hoM_closure(i) + prod(uncentMom.^gamma_p);
                    case 'log-normal'
                        Ii = hoI(j,:);
                        alpha_bar = convertI2alpha(Ii,n_s);
                        hoM_uncent = convertUncent2Cent(Ii,alpha_bar,mu,n_s);
                        lognormUncentMom = exp(alpha_bar*munorm + 1/2*alpha_bar*covarnormMat*alpha_bar');
                        eq = hoM_uncent - lognormUncentMom;
                        hoM_closure(i) = solve(eq,hoM_used(i));
                    otherwise
                        error('This option is not available.');
                end
            end
        end
    end
end

while (~isempty(setdiff(symvar(hoM_closure),var_system)))
    hoM_closure = mysubs(hoM_closure,hoM_used,hoM_closure);
end

% Substitutions
while (~isempty(setdiff(symvar(dMdt),var_system)))
    dMdt = mysubs(dMdt,hoM_used,hoM_closure);
end

while (~isempty(setdiff(symvar(H),var_system,'legacy')))
    H = mysubs(H,hoM_used,hoM_closure);
end

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

if isfield(options,'scale')
    
    if strcmp(options.scale,'concentration')
        scale = system.volumes;
        scale = mysubs(scale,system.parameter.variable,c);
        scale = mysubs(scale,system.kappa.variable,c_kappa);
        scale_mu = sym(zeros(n_s,1));
        if isfield(system,'compartments')
            for i = 1:length(scale)
                ind_comp{i} = find(strcmp(system.state.compartment,system.compartments{i}));
                scale_mu(ind_comp{i}) = scale(i);
            end
        else
            ind_comp{1} = 1:system.state.number;
            scale_mu(:) = scale;
        end
        n_C = n_I - n_s;
        scale_C = sym(zeros(n_C,1));
        for iC = 1:n_C
            tmp_ind = C_ind(iC,:);
            tmp_ind = tmp_ind(tmp_ind~=0);
            scale_C(iC) = prod(scale_mu(tmp_ind));
        end
        scale_M = [scale_mu;scale_C];
        
        f = mysubs(f,M,M.*scale_M);
        f = f./scale_M;
        M0 = M0./scale_M;
        M0 = simplify(M0);
        
    elseif strcmp(options.scale,'absolute')
        %         f = simplify(f);
    else
        error('The specified scaling is not available!');
    end
else
    %     f=simplify(f);
end


%% ASSEMBLE OUTPUT
MM.type  = 'centered';
MM.order = options.moment_order;
MM.closure = options.moment_closure;
MM.system.stoichiometry = S;
MM.sym.state.order = I;
MM.sym.state.moments = M;
MM.sym.state.ho_moments = hoM_used;
MM.sym.state.ho_moments_closure = hoM_closure;
% MM.sym.state.derivative = dMdt;
MM.sym.state.derivative = f;
MM.sym.state.M0 = M0;
if isfield(system,'output')
    MM.sym.output.order = Iy;
    MM.sym.output.moments = My;
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