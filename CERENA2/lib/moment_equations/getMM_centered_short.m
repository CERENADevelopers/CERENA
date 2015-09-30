% function MM = getMM_centered(system,options)
function MM = getMM_centered_short(varargin)

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

%% INITIALIZATION (1)
n_s  = length(system.state.variable);
if isfield(system,'output')
    n_sy = length(system.output.variable);
end
n_r = length(system.reaction);
n_c = length(system.parameter.variable);

%% INITIALIZATION (2)
% Generation of I-index
[I ,~] = getMomentIndexSet(n_s ,options.moment_order+2);
hoI = I(find(sum(I~=0,2)> options.moment_order),:);
I   = I(find(sum(I~=0,2)<=options.moment_order),:);
alpha   = convertI2alpha(  I,n_s);
hoalpha = convertI2alpha(hoI,n_s);
mu_ind = I(find(sum(I>=1,2) == 1),:);
C_ind  = I(find((2 <= sum(I>=1,2)).*(sum(I>=1,2) <= options.moment_order)),3:end);
covar_ind  = I(2 == sum(I>=1,2),:);
var_ind  = covar_ind(covar_ind(:,end-1)==covar_ind(:,end),:);
% Mean
mu  = getMu(mu_ind);
% Moments with order >= 2
C  = getC(C_ind);
% All moments
M  = [mu ;C ];
n_I = length(M);
hoM = getC(hoI);
n_hoI = length(hoM);

%outputs
if isfield(system,'output')
    [Iy,n_Iy] = getMomentIndexSet(n_sy,options.moment_order_output);
    % Mean
    muy = getMu(Iy(1:n_sy,:),'y');
    % Moments with order >= 2
    Cy = getC(Iy(n_sy+1:end,:),'y');
    % All moments
    My = [muy;Cy];
end

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

%% REACTION STOICHIOMETRY
S_e = system.eductStoichiometry;
S_p = system.productStoichiometry;
% Overall reaction stoichiometry
S = system.stoichiometry;

% Initial conditions
if options.moment_order == 1
    M0 = [system.state.mu0];
    M0 = mysubs(M0,system.parameter.variable,c);
    M0 = sym(M0);
elseif options.moment_order > 1
    M0 = sym(zeros(n_I,1));
    M0(1:length(system.state.mu0)+length(system.state.C0)) = [system.state.mu0;system.state.C0];
    M0 = mysubs(M0,system.parameter.variable,c);
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
        % Propensities
        h{o} = system.output.function(o);
        h{o} = mysubs(h{o},system.parameter.variable,c);
        exph{o} = mysubs(h{o},system.state.variable,X);
        h{o} = mysubs(h{o},system.state.variable,mu);
        dhdx{o} = simplify(jacobian(h{o},mu));        
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
    % Initialize
    H = sym(zeros(n_Iy,1));
    % Loop: Output moments
    for i = 1:length(H)
        % Assignment of I-index and alpha-index
        Iyi = Iy(i,find(Iy(i,:)~=0)); J = Iyi(end);
        alphaiy = [];
        for k = 1:n_sy
            alphaiy(k) = sum(Iyi == k);
        end
        %% MEAN
        if sum(alphaiy) == 1
            H(i) = dhdx{J}*mu; % H(i) = termy{J}(1);
            %% MOMENTS WITH ORDER >= 2
        else
            expH = [];
            for k=1:n_sy
                expH = [expH; exph{k} - h{k}];
            end
            expH = prod(expH(Iyi));
            expH = expand(expH);
            H(i) = mysubs(expH,XM(end:-1:1),uncentM(end:-1:1));
            
        end
        
    end
    H = simplify(H);
    disp('output map done.')
else
    H = [];
end

%% EXTRACTING IMPORTANT COVARIANCES
M2 = M;
hoM1 = [];
hoM_closure1 = [];
if isfield(options,'covReduc')
    if strcmp(options.covReduc,'true')
        ind_var = find((I(:,3)==I(:,4)).*(sum(I~=0,2)==2));
        ind_mean = find(sum(I>=1,2) == 1);
        indEq = [ind_mean;ind_var];
        % Extract the moments that show up in the equations for mean and variance
        MeanVarEq = dMdt(indEq);
        symVar = symvar(MeanVarEq);
        covar = transpose(setdiff(symVar,[c;mu]));
        covar = setdiff(covar,hoM);
        % exlude higher
        if isfield(system,'input')
            covar = setdiff(covar,system.input.variable);
        end
        if isfield(system,'addCovar')
            covar = [system.addCovar;covar];
        end
        M_reduced = unique([M(indEq);covar]);
        
        for i=1:length(M_reduced)
            ind_M(i) = find(ismember(M,M_reduced(i)));
        end
        ind_M = sort(ind_M);
        ind_hoM = sort(setdiff(1:n_I,ind_M));
        MAll = M;
        M0All = M0;
        IAll = I;
        dMdtAll = dMdt;
        hoMAll = hoM;
        hoIAll = hoI;
        M = M(ind_M);
        M0 = M0(ind_M);
        dMdt = dMdt(ind_M);
        I = I(ind_M,:);
        hoM = [MAll(ind_hoM);hoM];
        hoI = [IAll(ind_hoM,:);hoI];
        if strcmp(options.covRed_closure,'special')
            hoM1 = MAll(ind_hoM);
            hoI1 = IAll(ind_hoM,:);
            hoM_closure1 = covRedMomClosure(hoI1,system);
            while ~isempty(setdiff(symvar(hoM_closure1),M))
                hoM_closure1 = mysubs(hoM_closure1,hoM1,hoM_closure1);
            end
            M2 = MAll;
        elseif strcmp(options.covRed_closure,'normal')
            M2 = M;
        end
    end
end
%% MOMENT CLOSURE
% Determine moments of order > options.moment_order
hoM_used2 = transpose(setdiff(symvar(dMdt),[c;M2;system.time]));
if isfield(system,'input')
    hoM_used2 = setdiff(symvar(hoM_used2),system.input.variable);
end
hoM_closure2 = sym([]);
if ~isempty(hoM_used2)
    disp('Moment-closure used!')
    if strcmp(options.moment_closure,'low-dispersion')
        hoM_closure2 = sym(zeros(size(hoM_used2)));
        
    elseif strcmp(options.moment_closure,'derivative-matching')
        I_alpha_C = convertI2alpha(C_ind,n_s);
        I_alpha_mu = convertI2alpha(mu_ind,n_s);
        I_alpha = [I_alpha_mu;I_alpha_C];
        uncentMom = convertUncent2Cent(C_ind,I_alpha_C,mu,n_s);
        uncentMom = [mu;uncentMom];
    elseif strcmp(options.moment_closure,'log-normal')
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
    end
    
    % Closure
    
    for i = 1:length(hoM_used2)
        j = find(hoM==hoM_used2(i));
        Ii = hoI(j,(find(hoI(j,:)~=0)));
        switch options.moment_closure
            case 'low-dispersion'
                
            case 'mean-field'
                uI = unique(Ii);
                hoM_parts = sym(zeros(length(uI),1));
                for j = 1:length(uI)
                    Iij = Ii(Ii==uI(j));
                    hoM_parts(j) = getC(Iij);
                end
                hoM_closure2(i) = prod(hoM_parts);
                %             case 'low-dispersion'
                %                 hoM_closure = sym(zeros(size(hoM_used)));
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
                hoM_closure2(i) = solve(K,hoM_used2(i));
            case 'derivative-matching'
                Ii = hoI(j,:);
                alpha_bar = convertI2alpha(Ii,n_s);
                gamma_p = getDerMatchExponent(alpha_bar,I_alpha);
                hoM_uncent = convertUncent2Cent(Ii,alpha_bar,mu,n_s);
                eq = hoM_uncent - prod(uncentMom.^gamma_p);
                hoM_closure2(i) = solve(eq,hoM_used2(i));
            case 'log-normal'
                Ii = hoI(j,:);
                alpha_bar = convertI2alpha(Ii,n_s);
                hoM_uncent = convertUncent2Cent(Ii,alpha_bar,mu,n_s);
                lognormUncentMom = exp(alpha_bar*munorm + 1/2*alpha_bar*covarnormMat*alpha_bar');
                eq = hoM_uncent - lognormUncentMom;
                hoM_closure2(i) = solve(eq,hoM_used2(i));
            otherwise
                error('This option is not available.');
        end
    end
end

hoM_used = [transpose(hoM1),hoM_used2];
hoM_closure = [transpose(hoM_closure1),hoM_closure2];
if isfield(system,'input')
    while (~isempty(setdiff(symvar(hoM_closure),[c;M;system.time;system.input.variable])))
        hoM_closure = mysubs(hoM_closure,hoM_used,hoM_closure);
    end
else
    while (~isempty(setdiff(symvar(hoM_closure),[c;M;system.time])))
        hoM_closure = mysubs(hoM_closure,hoM_used,hoM_closure);
    end
end
% Substitutions
if isfield(system,'input')
    while (~isempty(setdiff(symvar(dMdt),[c;M;system.time;system.input.variable])))
        dMdt = mysubs(dMdt,hoM_used,hoM_closure);
    end
else
    while (~isempty(setdiff(symvar(dMdt),[c;M;system.time])))
        dMdt = mysubs(dMdt,hoM_used,hoM_closure);
    end
end
if isfield(system,'input')
    while (~isempty(setdiff(symvar(H),[c;M;system.time;system.input.variable],'legacy')))
        H = mysubs(H,hoM_used,hoM_closure);
    end
else
    while (~isempty(setdiff(symvar(H),[c;M;system.time],'legacy')))
        H = mysubs(H,hoM_used,hoM_closure);
    end
end

% dMdt = simplify(dMdt);
disp('moment closure done.')
%% EVALUATION OF JACOBIAN
f = dMdt;
if isfield(system,'input')
    f = mysubs(f,system.input.variable,system.input.function);
    f = mysubs(f,system.state.variable,mu);
    f = mysubs(f,system.parameter.variable,c);
end

if isfield(options,'scale')
    
    if strcmp(options.scale,'concentration')
        scale = system.volumes;
        scale = mysubs(scale,system.parameter.variable,c);
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
        f = simplify(f);

        % rescale initial conditions
        M0 = M0./scale_M;
        M0 = simplify(M0);
        
    elseif strcmp(options.scale,'absolute')
        f = simplify(f);
    else
        error('The specified scaling is not available!');
    end
else
    f=simplify(f);
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
end
if isfield(options,'covReduc')
    if strcmp(options.covReduc,'true')
        MM.sym.state.reducedInd = ind_M;
    end
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