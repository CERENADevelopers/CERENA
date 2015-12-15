function cMFSP = getcMM_centered(varargin)

%% CHECK / ASSIGN INPUTS
% Check required inputs
if nargin >= 3
    if ~isempty(varargin{1})
        % Assign inputs
        system = varargin{1};
        xmin = varargin{2};
        xmax = varargin{3};
    else
        % Error message:
        error('The first input must be non-empty!');
    end
else
    % Error message:
    error('At least three input is required!');
end

% Set default bounds
if isempty(xmin)
    xmin = floor(system.state.xmin);
end
if isempty(xmax)
    xmax = floor(system.state.xmax);
end

% Check for additional constraints g(x) == 1
g = [];
if nargin >= 4
    g = varargin{4};
end
if isempty(g)
    g = @(x) 1;
end

% Options
options.moment_order = 2;
options.moment_closure = 'low-dispersion';
options.filename = 'MCM_res';
options.filepath = '';
options.sensitivity = 'off';
if nargin >= 5
    options = setdefault(varargin{5},options);
end
if ~isfield(options,'moment_order_max') || isempty(options.moment_order_max)
    options.moment_order_max = options.moment_order + 2;
end

%% INITIALIZATION
n_s = length(system.state.variable);
n_r = length(system.reaction);
n_c = length(system.parameter.variable);
n_kappa = length(system.kappa.variable);

%% INITIALIZE ROUTINE
% Boundaries
xmin = max(max([floor(xmin),floor(system.state.xmin)],[],2),0);
xmax =     min([ ceil(xmax), ceil(system.state.xmax)],[],2);

% Idenification of species
y_ind = find(strcmp(system.state.type,'stochastic')); n_y = length(y_ind);
z_ind = find(strcmp(system.state.type,'moment'));     n_z = length(z_ind);
% Reindexingy
system.state.variable = system.state.variable([y_ind;z_ind]);
system.state.type     = system.state.type([y_ind;z_ind]);
system.state.name     = system.state.name([y_ind;z_ind]);
system.state.xmin     = system.state.xmin([y_ind;z_ind]);
system.state.xmax     = system.state.xmax([y_ind;z_ind]);
system.state.mu0      = system.state.mu0([y_ind;z_ind]);
% Extract initial conditions for the moments of 'moment'-species:
[ind_M0,~] = getMomentIndexSet(n_s,2);
ind_C0  = ind_M0(2 == sum(ind_M0>=1,2),:);
n_C0 = size(ind_C0,1);
ind_C0_z = [];
for i = 1:n_C0
    if all(ismember(ind_C0(i,:),z_ind))
        ind_C0_z = [ind_C0_z;i];
    end
end
system.state.C0 = system.state.C0(ind_C0_z);

% Trunctation
ymin = xmin(y_ind);
ymax = xmax(y_ind);
max_index = prod(ymax-ymin+1);
[~,ind] = sort([y_ind;z_ind]);
g = @(X) g(X(ind)) && system.state.constraint(X(ind));

%% GENERATE STATE INDEX SET:
% Note: The states of the FSP are ordered according to the index set.
% Initialize index matrix:
y_index = zeros(max_index,n_y);
% Initialize state
ind = ymin;
% Construct index matrix
for i = 1:prod(ymax-ymin+1)
    % Assign index
    y_index(i,:) = ind';
    % Uypdate current index:
    j = min(setdiff(1:n_y,find(ind == ymax))); % index which is updated
    ind(j) = ind(j) + 1;
    ind(1:j-1) = ymin(1:j-1); % reset of all previous indexes
end

% Truncate index such that only the indexes which fulfill g(x) <= 0
% are contained:
ind = [];
% Construct index matrix
for i = 1:max_index
    % Evaluate g
    I = [y_index(i,:)';zeros(n_z,1)];
    if g(I)
        ind = [ind;i];
    end
end
% Assign indexes which fulfill constraint
y_index = y_index(ind,:);

% Number of states of FSP:
n_y_index = size(y_index,1);

%% CONSTRUCTION OF STATE VECTOR AND PARAMETER VECTOR
% Conditional expectation
[M_ind,n_M] = getMomentIndexSet(n_z,options.moment_order_max);
% moments described by model
mu_ind = M_ind(find(sum(M_ind>=1,2) == 1),:);
C_ind  = M_ind(find((2 <= sum(M_ind>=1,2)).*(sum(M_ind>=1,2) <= options.moment_order)),options.moment_order_max-options.moment_order+1:end);
n_C = size(C_ind,1);
% higher order moments
hoC_ind = M_ind(find(sum(M_ind>=1,2) >  options.moment_order),:);
n_hoC = size(hoC_ind,1);

% State vector
M = sym(zeros((1+n_z+n_C)*n_y_index,1));
% Loop: discrete states
for iy = 1:n_y_index
    % Name of discrete state
    y_index_str{iy} = ['_y' strrep(num2str(y_index(iy,:),'_%d'),' ','')];
    
    % Marginal probability
    p(iy) = sym(['p' y_index_str{iy}]);
    % Conditional mean
    mu{iy}  = getMu(mu_ind,y_index_str{iy});
    % Conditional central moments of order >= 2
    C{iy}   = getC(C_ind  ,y_index_str{iy});
    hoC{iy} = getC(hoC_ind,y_index_str{iy});
    
    % Vector of states
    M(iy) = p(iy);
    M((n_z+n_C)*(iy-1)+n_y_index    +[1:n_z]) = mu{iy};
    M((n_z+n_C)*(iy-1)+n_y_index+n_z+[1:n_C]) = C{iy};
    
    hoM((iy-1)*n_hoC+[1:n_hoC]) = hoC{iy};
end

% Derivative vector
dMdt = sym(zeros((1+n_z+n_C)*n_y_index,1));
for i = 1:length(M)
    dMdt(i) = sym(['d' char(M(i)) 'dt']);
end

%outputs
if isfield(system,'output')
    n_so = length(system.output.variable);
    [Io,n_Io] = getMomentIndexSet(n_so,options.moment_order_output);
    % Mean
    muo = getMu(Io(1:n_so,:),'o');
    % Moments with order >= 2
    Co = getC(Io(n_so+1:end,:),'o');
    % All moments
    Mo = [muo;Co];
end

% Parameters
c = sym(zeros(n_c,1));
for i = 1:n_c
    c(i) = sym(['c' num2str(i,'%d')]);
end
% Constant parameters (kappa)
c_kappa = sym(zeros(n_kappa,1));
for i = 1:n_kappa
    c_kappa(i) = sym(['c_kappa' num2str(i,'%d')]);
end
%% REACTION STOICHIOMETRY
S_e = zeros(n_s,n_r);
S_p = zeros(n_s,n_r);
% Loop: reatcions
for j = 1:n_r
    % Loop: educts
    for i = 1:length(system.reaction(j).educt)
        ind = find(system.state.variable == system.reaction(j).educt(i));
        S_e(ind,j) = S_e(ind,j) + 1;
    end
    % Loop: products
    for i = 1:length(system.reaction(j).product)
        ind = find(system.state.variable == system.reaction(j).product(i));
        S_p(ind,j) = S_p(ind,j) + 1;
    end
end
% Overall reaction stoichiometry
S = S_p - S_e;

% Check reaction order
if max(sum(S_e,1)) >= 3
    error('At most bimolecular reactions are allowed.');
end

% Individual reaction stoichiometry
S_y  = S(1:n_y,:);
if isempty(S_y)
    S_y = 0;
end
S_ey = S_e(1:n_y,:);
if isempty(S_ey)
    S_ey = 0;
end
S_z  = S(n_y+1:end,:);
if isempty(S_z)
    S_z = 0;
end
S_ez = S_e(n_y+1:end,:);
if isempty(S_ez)
    S_ez = 0;
end

%% Initial conditions
M0 = sym(1e-10+zeros((1+n_z+n_C)*n_y_index,1));
y_index_0 = transpose(system.state.fmu0(1:n_y));
ind_y_0 = find(ismember(y_index , y_index_0,'rows'));
M0(ind_y_0) = 1-(n_y-1)*1e-10;
for iy = 1:n_y_index
    if options.moment_order >= 2
        M0(n_y_index + (iy-1)*(n_z+n_C)+[1:n_z+length(system.state.C0)]) = [system.state.mu0(n_y+1:end);system.state.C0];
    else
        M0(n_y_index + (iy-1)*n_z+[1:n_z]) = system.state.mu0(n_y+1:end);
    end
end
M0 = mysubs(M0,system.parameter.variable,c);
M0 = mysubs(M0,system.kappa.variable,c_kappa);
M0 = sym(M0);
%% ASSIGN PROPENSITIES AND COMPUTE DERIVATIVES
% Moment states
z = sym(zeros(n_z,1));
for i = 1:n_z
    z(i) = sym(['z' num2str(i,'%d')]);
end
% Loop: reactions
for j = 1:n_r
    % c
    rate(j) = mysubs(system.reaction(j).parameter,system.parameter.variable,c);
    rate(j) = mysubs(rate(j),system.kappa.variable,c_kappa);
    term{j} = sym([]);
    L{j} = [];
    % h_y(y)
    switch max(S_ey(:,j))
        case 0
            hy{j} = @(y) 1;
        case 1
            ind = find(S_ey(:,j)==1);
            if length(ind) == 1
                hy{j} = @(y) y(ind);
            else
                hy{j} = @(y) y(ind(1))*y(ind(2));
            end
        case 2
            ind = find(S_ey(:,j)==2);
            hy{j} = @(y) y(ind)*(y(ind)-1)/2;
    end
    
    % h_z(z)
    switch max(S_ez(:,j))
        case 0
            hz{j} = 1;
        case 1
            ind = find(S_ez(:,j)==1);
            if length(ind) == 1
                hz{j} = z(ind);
            else
                hz{j} = z(ind(1))*z(ind(2));
            end
        case 2
            ind = find(S_ez(:,j)==2);
            hz{j} = z(ind)*(z(ind)-1)/2;
    end
    hz{j} = sym(hz{j});
    dhzdz{j} = simplify(jacobian(hz{j},z));
    d2hzdz2{j} = simplify(hessian(hz{j},z));
    
    %
    term{j}(1) = hz{j};
    L{j}(1,1:2) = 0;
    for k = 1:n_z
        term{j}(end+1) = dhzdz{j}(k);
        L{j}(end+1,1:2) = [0,k];
    end
    for k = 1:n_z
        for l = k:n_z
            if k == l
                term{j}(end+1) = 0.5*d2hzdz2{j}(k,l);
            else
                term{j}(end+1) =     d2hzdz2{j}(k,l);
            end
            L{j}(end+1,1:2) = [k,l];
        end
    end
    
    % Reduction to non-zero terms
    ind = [];
    for k = 1:length(term{j})
        if ~strcmp(char(term{j}(k)),'0')
            ind = [ind,k];
        end
    end
    term{j} = term{j}(ind);
    L{j} = L{j}(ind,:);
    
    % Function to evaluate expectation
    E_CI_hz{j} = @(I,str) mysubs(term{j},z,getMu(mu_ind,str)) * getC(sort([L{j},repmat(I,size(L{j},1),1)],2),str);
end

%% CONSTRUCTION OF DERIVATIVES
% Stochastic states
dpdt   = sym(zeros(    n_y_index,1));
pdmudt = sym(zeros(n_z*n_y_index,1));
pdCdt  = sym(zeros(n_C*n_y_index,1));
% Loop: index
for iy = 1:n_y_index
    
    %% MARGINAL PROBABILITIES
    % Loop: reactions
    for j = 1:n_r
        % Sum 1
        y_  = y_index(iy,:)'-S_y(:,j);
        iy_ = getRowInMat(y_index,y_);
        % Check if initial state is in reachable set of CME
        if ~isempty(iy_)
            sy_ = y_index_str{iy_};
            dpdt(iy) = dpdt(iy) + rate(j)*hy{j}(y_)*E_CI_hz{j}([],sy_)*p(iy_);
        end
        
        % Sum 2
        y  = y_index(iy,:)';
        sy = y_index_str{iy};
        dpdt(iy) = dpdt(iy) - rate(j)*hy{j}(y)*E_CI_hz{j}([],sy)*p(iy);
    end
    
    %% CONDITIONAL MEANS
    % Loop: reactions
    for j = 1:n_r
        % Sum 1
        y_  = y_index(iy,:)'-S_y(:,j);
        iy_ = getRowInMat(y_index,y_);
        % Check if initial state is in reachable set of CME
        if ~isempty(iy_)
            sy_ = y_index_str{iy_};
            % Loop: conditional means
            for iz = 1:n_z
                pdmudt(n_z*(iy-1)+iz) = pdmudt(n_z*(iy-1)+iz) ...
                    + rate(j)*hy{j}(y_)*(E_CI_hz{j}(iz,sy_)+(getMu(iz,sy_)+S_z(iz,j))*E_CI_hz{j}([],sy_))*p(iy_);
            end
        end
        
        % Sum 2
        y  = y_index(iy,:)';
        sy = y_index_str{iy};
        % Loop: conditional means
        for iz = 1:n_z
            pdmudt(n_z*(iy-1)+iz) = pdmudt(n_z*(iy-1)+iz) ...
                - rate(j)*hy{j}(y)*(E_CI_hz{j}(iz,sy)+getMu(iz,sy)*E_CI_hz{j}([],sy))*p(iy);
        end
    end
    
    % Component 3
    for iz = 1:n_z
        pdmudt(n_z*(iy-1)+iz) = pdmudt(n_z*(iy-1)+iz) ...
            - getMu(iz,sy)*dpdt(iy);
    end
    
    %% CONDITIONAL CENTRAL MOMENTS
    % Loop: reactions
    for j = 1:n_r
        % Sum 1
        y_  = y_index(iy,:)'-S_y(:,j);
        iy_ = getRowInMat(y_index,y_);
        % Check if initial state is in reachable set of CME
        if ~isempty(iy_)
            sy_ = y_index_str{iy_};
            % Loop: conditional central moments
            for iC = 1:n_C
                % Assignment of I-index and alpha-index
                I_I = C_ind(iC,find(C_ind(iC,:)~=0));
                I_alpha = getAlphaIndex(I_I,n_z);
                % Generation of combinations
                [K_alpha,K_I] = getLset(I_I,n_z,'full');
                % Loop: combination
                for k = 1:size(K_alpha,1)
                    k_I = K_I(k,find(K_I(k,:)~=0));
                    k_alpha = K_alpha(k,:);
                    % Prefactor
                    s = nchoosek_vec(I_alpha,k_alpha)*prod((mu{iy_}-mu{iy}+S_z(:,j)).^(I_alpha(:)-k_alpha(:)));
                    pdCdt(n_C*(iy-1)+iC) = pdCdt(n_C*(iy-1)+iC) ...
                        + rate(j)*s*hy{j}(y_)*E_CI_hz{j}(k_I,sy_)*p(iy_);
                end
            end
        end
        
        % Sum 2
        y  = y_index(iy,:)';
        sy = y_index_str{iy};
        % Loop: conditional central moments
        for iC = 1:n_C
            pdCdt(n_C*(iy-1)+iC) = pdCdt(n_C*(iy-1)+iC) ...
                - rate(j)*hy{j}(y)*E_CI_hz{j}(C_ind(iC,:),sy)*p(iy);
        end
    end
    
    % Sum 3
    % Loop: conditional central moments
    for iC = 1:n_C
        % Assignment of I-index and alpha-index
        I_I = C_ind(iC,:);
        I_alpha = getAlphaIndex(I_I,n_z);
        % Set of variable which are contained
        Iz = find(I_alpha >= 1);
        % Loop: high-abundance species with I_i >= 1
        for iz = Iz(:)'
            iz_ = find(I_I==iz,1);
            I_I_eiz = I_I([1:iz_-1,iz_+1:end]);
            pdCdt(n_C*(iy-1)+iC) = pdCdt(n_C*(iy-1)+iC) ...
                - I_alpha(iz)*getC(I_I_eiz,sy)*pdmudt(n_z*(iy-1)+iz)*p(iy);
        end
    end
    
    % Sum 4
    % Loop: conditional central moments
    for iC = 1:n_C
        pdCdt(n_C*(iy-1)+iC) = pdCdt(n_C*(iy-1)+iC) ...
            - getC(C_ind(iC,:),sy)*dpdt(iy);
    end
    
end

% Simplification of right-hand side
dpdt = simplify(dpdt);
pdmudt = simplify(pdmudt);
pdCdt = simplify(pdCdt);


%% CONSTRUCT OUTPUT MAP
% Overall moments of species
[H1,H_ind1] = getMomentsBarMu_sym(y_index,y_ind,z_ind,C_ind,p,mu,C,options.moment_order);
H1 = simplify(H1);
if isfield(system,'output')
    % Moments of output
    [H2,My,Iy,n_Iy,M_H2,I_H2] = getOutputMoments(system,H_ind1,[],[],[],[],c,c_kappa,options);
    H_ind2 = Iy;
    H2 = mysubs(H2,M_H2,H1);
else
    H_ind2 = [];
    H2 = [];
end
H = [H1,transpose(H2)];
disp('output map done.')

%% MOMENT CLOSURE
% Higher order moments contained in equations
hoC_used = transpose(setdiff(symvar([dpdt;pdmudt;pdCdt]),[c;c_kappa;M;system.time]));
if isfield(system,'input')
    hoC_used = setdiff(symvar(hoC_used),system.input.variable);
end
if isfield(system,'output')
    hoC_used = setdiff(symvar(hoC_used),system.output.variable);
end
choC = sym(zeros(length(hoC_used),1));
ind_hoC_used = 1:length(hoC_used);
if ~isempty(hoC_used)
    disp('Moment-closure used!')
    %%%%   USER-DEFINED CLOSURE
    if strcmp(options.moment_closure,'user-defined')
        ind_hoC_used = [];
        for i = 1:length(hoC_used)
            % Properties
            for iy = 1:length(hoC)
                if ~isempty(find(hoC{iy} == hoC_used(i)))
                    indn = find(hoC{iy} == hoC_used(i));
                    I = hoC_ind(indn,(find(hoC_ind(indn,:)~=0)));
                    iy_hoC = iy;
                end
            end
            nargout_ud = nargout('ud_closure');
            if nargout_ud == 1
                choC_user = ud_closure(I,options.moment_order,system);
                symvar_choc = symvar(choC_user);
                for j = 1:length(symvar_choc)
                    choC_user = subs(choC_user, symvar_choc(j),sym([char(symvar_choc(j)),'_',y_index_str{iy_hoC}]));
                end
                choC(i) = choC_user;
                %                 ind_hoC_used = [];
            elseif nargout_ud == 2
                [choC_user,flag_closure(i)] = ud_closure(I,options.moment_order,system);
                if flag_closure(i) ~= 0
                    symvar_choc = symvar(choC_user);
                    for j = 1:length(symvar_choc)
                        choC_user = subs(choC_user, symvar_choc(j),sym([char(symvar_choc(j)),'_',y_index_str{iy_hoC}]));
                    end
                    choC(i) = choC_user;
                else
                    ind_hoC_used = [ind_hoC_used,i];
                end
                options.moment_closure = 'zero-cumulants';
            elseif nargout_ud >= 3
                [choC_user,flag_closure(i),moment_closure_ud] = ud_closure(I,options.moment_order,system);
                if flag_closure(i) ~= 0
                    symvar_choc = symvar(choC_user);
                    for j = 1:length(symvar_choc)
                        choC_user = subs(choC_user, symvar_choc(j),sym([char(symvar_choc(j)),y_index_str{iy_hoC}]));
                    end
                    choC(i) = choC_user;
                else
                    ind_hoC_used = [ind_hoC_used,i];
                end
                options.moment_closure = moment_closure_ud;
            end
        end
        
    end
    % Loop: higher order moments
    %     for i = 1:length(hoC_used)
    for i = 1:ind_hoC_used
        % Properties
        for iy = 1:length(hoC)
            if ~isempty(find(hoC{iy} == hoC_used(i)))
                indn = find(hoC{iy} == hoC_used(i));
                I = hoC_ind(indn,(find(hoC_ind(indn,:)~=0)));
                iy_hoC = iy;
            end
        end
        switch options.moment_closure
            case 'mean-field'
                uI = unique(I);
                hoC_parts = sym(zeros(length(uI),1));
                for j = 1:length(uI)
                    Ii = I(I==uI(j));
                    if length(Ii)==1
                        hoC_parts(j) = getMu(Ii,y_index_str{iy_hoC});
                    elseif length(Ii) > options.moment_order
                        Ii_1 = floor(length(Ii)/options.moment_order);
                        Ii_2 = mod(length(Ii),options.moment_order);
                        if Ii_2 <= 1
                            hoC_parts(j) = (getC(Ii(1:options.moment_order),y_index_str{iy_hoC}))^Ii_1 * getMu(Ii(1),y_index_str{iy_hoC});
                        else
                            hoC_parts(j) = (getC(Ii(1:options.moment_order),y_index_str{iy_hoC}))^Ii_1 * getMu(Ii(1:Ii_2),y_index_str{iy_hoC});
                        end
                    else
                        hoC_parts(j) = getC(Ii,y_index_str{iy_hoC});
                    end
                end
                choC(i) = prod(hoC_parts);
            case 'low-dispersion'
                choC(i) = 0;
            case 'derivative-matching'
                I = hoC_ind(indn,:);
                alpha_bar = convertI2alpha(I,n_z);
                I_alpha_C = convertI2alpha(C_ind,n_z);
                I_alpha_mu = convertI2alpha(mu_ind,n_z);
                I_alpha = [I_alpha_mu;I_alpha_C];
                I_I = C_ind;
                gamma_p = getDerMatchExponent(alpha_bar,I_alpha);
                uncentMom = convertUncent2Cent(I_I,I_alpha_C,mu,n_z,y_index_str,iy_hoC);
                choC_uncent = convertUncent2Cent(I,alpha_bar,mu,n_z,y_index_str,iy_hoC);
                uncentMom = [mu{iy_hoC};uncentMom];
                eq = choC_uncent;
                choC(i) = solve(eq,hoC_used(i));
                choC(i) = choC(i) + prod(uncentMom.^gamma_p);
            case 'zero-cumulants'
                % generate all possible partitions in YI
                B = partitions(I);
                % generating cumulant corresponding to this moment
                K = sym(0);
                for b = 1:length(B)
                    parts = [];
                    for ib = 1:length(B{b})
                        alpha_ib = convertI2alpha(B{b}{ib},n_z);
                        uncentMom_ib = convertUncent2Cent(B{b}{ib},alpha_ib,mu,n_z,y_index_str,iy_hoC);
                        parts = [parts,uncentMom_ib];
                    end
                    npart = length(parts);
                    K = K + factorial(npart-1) * (-1)^(npart-1) * prod(parts);
                end
                choC(i) = solve(K,hoC_used(i));
            otherwise
                error('This option is not available.');
        end
    end
    if isfield(system,'input')
        while (~isempty(setdiff(symvar(choC),[c;c_kappa;M;system.time;system.input.variable])))
            choC = mysubs(choC,hoC_used,choC);
        end
    else
        while (~isempty(setdiff(symvar(choC),[c;c_kappa;M;system.time])))
            choC = mysubs(choC,hoC_used,choC);
        end
    end
    % Substitution
    if isfield(system,'input')
        while (~isempty(setdiff(symvar(dpdt),[c;c_kappa;M;system.time;system.input.variable])))
            dpdt  = mysubs(dpdt,hoC_used,choC);
        end
    else
        while (~isempty(setdiff(symvar(dpdt),[c;c_kappa;M;system.time])))
            dpdt  = mysubs(dpdt,hoC_used,choC);
        end
    end
    if isfield(system,'input')
        while (~isempty(setdiff(symvar(pdmudt),[c;c_kappa;M;system.time;system.input.variable])))
            pdmudt = mysubs(pdmudt,hoC_used,choC);
        end
    else
        while (~isempty(setdiff(symvar(pdmudt),[c;c_kappa;M;system.time])))
            pdmudt = mysubs(pdmudt,hoC_used,choC);
        end
    end
    %     dpdt   = simplify(dpdt);
    %     pdmudt = simplify(pdmudt);
    if n_C >= 1
        if isfield(system,'input')
            while (~isempty(setdiff(symvar(pdCdt),[c;c_kappa;M;system.time;system.input.variable])))
                pdCdt = mysubs(pdCdt,hoC_used,choC);
            end
        else
            while (~isempty(setdiff(symvar(pdCdt),[c;c_kappa;M;system.time])))
                pdCdt = mysubs(pdCdt,hoC_used,choC);
            end
        end
        %         pdCdt  = simplify(pdCdt);
    end
    
    if isfield(system,'input')
        while (~isempty(setdiff(symvar(H),[c;c_kappa;M;system.time;system.input.variable],'legacy')))
            H = mysubs(H,hoM_used,hoM_closure);
        end
    else
        while (~isempty(setdiff(symvar(H),[c;c_kappa;M;system.time],'legacy')))
            H = mysubs(H,hoM_used,hoM_closure);
        end
    end
end
%% Scale
if isfield(options,'scale')
    
    if strcmp(options.scale,'concentration')
        scale = system.volumes;
        scale = mysubs(scale,system.parameter.variable,c);
        scale = mysubs(scale,system.kappa.variable,c_kappa);
        scale_mu = sym(zeros(n_z,1));
        if isfield(system,'compartments')
            for i = 1:length(scale)
                tmp = find(strcmp(system.state.compartment,system.compartments{i}));
                ind_comp{i} = tmp(strcmp(system.state.type(tmp),'moment'));
                ind_comp{i} = ind_comp{i} - n_y;
                scale_mu(ind_comp{i}) = scale(i);
            end
        else
            ind_comp{1} = 1:n_z;
            scale_mu(:) = scale;
        end
        
        scale_C = sym(zeros(n_C,1));
        
        for iC = 1:n_C
            tmp_ind = C_ind(iC,:);
            tmp_ind = tmp_ind(tmp_ind~=0);
            scale_C(iC) = prod(scale_mu(tmp_ind));
        end
        
        for iy = 1:n_y_index
            dpdt   = mysubs(dpdt  ,mu{iy},mu{iy}.*scale_mu);
            pdmudt = mysubs(pdmudt,mu{iy},mu{iy}.*scale_mu);
            pdCdt  = mysubs(pdCdt ,mu{iy},mu{iy}.*scale_mu);
            dpdt   = mysubs(dpdt  ,C{iy},C{iy}.*scale_C);
            pdmudt = mysubs(pdmudt,C{iy},C{iy}.*scale_C);
            pdCdt  = mysubs(pdCdt ,C{iy},C{iy}.*scale_C);
        end
        for iy = 1:n_y_index
            pdmudt((iy-1)*n_z+[1:n_z]) = pdmudt((iy-1)*n_z+[1:n_z])./scale_mu;
            pdCdt((iy-1)*n_C+[1:n_C]) = pdCdt((iy-1)*n_C+[1:n_C])./scale_C;
        end
        
        % rescale initial conditions
        for iy = 1:n_y_index
            M0((iy-1)*(n_z+n_C) + n_y_index + [1:n_z]) = M0((iy-1)*(n_z+n_C) + n_y_index + [1:n_z])./scale_mu;
            M0((iy-1)*(n_z+n_C) + n_y_index + n_z + [1:n_C] ) = M0((iy-1)*(n_z+n_C) + n_y_index + n_z + [1:n_C])./scale_C;
        end
        
        M0 = simplify(M0);
        
    elseif strcmp(options.scale,'absolute')
        %         dpdt = simplify(dpdt);
        %         pdmudt = simplify(pdmudt);
        %         pdCdt = simplify(pdCdt);
    else
        error('The specified scaling is not available!');
    end
else
    %     dpdt = simplify(dpdt);
    %     pdmudt = simplify(pdmudt);
    %     pdCdt = simplify(pdCdt);
end
%% CONSTRUCTION OF VECTOR FIELD, MASS MATRIX, AND RESIDUAL FUNCTION - STATES
VFsym = dpdt;
MMsym = sym(ones(n_y_index,1));
for iy = 1:n_y_index
    VFsym((n_z+n_C)*(iy-1)+n_y_index    +[1:n_z]) = pdmudt(n_z*(iy-1)+[1:n_z]);
    VFsym((n_z+n_C)*(iy-1)+n_y_index+n_z+[1:n_C]) = pdCdt(n_C*(iy-1)+[1:n_C]);
    MMsym((n_z+n_C)*(iy-1)+n_y_index +[1:n_z+n_C]) = p(iy);
end
if isfield(system,'input')
    VFsym = mysubs(VFsym,system.input.variable,system.input.function);
    VFsym = mysubs(VFsym,system.parameter.variable,c);
    VFsym = mysubs(VFsym,system.kappa.variable,c_kappa);
    MMsym = mysubs(MMsym,system.input.variable,system.input.function);
    MMsym = mysubs(MMsym,system.parameter.variable,c);
    MMsym = mysubs(MMsym,system.kappa.variable,c_kappa);
end
VFsym = simplify(VFsym);
VFsym = collect(VFsym);

RESsym = -MMsym.*dMdt + VFsym;

%% ASSEMBLE OUTPUT
cMFSP.type  = 'centered';
cMFSP.order = options.moment_order;
cMFSP.closure = options.moment_closure;
% Stochastic state
cMFSP.state.stochatic.name = system.state.name(1:n_y);
cMFSP.state.stochatic.variable = system.state.variable(1:n_y);
cMFSP.state.stochatic.xmin = system.state.xmin(1:n_y);
cMFSP.state.stochatic.xmax = system.state.xmax(1:n_y);
cMFSP.state.stochatic.state_index = y_ind;
cMFSP.state.stochatic.FSP_index = y_index;
% Expectation state
cMFSP.state.expectation.name = system.state.name(n_y+1:end);
cMFSP.state.expectation.variable = system.state.variable(n_y+1:end);
cMFSP.state.expectation.xmin = system.state.xmin(n_y+1:end);
cMFSP.state.expectation.xmax = system.state.xmax(n_y+1:end);
cMFSP.state.expectation.state_index = z_ind;
cMFSP.state.expectation.C_index = C_ind;
% Symbolic
cMFSP.state.sym.p = p;
cMFSP.state.sym.mu = mu;
cMFSP.state.sym.C = C;
cMFSP.state.sym.hoC = hoC_used;
cMFSP.state.sym.hoC_closure = choC;
cMFSP.state.sym.M = M;
cMFSP.state.sym.dMdt = dMdt;
cMFSP.state.sym.M0 = M0;

% Symbolic
cMFSP.parameters.sym.c = c;
cMFSP.kappa.sym.c_kappa = c_kappa;
cMFSP.system.stoichiometry = S;

% State of cMFSP
cMFSP.derivatives.sym.dpdt = dpdt;
cMFSP.derivatives.sym.pdmudt = pdmudt;
cMFSP.derivatives.sym.pdCdt = pdCdt;
cMFSP.derivatives.sym.VF = VFsym;
cMFSP.derivatives.sym.mass_matrix = MMsym;
cMFSP.state_moments.order = H_ind1;
cMFSP.output.order = H_ind2;
cMFSP.output.function = H;
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
