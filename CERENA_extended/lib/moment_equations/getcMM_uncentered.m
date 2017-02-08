function cMFSP = getcMM_uncentered(varargin)

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
options.moment_closure = 'mean field';
options.filename = 'MCM_res';
options.filepath = '';
options.sensitivity = 'off';
if nargin >= 5
    options = setdefault(varargin{5},options);
end
if ~isfield(options,'taylor_order')
    options.taylor_order = [];
end
%% INITIALIZATION
n_s = length(system.state.variable);
n_r = length(system.reaction);
n_c = length(system.parameter.variable);

%% INITIALIZE ROUTINE
% Boundaries
xmin = max(max([floor(xmin),floor(system.state.xmin)],[],2),0);
xmax =     min([ ceil(xmax), ceil(system.state.xmax)],[],2);

% Idenification of species
z_ind = find(strcmp(system.state.type,'stochastic')); n_z = length(z_ind);
y_ind = find(strcmp(system.state.type,'moment'));     n_y = length(y_ind);
% Reindexing
system.state.variable = system.state.variable([z_ind;y_ind]);
system.state.type     = system.state.type([z_ind;y_ind]);
system.state.name     = system.state.name([z_ind;y_ind]);
system.state.xmin     = system.state.xmin([z_ind;y_ind]);
system.state.xmax     = system.state.xmax([z_ind;y_ind]);

% Trunctation
xmin = xmin(z_ind);
xmax = xmax(z_ind);
max_index = prod(xmax-xmin+1);
[~,ind] = sort([z_ind;y_ind]);
g = @(X) g(X(ind)) && system.state.constraint(X(ind));

%% GENERATE STATE INDEX SET:
% Note: The states of the FSP are ordered according to the index set.

% Initialize index matrix:
z_index = zeros(max_index,n_z);
% Initialize state
ind = xmin;
% Construct index matrix
for i = 1:prod(xmax-xmin+1)
    % Assign index
    z_index(i,:) = ind';
    % Update current index:
    j = min(setdiff(1:n_z,find(ind == xmax))); % index which is updated
    ind(j) = ind(j) + 1;
    ind(1:j-1) = xmin(1:j-1); % reset of all previous indexes
end

% Truncate index such that only the indexes which fulfill g(x) <= 0
% are contained:
ind = [];
% Construct index matri~x
for i = 1:max_index
    % Evaluate g
    I = [z_index(i,:)';zeros(n_y,1)];
    if g(I)
        ind = [ind;i];
    end
end
% Assign indexes which fulfill constraint
z_index = z_index(ind,:);

% Number of states of FSP:
n_z_index = size(z_index,1);

%% STATE VECTORS
% Stochastic states
Z = sym(zeros(n_z,1));
for i = 1:n_z
    Z(i) = sym(['Z' num2str(i,'%d')]);
end

% States for which the momements are computed
Y = sym(zeros(n_y,1));
for i = 1:n_y
    Y(i) = sym(['Y' num2str(i,'%d')]);
end

% State vector
X = [Z;Y];

%% PARAMETER VECTOR
c = sym(zeros(n_c,1));
% Loop: parameters
for i = 1:n_c
    c(i) = sym(['c' num2str(i,'%d')]);
end

%% REACTION STOICHIOMETRY
S_e = zeros(n_s,n_r);
S_p = zeros(n_s,n_r);
% Loop: reatcions
for r = 1:n_r
    % Loop: educts
    for j = 1:length(system.reaction(r).educt)
        ind = find(system.state.variable == system.reaction(r).educt(j));
        S_e(ind,r) = S_e(ind,r) + 1;
    end
    % Loop: products
    for j = 1:length(system.reaction(r).product)
        ind = find(system.state.variable == system.reaction(r).product(j));
        S_p(ind,r) = S_p(ind,r) + 1;
    end
end
% Overall reaction stoichiometry
S = S_p - S_e;

% Individual reaction stoichiometry
S_z  = S(1:n_z,:);
S_y  = S(n_z+1:end,:);
%% PROPENSITIES
w = sym(zeros(n_r,1));
% ro_max = 0;
% Loop: reactions
for r = 1:n_r
    w(r) = subs(subs(system.reaction(r).propensity,...
        system.state.variable,[Z;Y],0),...
        system.parameter.variable,c,0);
    %     ro_max = max(ro_max,sum(getMonomeDegree(w(r),Y)));
end
if options.moment_order <4
    ro_max = 2;
else
    ro_max = options.moment_order - 1;
end
%% CONSTRUCTION OF STATE VECTOR p(z|t) AND E(x^I|z,t)
% Stochastic states
p = sym(zeros(n_z_index,1));
% Loop: index
for i = 1:n_z_index
    istr = strrep(num2str(z_index(i,:),'%d_'),' ','');
    p(i) = sym(['p_z_' istr(1:end-1)]);
end
% Conditional expectation
[E_ind,n_E] = getMomentIndexSet(n_y,options.moment_order+ro_max);
% higer order moments
hoE_ind = E_ind(find(sum(E_ind>=1,2) >  options.moment_order),:);
n_hoE = size(hoE_ind,1);
% moments described by model
E_ind = E_ind(find(sum(E_ind>=1,2) <= options.moment_order),:);
n_E = size(E_ind,1);

E_funY = sym(zeros(n_E,1));
E  = sym(zeros(n_E*n_z_index,1));
% Loop: index
for e = 1:n_E
    E_funY(e) = prod(Y(E_ind(e,find(E_ind(e,:)~=0))));
    for i = 1:n_z_index
        istr = strrep(num2str(z_index(i,:),'%d_'),' ','');
%         E(n_E*(i-1)+e) = sym(['E_' num2str(E_ind(e,find(E_ind(e,:)~=0)),'%d') '_z_' istr(1:end-1)]);
        E(n_E*(i-1)+e) = sym(['E_' strrep(num2str(E_ind(e,find(E_ind(e,:)~=0))),'  ','_') '_z_' istr(1:end-1)]);
    end
end

hoE_funY = sym(zeros(n_hoE,1));
hoE  = sym(zeros(n_hoE*n_z_index,1));
% Loop: index
for e = 1:n_hoE
    hoE_funY(e) = prod(Y(hoE_ind(e,find(hoE_ind(e,:)~=0))));
    for i = 1:n_z_index
        istr = strrep(num2str(z_index(i,:),'%d_'),' ','');
        hoE(n_hoE*(i-1)+e) = sym(['E_' num2str(hoE_ind(e,find(hoE_ind(e,:)~=0)),'%d') '_z_' istr(1:end-1)]);
    end
end

%% CONSTRUCTION OF DERIVATIVES dp(z|t)/dt AND dEp(z|t)/dt
% Stochastic states
dpdt = sym(zeros(n_z_index,1));
dEpdt = sym(zeros(n_E*n_z_index,1));
% Loop: index
for i = 1:n_z_index
    for r = 1:n_r
        if 1 %S(i,r) ~= 0
            % Influx
            Z_ind_ = z_index(i,:)'-S_z(:,r);
            % Check if initial state is in reachable set of CME
            if min(system.state.xmin(1:n_z) <= Z_ind_-S_e(1:n_z,r)) && ...
                    min(system.state.xmax(1:n_z) >= Z_ind_-S_e(1:n_z,r))
                % Check if state is contained in FSP
                j = find((bsxfun(@minus,z_index,Z_ind_') == 0)*ones(n_z,1) == n_z);
                if ~isempty(j)
                    % Construct dpdt
                    if isfield(system.reaction(r),'propensityType')
                        switch system.reaction(r).propensityType
                            case 'polynomial'
                                flux = w(r)*p(j);
                                flux = subs(flux,Z,Z_ind_,0);
                            case 'non-polynomial'
                                wr = w(r);
                                wr = subs(wr,Z,Z_ind_,0);
                                wr = getTaylorExpansion(wr,Y,E(n_E*(j-1)+1:n_E*(j-1)+n_y),options.taylor_order);
                                flux = wr*p(j);
                        end
                    else % interpret as polynomial
                        flux = w(r)*p(j);
                        flux = subs(flux,Z,Z_ind_,0);
                    end
                    flux = expand(sym(flux));
                    flux = mysubs(flux,hoE_funY(end:-1:1),hoE(n_hoE*j:-1:n_hoE*(j-1)+1));
                    flux = mysubs(flux,  E_funY(end:-1:1),  E(  n_E*j:-1:  n_E*(j-1)+1));
                    dpdt(i) = dpdt(i) + flux;
                    % Construct dEpdt
                    for e = 1:n_E
                        if isfield(system.reaction(r),'propensityType')
                            switch system.reaction(r).propensityType
                                case 'polynomial'
                                    flux = expand(subs(E_funY(e),Y,Y+S_y(:,r),0)*w(r)*p(j));
                                    flux = subs(flux,Z,Z_ind_,0);
                                case 'non-polynomial'
                                    flux = expand(subs(E_funY(e),Y,Y+S_y(:,r),0)*wr*p(j));
                            end
                        else %interpret as polynomial
                            flux = expand(subs(E_funY(e),Y,Y+S_y(:,r),0)*w(r)*p(j));
                            flux = subs(flux,Z,Z_ind_,0);
                        end
                        flux = expand(sym(flux));
                        flux = mysubs(flux,hoE_funY(end:-1:1),hoE(n_hoE*j:-1:n_hoE*(j-1)+1));
                        flux = mysubs(flux,  E_funY(end:-1:1),  E(  n_E*j:-1:  n_E*(j-1)+1));
                        dEpdt(n_E*(i-1)+e) = dEpdt(n_E*(i-1)+e) + flux;
                    end
                end
            end
            
            % Outflux
            Z_ind_ = z_index(i,:)';
            % Check if initial state is in reachable set of CME
            if min(system.state.xmin(1:n_z) <= Z_ind_-S_e(1:n_z,r)) && ...
                    min(system.state.xmax(1:n_z) >= Z_ind_-S_e(1:n_z,r))
                % Construct dpdt
                if isfield(system.reaction(r),'propensityType')
                    switch system.reaction(r).propensityType
                        case 'polynomial'
                            flux = w(r)*p(i);
                            flux = subs(flux,Z,Z_ind_,0);
                        case 'non-polynomial'
                            wr = w(r);
                            wr = subs(wr,Z,Z_ind_,0);
                            wr = getTaylorExpansion(wr,Y,E(n_E*(i-1)+1:n_E*(i-1)+n_y),options.taylor_order);
                            flux = wr*p(i);
                    end
                else %interpret as polynomial
                    flux = w(r)*p(i);
                    flux = subs(flux,Z,Z_ind_,0);
                end
                flux = expand(sym(flux));
                flux = mysubs(flux,hoE_funY(end:-1:1),hoE(n_hoE*i:-1:n_hoE*(i-1)+1));
                flux = mysubs(flux,  E_funY(end:-1:1),  E(  n_E*i:-1:  n_E*(i-1)+1));
                dpdt(i) = dpdt(i) - flux;
                % Construct dEpdt
                for e = 1:n_E
                    if isfield(system.reaction(r),'propensityType')
                        switch system.reaction(r).propensityType
                            case 'polynomial'
                                flux = expand(E_funY(e)*w(r)*p(i));
                                flux = subs(flux,Z,Z_ind_,0);
                            case 'non-polynomial'
                                flux = expand(E_funY(e)*wr*p(i));
                        end
                    else %interpret as polynomial
                        flux = expand(E_funY(e)*w(r)*p(i));
                        flux = subs(flux,Z,Z_ind_,0);
                    end
                    flux = expand(sym(flux));
                    flux = mysubs(flux,hoE_funY(end:-1:1),hoE(n_hoE*i:-1:n_hoE*(i-1)+1));
                    flux = mysubs(flux,  E_funY(end:-1:1),  E(  n_E*i:-1:  n_E*(i-1)+1));
                    dEpdt(n_E*(i-1)+e) = dEpdt(n_E*(i-1)+e) - flux;
                end
            end
        end
    end
end
dpdt = simplify(dpdt);
dEpdt = simplify(dEpdt);

%% MOMENT CLOSURE
% Higher order moments contained in equations
hoE_used = transpose(setdiff(symvar([dpdt;dEpdt]),[p;c;E;system.time]));
if isfield(system,'input')
    hoE_used = setdiff(symvar(hoE_used),system.input.variable);
end
if ~isempty(hoE_used)
    disp('Moment-closure used!')    
    choE = sym(zeros(length(hoE_used),1));
    % Loop: higher order moments
    for i = 1:length(hoE_used)
        ind = find(hoE == hoE_used(i));
        j = floor((ind-1)/n_hoE) + 1;
        k = ind - n_hoE*(j-1);
        I = hoE_ind(k,(find(hoE_ind(k,:)~=0)));
        switch options.moment_closure
            case 'mean-field'
                %             I = hoE_ind(k,(find(hoE_ind(k,:)~=0)));
                choE(i) = subs(hoE_funY(k),E_funY,E(n_E*(j-1)+1:n_E*j),0);
            case 'low-dispersion'
                %             I = hoE_ind(k,(find(hoE_ind(k,:)~=0)));
                YI = Y(I');
                EI = subs(YI,E_funY(end:-1:1),E(n_E*j:-1:n_E*(j-1)+1),0);
                eq = mysubs(expand(prod(YI - EI)),...
                    [hoE_funY(end:-1:1);E_funY(end:-1:1)],...
                    [hoE(n_hoE*j:-1:n_hoE*(j-1)+1);E(n_E*j:-1:n_E*(j-1)+1)]);
                choE(i) = solve(eq,hoE_used(i));
                choE(i) = simplify(choE(i));
                %             choE(i) = subs(choE(i),hoE_used(1:i-1),choE(1:i-1));
            case 'derivative-matching'
                I = hoE_ind(k,:);
                alpha_bar = convertI2alpha(I,n_y);
                alpha = convertI2alpha(E_ind,n_y);
                
                gamma_p = getDerMatchExponent(alpha_bar,alpha);
                choE(i) = prod(E(n_E*(j-1)+1:n_E*j).^gamma_p);
                
            case 'zero-cumulants'
                YI = Y(I');
                % generate all possible partitions in YI
                B = partitions(YI);
                % generating cumulant corresponding to this moment
                K = sym(0);
                for b = 1:length(B)
                    parts = [];
                    for ib = 1:length(B{b})
                        parts = [parts,prod(B{b}{ib})];
                    end
                    parts = mysubs(parts,...
                        [hoE_funY(end:-1:1);E_funY(end:-1:1)],...
                        [hoE(n_hoE*j:-1:n_hoE*(j-1)+1);E(n_E*j:-1:n_E*(j-1)+1)]);
                    npart = length(parts);
                    K = K + factorial(npart-1) * (-1)^(npart-1) * prod(parts);
                end
                choE(i) = solve(K,hoE_used(i));
                
            otherwise
                error('This option is not available.');
        end
    end
    if isfield(system,'input')
        while (~isempty(setdiff(symvar(choE),[p;c;E;system.time;system.input.variable])))
            choE = subs(choE,hoE_used,choE,0);
        end
    else
        while (~isempty(setdiff(symvar(choE),[p;c;E;system.time])))
            choE = subs(choE,hoE_used,choE,0);
        end
    end
    % Substitution
    if isfield(system,'input')
        while (~isempty(setdiff(symvar(dpdt),[p;c;E;system.time;system.input.variable])))
            dpdt  = subs(dpdt,hoE_used,choE,0);
        end
    else
        while (~isempty(setdiff(symvar(dpdt),[p;c;E;system.time])))
            dpdt  = subs(dpdt,hoE_used,choE,0);
        end
    end
    if isfield(system,'input')
        while (~isempty(setdiff(symvar(dEpdt),[p;c;E;system.time;system.input.variable])))
            dEpdt = subs(dEpdt,hoE_used,choE,0);
        end
    else
        while (~isempty(setdiff(symvar(dEpdt),[p;c;E;system.time])))
            dEpdt = subs(dEpdt,hoE_used,choE,0);
        end
    end
end
dpdt = simplify(dpdt);
dEpdt = simplify(dEpdt);
% Derivative vector
dMdt = sym(zeros((1+n_E)*n_z_index,1));
for i = 1:length(p)
    dMdt(i) = sym(['d' char(p(i)) 'dt']);
end
for i = length(p)+1:length(dMdt)
    dMdt(i) = sym(['d' char(E(i-length(p))) 'dt']);
end
%% CONSTRUCTION OF VECTOR FIELD AND MASS MATRIX
VFsym = dpdt;
MMsym = sym(ones(n_z_index,1));
for i = 1:n_z_index
    VFsym(n_z_index + [(i-1)*n_E+1:i*n_E]) = dEpdt((i-1)*n_E+1:i*n_E) - E((i-1)*n_E+1:i*n_E).*dpdt(i);
    MMsym(n_z_index + [(i-1)*n_E+1:i*n_E]) = p(i);
end
if isfield(system,'input')
    VFsym = subs(VFsym,system.input.variable,system.input.function,0);
    VFsym = subs(VFsym,system.parameter.variable,c,0);
    MMsym = subs(MMsym,system.input.variable,system.input.function,0);
    MMsym = subs(MMsym,system.parameter.variable,c,0);
end
VFsym = simplify(VFsym);
VFsym = collect(VFsym);

RESsym = -MMsym.*dMdt + VFsym;
%% CONVERSION TO EXECUTABLE FUNCTION
str_VF = '@(x,theta) [';
str_MM = '@(x,theta) [';
str_RES = '@(x,dx,theta) [';
for i = 1:length(VFsym)-1
    str_VF = [str_VF char(VFsym(i)) ';'];
    str_MM = [str_MM char(MMsym(i)) ';'];
    str_RES = [str_RES char(RESsym(i)) ';'];
end

str_VF = [str_VF char(VFsym(end)) ']'];
% str_MM = [str_MM char(MMsym(end)) '],0,' num2str(length(VFsym),'%d') ',' num2str(length(VFsym),'%d') ')'];
str_MM = [str_MM char(MMsym(end)) ']'];
str_RES = [str_RES char(RESsym(end)) ']'];
% Exchange variables (backwards to avoid errors)
for i = length(E):-1:1
    str_VF = strrep(str_VF,char(E(i)),['x(' num2str(i+n_z_index,'%d') ')']);
    str_MM = strrep(str_MM,char(E(i)),['x(' num2str(i+n_z_index,'%d') ')']);
    str_RES = strrep(str_RES,char(dMdt(i+n_z_index)),['dx(' num2str(i+n_z_index,'%d') ')']);
    str_RES = strrep(str_RES,char(E(i)),['x(' num2str(i+n_z_index,'%d') ')']);
    
end
for i = length(p):-1:1
    str_VF = strrep(str_VF,char(p(i)),['x(' num2str(i,'%d') ')']);
    str_MM = strrep(str_MM,char(p(i)),['x(' num2str(i,'%d') ') + 1e-6']);
    str_RES = strrep(str_RES,char(dMdt(i)),['dx(' num2str(i,'%d') ')']);
    str_RES = strrep(str_RES,char(p(i)),['x(' num2str(i,'%d') ')']);
end
% Exchange parameters (backwards to avoid errors)
for i = n_c:-1:1
    str_VF = strrep(str_VF,char(c(i)),['theta(' num2str(i,'%d') ')']);
    str_MM = strrep(str_MM,char(c(i)),['theta(' num2str(i,'%d') ')']);
    str_RES = strrep(str_RES,char(c(i)),['theta(' num2str(i,'%d') ')']);
end
% Construct function
VF = eval(str_VF);
MM = eval(str_MM);
RES = eval(str_RES);

%% CONSTRUCT FUNCTION HANDLE
clear([options.filepath options.filename '.m']);
% Open file
fid = fopen([options.filepath options.filename '.m'],'w');
% Construct string
str_FUN = ['function [res,flag,new_data] = ' options.filename '(t,x,dx,data) \n\n'...
    'theta = data.params;\n\n'...
    'res =' strrep(str_RES(14:end),';',';...\n       ') ';\n\n'...
    'flag = 0;\n'...
    'new_data = [];'];
% Write to file
fprintf(fid,str_FUN);
fclose(fid);
% Rehash to ensure that function is known / used
rehash
%% CONSTRUCT FUNCTION HANDLE - SIMULATION
clear([options.filepath 'simUncentered_' options.filename '.m']);
% Open file
fid = fopen([options.filepath 'simUncentered_' options.filename '.m'],'w');
% Construct string
str_FUN = ['function [p_cMM,cmu_cMM,cC_cMM] = simUncentered_' options.filename '(cMM,t,theta,x0,dx0) \n\n'...
    'data.params = theta;\n\n'...
    '%% Determine dimension\n'...
    'n_y = size(cMM.state.stochatic.FSP_index,1);\n'...
    'n_z = length(cMM.state.expectation.state_index);\n'...
    'n_C = size(cMM.state.expectation.E_index,1);\n'...
    '%% Set solver options\n'...
    'options = IDASetOptions(''RelTol'',1e-10,...\n'...
    '''AbsTol'',1e-10,...\n'...
    '''VariableTypes'',ones(n_y*(1+n_C),1),...\n'...
    '''suppressAlgVars'',''on'',...\n'...
    '''MaxNumSteps'', 10^4,...\n'...
    '''LinearSolver'',''Dense'');\n'...
    '%% Initiallize IDAS and solve ODE\n'...
    'IDAInit(@(t,x,dx) ' options.filename '(t,x,dx,data),0,x0,dx0,options);\n'...
    '[status,ty,y] = IDASolve(t(2:end),''Normal'');\n'...
    'cM = [x0'';y''];\n'...
    '%% Free memory\n'...
    'IDAFree;\n\n'...
    '%% Reorder result\n'...
    'for iy = 1:n_y\n'...
    '       p_cMM{iy}   = cM(:,iy);\n'...
    '       cmu_cMM{iy} = cM(:,(n_C)*(iy-1)+n_y    +[1:n_z]);\n'...
    '       cC_cMM{iy}  = cM(:,(n_C)*(iy-1)+n_y+[n_z+1:n_C]);\n'...
    'end\n'];
% Write to file
fprintf(fid,str_FUN);
fclose(fid);
% Rehash to ensure that function is known / used
rehash
%% ASSEMBLE OUTPUT
% Stochastic state
cMFSP.state.stochatic.name = system.state.name(1:n_z);
cMFSP.state.stochatic.variable = system.state.variable(1:n_z);
cMFSP.state.stochatic.xmin = system.state.xmin(1:n_z);
cMFSP.state.stochatic.xmax = system.state.xmax(1:n_z);
cMFSP.state.stochatic.state_index = z_ind;
cMFSP.state.stochatic.FSP_index = z_index;
% Expectation state
cMFSP.state.expectation.name = system.state.name(n_z+1:end);
cMFSP.state.expectation.variable = system.state.variable(n_z+1:end);
cMFSP.state.expectation.xmin = system.state.xmin(n_z+1:end);
cMFSP.state.expectation.xmax = system.state.xmax(n_z+1:end);
cMFSP.state.expectation.state_index = y_ind;
cMFSP.state.expectation.E_index = E_ind;
% Symbolic
cMFSP.state.sym.stochatic = p;
cMFSP.state.sym.expectation = E;
cMFSP.state.sym.moment = [p;E];
cMFSP.state.sym.ho_moment = hoE_used;
cMFSP.state.sym.ho_moment_closure = choE;
cMFSP.state.sym.dMdt = dMdt;

% Symbolic
cMFSP.parameters.sym.c = c;
cMFSP.system.stoichiometry = S;
% State of cMFSP
cMFSP.derivatives.sym.dpdt = dpdt;
cMFSP.derivatives.sym.dEpdt = dEpdt;
cMFSP.derivatives.sym.VF = VFsym;
cMFSP.derivatives.sym.mass_matrix = MMsym;
cMFSP.derivatives.sym.res = RESsym;
cMFSP.derivatives.fun.VF = VF;
cMFSP.derivatives.fun.mass_matrix = MM;
cMFSP.derivatives.fun.res = RES;

