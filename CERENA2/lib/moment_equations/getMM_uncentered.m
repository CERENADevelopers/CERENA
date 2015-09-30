% function MM = getMM_uncentered(system,options)
function MM = getMM_uncentered(varargin)

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
options.filename = 'simUncentered_MM';
if nargin >= 2
    options = setdefault(varargin{2},options);
end
if ~isfield(options,'taylor_order')
    options.taylor_order = [];
end
%% INITIALIZATION (1)
n_s = length(system.state.variable);
if isfield(system,'output')
    n_sy = length(system.output.variable);
end
n_r = length(system.reaction);
n_c = length(system.parameter.variable);

%% INITIALIZATION (2)
if options.moment_order <4
    ro_max = 2;
else
    ro_max = options.moment_order - 1;
end
% Generation of I-index
[I,n_I] = getMomentIndexSet(n_s,options.moment_order+ro_max);
hoI = I(find(sum(I~=0,2)> options.moment_order),:);
I   = I(find(sum(I~=0,2)<=options.moment_order),:);
% Generation of alpha-index
alpha   = convertI2alpha(  I,n_s);
hoalpha = convertI2alpha(hoI,n_s);
% Moments
M   = getM(I);
hoM = getM(hoI);
% Dimensionalities
n_M = length(M);
n_hoM = length(hoM);

%outputs
%[Iyi,n_Iyi] = getMomentIndexSet(1,options.moment_order);
if isfield(system,'output')
    Iy = zeros(n_sy*options.moment_order_output,options.moment_order_output);
    for i = 1:n_sy
        for j = 1:options.moment_order_output
            Iy(n_sy*(j-1)+i,end:-1:(end-j+1)) = i*ones(1,j);
        end
    end
    n_Iy = size(Iy,1);
    % All moments
    My = getM(Iy);
end

% Parameters
c = sym(zeros(n_c,1));
for i = 1:n_c
    c(i) = sym(['c' num2str(i,'%d')]);
end

% States
X = sym(zeros(n_s,1));
% Loop: states
for i = 1:n_s
    X(i) = sym(['X' num2str(i,'%d')]);
end
XM   = prod(repmat(transpose(X),  n_M,1).^alpha  ,2);
hoXM = prod(repmat(transpose(X),n_hoM,1).^hoalpha,2);
% Initial conditions
if options.moment_order == 1
    M0 = [system.state.mu0];
    M0 = subs(M0,system.parameter.variable,c,0);
elseif options.moment_order > 1
    M0 = [system.state.mu0;system.state.C0];
    M0 = sym(subs(M0,system.parameter.variable,c,0));
end
%% PROPENSITIES
for r = 1:n_r
    % Propensities
    w{r} = system.reaction(r).propensity;
    w{r} = subs(w{r},system.parameter.variable,c,0);
    w{r} = subs(w{r},system.state.variable,X,0);
    if isfield(system.reaction(r),'propensityType')
        switch system.reaction(r).propensityType
            case 'polynomial'
                w{r} = w{r};
            case 'non-polynomial'
                w{r} = getTaylorExpansion(w{r},X,M,options.taylor_order);
        end
    else %interpret as polynomial
        w{r} = w{r};
    end
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
%% ASSIGN OUTPUT MAP AND COMPUTE DERIVATIVES
if isfield(system,'output')
    for o = 1:n_sy
        h{o} = system.output.function(o);
        h{o} = subs(h{o},system.parameter.variable,c,0);
        h{o} = subs(h{o},system.state.variable,X,0);
    end
end
%% CONSTRUCTION OF MOMENT EQUATION
dMdt = sym(zeros(n_M,1));
% Loop: moments
for i = 1:n_M
    % Assignment of I-index
    Ii = I(i,find(I(i,:)~=0));
    % Loop: reactions
    for r = 1:n_r
        dMdt(i) = dMdt(i) + expand((prod(X(Ii')+S(Ii,r)) - prod(X(Ii)))*w{r});
    end
end

%% SUBSTITUDE EXPECTED VALUES
dMdt = mysubs(dMdt,hoXM(end:-1:1,:),hoM(end:-1:1));
dMdt = mysubs(dMdt,  XM(end:-1:1,:),  M(end:-1:1));
%% CONSTRUCT OUTPUT MAP
if isfield(system,'output')
    % Initialize
    H = sym(zeros(n_Iy,1));
    % Loop: Output moments
    for i = 1:length(H)
        % Assignment of I-index and alpha-index
        Iyi = Iy(i,find(Iy(i,:)~=0));
        alphaiy = zeros(n_sy,1);
        hM = sym(zeros(n_sy),1);
        for k = 1:n_sy
            alphaiy(k) = sum(Iyi == k);
            hM(k) = h{k}^alphaiy(k);
        end
        H(i) = subs(prod(hM),XM(end:-1:1,:),M(end:-1:1),0);
    end
    H = simplify(H);
    disp('output map done.')
else
    H = [];
end
%% MOMENT CLOSURE
% Determine moments of order > options.moment_order
hoM_used = transpose(setdiff(symvar(dMdt),[c;M;system.time]));
if isfield(system,'input')
    hoM_used = setdiff(symvar(hoM_used),system.input.variable);
end
if ~isempty(hoM_used)
    disp('Moment-closure used!')
    
    % Closure
    hoM_closure = sym([]);
    % Loop: higher order moments
    for i = 1:length(hoM_used)
        j = find(hoM==hoM_used(i));
        Ii = hoI(j,(find(hoI(j,:)~=0)));
        switch options.moment_closure
            case 'mean-field'
                hoM_closure(i) = prod(M(Ii));
            case 'low-dispersion'
                Xi = X(Ii);
                Mi = M(Ii);
                eq = expand(prod(Xi(:) - Mi(:)));
                eq = mysubs(eq,hoXM(end:-1:1),hoM(end:-1:1));
                eq = mysubs(eq,  XM(end:-1:1),  M(end:-1:1));
                hoM_closure(i) = solve(eq,hoM_used(i));
            case 'derivative-matching'
                %             I = hoE_ind(k,:);
                alpha_bar = convertI2alpha(Ii,n_s);
                alpha = convertI2alpha(I,n_s);
                
                
                gamma_p = getDerMatchExponent(alpha_bar,alpha);
                hoM_closure(i) = prod(M.^gamma_p);
                
            case 'zero-cumulants'
                XI = X(Ii');
                % generate all possible partitions in YI
                B = partitions(XI);
                % generating cumulant corresponding to this moment
                K = sym(0);
                for b = 1:length(B)
                    parts = [];
                    for ib = 1:length(B{b})
                        parts = [parts,prod(B{b}{ib})];
                    end
                    parts = mysubs(parts,...
                        [hoXM(end:-1:1);XM(end:-1:1)],...
                        [hoM(end:-1:1);M(end:-1:1)]);
                    npart = length(parts);
                    K = K + factorial(npart-1) * (-1)^(npart-1) * prod(parts);
                end
                hoM_closure(i) = solve(K,hoM_used(i));
            otherwise
                error('This option is not available.');
        end
    end
    if isfield(system,'input')
        while (~isempty(setdiff(symvar(hoM_closure),[c;M;system.time;system.input.variable])))
            hoM_closure = subs(hoM_closure,hoM_used,hoM_closure,0);
        end
    else
        while (~isempty(setdiff(symvar(hoM_closure),[c;M;system.time])))
            hoM_closure = subs(hoM_closure,hoM_used,hoM_closure,0);
        end
    end
    % Substitutions
    if isfield(system,'input')
        while (~isempty(setdiff(symvar(dMdt),[c;M;system.time;system.input.variable])))
            dMdt = subs(dMdt,hoM_used,hoM_closure,0);
        end
    else
        while (~isempty(setdiff(symvar(dMdt),[c;M;system.time])))
            dMdt = subs(dMdt,hoM_used,hoM_closure,0);
        end
    end
end
dMdt = simplify(dMdt);
%% EVALUATION OF JACOBIAN
f = dMdt;
if isfield(system,'input')
    f = subs(f,system.input.variable,system.input.function,0);
    f = subs(f,system.state.variable,X,0);
    f = subs(f,XM,M,0);
    f = subs(f,system.parameter.variable,c,0);
end
J = jacobian(f,M);
dfdc = jacobian(f,c);
dHdc = jacobian(H,c);
dHdx = jacobian(H,M);
dM0dc = jacobian(M0,c);
%% CONVERSION TO FUNCTION HANDLE
% Write function
str_f     = '@(t,x,theta) [';
str_M0    = '@(theta) [';
str_dfdc  = '@(t,x,dxdt,theta) [';
str_dM0dc = '@(theta) [';
str_J     = '@(t,x,theta) [';
for i1 = 1:length(M)
    % f and M0
    str_f = [str_f char(f(i1))];
    str_M0 = [str_M0 char(M0(i1))];
    % J
    for i2 = 1:length(M)-1
        str_J = [str_J char(J(i1,i2)) ','];
    end
    str_J = [str_J char(J(i1,length(M)))];
    % dfdc
    for i3 = 1:length(c)-1
        str_dfdc  = [str_dfdc  char(dfdc(i1,i3))  ','];
        str_dM0dc = [str_dM0dc char(dM0dc(i1,i3)) ','];
    end
    str_dfdc  = [str_dfdc  char(dfdc(i1,length(c))) ];
    str_dM0dc = [str_dM0dc char(dM0dc(i1,length(c)))];
    % Closing of expression
    if i1 < length(M)
        str_f     = [str_f     ';'];
        str_M0    = [str_M0    ';'];
        str_J     = [str_J     ';'];
        str_dfdc  = [str_dfdc  ';'];
        str_dM0dc = [str_dM0dc ';'];
    else
        str_f     = [str_f     ']'];
        str_M0    = [str_M0    ']'];
        str_J     = [str_J     ']'];
        str_dfdc  = [str_dfdc  ']'];
        str_dM0dc = [str_dM0dc ']'];
    end
end

str_H = '@(t,x,theta) [';
str_dHdx = '@(t,x,theta) [';
str_dHdc = '@(t,x,theta) [';
if isfield(system,'output')
    for i1 = 1:length(My)
        % H
        str_H = [str_H char(H(i1))];
        % dHdx
        for i2 = 1:length(M)-1
            str_dHdx = [str_dHdx char(dHdx(i1,i2)) ','];
        end
        str_dHdx = [str_dHdx char(dHdx(i1,length(My)))];
        % dHdc
        for i3 = 1:length(c)-1
            str_dHdc = [str_dHdc char(dHdc(i1,i3)) ','];
        end
        str_dHdc = [str_dHdc char(dHdc(i1,length(c)))];
        % Closing of expression
        if i1 < length(My)
            str_H = [str_H ','];
            str_dHdx = [str_dHdx ';'];
            str_dHdc = [str_dHdc ';'];
        else
            str_H = [str_H ']'];
            str_dHdx = [str_dHdx ']'];
            str_dHdc = [str_dHdc ']'];
        end
    end
else
    str_H = [str_H ']'];
    str_dHdx = [str_dHdx ']'];
    str_dHdc = [str_dHdc ']'];
end

% Exchange variables (backwards to avoid errors)
for i = length(M):-1:1
    str_f     = strrep(str_f    ,char(M(i)),['x('   num2str(i,'%d') ')']);
    str_J     = strrep(str_J    ,char(M(i)),['x('   num2str(i,'%d') ')']);
    str_dfdc  = strrep(str_dfdc ,char(M(i)),['x('   num2str(i,'%d') ')']);
    str_H     = strrep(str_H    ,char(M(i)),['x(:,' num2str(i,'%d') ')']);
    str_dHdx  = strrep(str_dHdx ,char(M(i)),['x(:,' num2str(i,'%d') ')']);
    str_dHdc  = strrep(str_dHdc ,char(M(i)),['x(:,' num2str(i,'%d') ')']);

end
% Exchange parameters (backwards to avoid errors)
for i = n_c:-1:1
    str_f     = strrep(str_f    ,char(c(i)),['theta(' num2str(i,'%d') ')']);
    str_M0    = strrep(str_M0   ,char(c(i)),['theta(' num2str(i,'%d') ')']);
    str_J     = strrep(str_J    ,char(c(i)),['theta(' num2str(i,'%d') ')']);
    str_dfdc  = strrep(str_dfdc ,char(c(i)),['theta(' num2str(i,'%d') ')']);
    str_dM0dc = strrep(str_dM0dc,char(c(i)),['theta(' num2str(i,'%d') ')']);
    str_H     = strrep(str_H    ,char(c(i)),['theta(' num2str(i,'%d') ')']);
    str_dHdx  = strrep(str_dHdx ,char(c(i)),['theta(' num2str(i,'%d') ')']);
    str_dHdc  = strrep(str_dHdc ,char(c(i)),['theta(' num2str(i,'%d') ')']);

end
% substitute time
str_f     = strrep(strrep(str_f    ,' ',''),char(system.time),'t');
str_M0    = strrep(strrep(str_M0   ,' ',''),char(system.time),'t');
str_J     = strrep(strrep(str_J    ,' ',''),char(system.time),'t');
str_dfdc  = strrep(strrep(str_dfdc ,' ',''),char(system.time),'t');
str_dM0dc = strrep(strrep(str_dM0dc,' ',''),char(system.time),'t');
str_H     = strrep(strrep(str_H    ,' ',''),char(system.time),'t');
str_dHdx  = strrep(strrep(str_dHdx ,' ',''),char(system.time),'t');
str_dHdc  = strrep(strrep(str_dHdc ,' ',''),char(system.time),'t');

% Construct function
MM.fun.f     = eval(str_f    );
MM.fun.M0    = eval(str_M0   );
MM.fun.J     = eval(str_J    );
MM.fun.dfdc  = eval(str_dfdc );
MM.fun.dM0dc = eval(str_dM0dc);
MM.fun.H     = eval(str_H    );
MM.fun.dHdx  = eval(str_dHdx );
MM.fun.dHdc  = eval(str_dHdc );

%% CONSTRUCT SIMULATION FILE
clear(['simUncentered_' options.filename '.m']);
% Open file
fid = fopen(['simUncentered_' options.filename '.m'],'w');
% Construct string
str_FUN_Sim = ['%% function [X,Y,SX,SY] =  simUncentered_' options.filename '(t,theta,x0) \n'...
    'function varargout = simUncentered_' options.filename '(varargin) \n\n'...
    't = varargin{1};\n'...
    'theta = varargin{2};\n'...
    'x0 = [];\n'...
    'if nargin>2\n'...
    '    x0 = varargin{3};\n'...
    'end\n'...
    'data.theta = theta;\n'...
    '%% Set solver options \n'...
    'options_CVode = CVodeSetOptions(''RelTol'',1e-6,...\n'...
    '                                ''AbsTol'',1e-6,...\n'...
    '                                ''MaxNumSteps'',10^6,...\n'...
    '                                ''JacobianFn'',@jacfn,...\n'...
    '                                ''Userdata'',data);\n'...
    'options_CVodes = CVodeSensSetOptions(''method'',''Simultaneous'',...\n'...
    '                                     ''ErrControl'',true,...\n'...
    '                                     ''ParamScales'',1:length(theta));\n'...
    '\n'...
    '%% Initial conditions\n'...
    'if isempty(x0)\n'...
    '    x0 = x0fun(theta);\n'...
    'end\n'...
    'if nargout >= 3\n'...
    '    sx0 = sx0fun(theta);\n'...
    'end\n'...
    '\n'...
    '%% Initialization of CVode\n'...
    'CVodeInit(@rhs,''BDF'',''Newton'',0,x0,options_CVode);\n'...
    'if nargout >= 3\n'...
    '    CVodeSensInit(length(theta),@rhsS,sx0,options_CVodes);\n'...
    'end\n'...
    '\n'...
    '%% Simulation\n'...
    'if nargout <= 2\n'...
    '    [status,~,x] = CVode(t(2:end),''Normal'');\n'...
    '    X = [x0'';x''];\n'...
    '    Y = rhsO(t,X,theta);\n'...
    'else\n'...
    '    [status,~,x,sx] = CVode(t(2:end),''Normal'');\n'...
    '    X = [x0'';x''];\n'...
    '    Y = rhsO(t,X,theta);\n'...
    '    SX = zeros(length(t),length(x0),length(theta));\n'...
    '    SX(1,:,:) = sx0;\n'...
    '    SX(2:end,:,:) = permute(sx,[3,1,2]);\n'...
    '    SY = rhsOS(t,X,SX,Y,theta);\n'...
    'end\n'...
    '\n'...
    '%% Free memory\n'...
    'CVodeFree;\n'...
    '\n'...
    '%% Assign output\n'...
    'varargout{1} = X;\n'...
    'if nargout >= 2\n'...
    '    varargout{2} = Y;\n'...
    'end\n'...
    'if nargout >= 3\n'...
    '    varargout{3} = SX;\n'...
    'end\n'...
    'if nargout >= 4\n'...
    '    varargout{4} = SY;\n'...
    'end\n'...
    'if nargout >= 5\n'...
    '    error(''Too many output arguments.'');\n'...
    'end\n'...
    '\n\n'...
    '%%%% RIGHT-HAND SIDE\n'...
    'function [dxdt,flag,new_data] = rhs(t,x,data) \n\n'...
    'theta = data.theta;\n'...
    'dxdt = ' strrep(str_f(13:end),';',';...\n        ') ';\n\n'...
    'flag = 0;\n'...
    'new_data = [];'...
    '\n\n'...
    '%%%% RIGHT-HAND SIDE OF SENSITIVITIES\n'...
    'function [dsxdt,flag,new_data] = rhsS(t,x,dxdt,sx,data) \n\n'...
    'theta = data.theta;\n'...
    'J = jacfn(t,x,dxdt,data);\n'...
    'dfdtheta = ' strrep(str_dfdc(18:end),';',';...\n              ') ';\n\n'...
    'dsxdt = J*sx + dfdtheta;\n\n'...
    'flag = 0;\n'...
    'new_data = [];'...
    '\n\n'...
    '%%%% JACOBIAN\n'...
    'function [J,flag,new_data] = jacfn(t,x,dxdt,data) \n\n'...
    'theta = data.theta;\n'...
    'J = ' strrep(str_J(13:end),';',';...\n     ') ';\n'...
    'flag = 0;\n'...
    'new_data = [];'...
    '\n\n'...
    '%%%% OUTPUT MAP\n'...
    'function y = rhsO(t,x,theta) \n\n'...
    'y = ' str_H(13:end) ';\n'...
    '\n\n'...
    '%%%% OUTPUT MAP OF SENSITIVITIES\n'...
    'function sy = rhsOS(t,x,sx,y,theta) \n\n'...
    'sy = zeros(length(t),size(y,2),length(theta));\n'...
    'for k = 1:length(t)\n'...
    '    dHdx = ' strrep(strrep(str_dHdx(13:end),':','k'),';',';...\n            ') ';\n'...
    '    dHdtheta = ' strrep(strrep(str_dHdc(13:end),':','k'),';',';...\n                ') ';\n'...
    '    sy(k,:,:) = dHdx*squeeze(sx(k,:,:)) + dHdtheta;\n'...
    'end\n'...
    '\n\n'...
    '%%%% INITIAL CONDITIONS FOR STATE\n'...
    'function x0 = x0fun(theta) \n\n'...
    'x0 = ' strrep(str_M0(9:end),';',';...\n      ') ';\n'...
    '\n\n'...
    '%%%% INITIAL CONDITIONS FOR STATE SENSITIVITY\n'...
    'function sx0 = sx0fun(theta) \n\n'...
    'sx0 = ' strrep(str_dM0dc(9:end),';',';...\n       ') ';\n'...
    '\n\n'...
    ];
% Write to file
fprintf(fid,str_FUN_Sim);
fclose(fid);

% Rehash to ensure that function is known / used
rehash

%% ASSEMBLE OUTPUT
MM.type  = 'uncentered';
MM.order = options.moment_order;
MM.closure = options.moment_closure;
MM.system.stoichiometry = S;
MM.sym.state.order = I;
MM.sym.state.moments = M;
MM.sym.state.ho_moments = hoM_used;
MM.sym.state.ho_moments_closure = hoM_closure;
MM.sym.state.derivative = dMdt;
MM.sym.state.M0 = M0;
if isfield(system,'output')
    MM.sym.output.order = Iy;
    MM.sym.output.moments = My;
    MM.sym.output.function = H;
end
