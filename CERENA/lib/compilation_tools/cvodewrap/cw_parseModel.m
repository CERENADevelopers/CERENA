% parseModel parses the model definition struct which contains the symbolic description of the ODE system and
% symbolically derives all expressions required by CVODES
%
% USAGE:
% ======
% [struct] = parseModel(modelstruct)
%
% INPUTS:
% =======
% modelstruct | is the struct which contains the symbolic description of the ODE system
%
% OUTPUS:
% =======
% struct | is the struct which can be subsequently used by generateC and generateM

function [struct] = cw_parseModel(struct)

if(~isstruct(struct))
    error('The provided symfun did not return a struct')
end

if(~isfield(struct,'sym'))
    error('Model struct is missing the sym field (.sym)!')
end

% check whether sym is properly defined
if(~isfield(struct.sym,'x'))
    error('Model struct is missing the definition of the state vector x (.sym.x)!')
end
if(~isfield(struct.sym,'xdot'))
    error('Model struct is missing the definition of the right hand side xdot (.sym.xdot)!')
end
if(~isfield(struct.sym,'p'))
    error('Model struct is missing the definition of the parameter vector p (.sym.p)!')
end
if(~isfield(struct.sym,'x0'))
    error('Model struct is missing the definition of the vector of initial conditions x0 (.sym.x0)!')
end
if(~isfield(struct.sym,'y'))
    error('Model struct is missing the definition of the vector of observables y (.sym.y)!')
end
if(~isfield(struct,'atol'))
    struct.atol = 1e-8;
end
if(~isfield(struct,'rtol'))
    struct.rtol = 1e-8;
end
if(~isfield(struct,'maxsteps'))
    struct.maxsteps = 1e4;
end
if(~isfield(struct,'noreorder'))
    struct.noreorder = 0;
end
if(~isfield(struct,'noadjoint'))
    struct.noadjoint = 0;
end
if(~isfield(struct,'noforward'))
    struct.noforward = 0;
end
if(size(struct.sym.x,1)<size(struct.sym.x,2))
    struct.sym.x = transpose(struct.sym.x);
end
if(size(struct.sym.xdot,1)<size(struct.sym.xdot,2))
    struct.sym.xdot = transpose(struct.sym.xdot);
end
if(size(struct.sym.x0,1)<size(struct.sym.x0,2))
    struct.sym.x0 = transpose(struct.sym.x0);
end
if(~all([size(struct.sym.x,2)==size(struct.sym.xdot,2),size(struct.sym.xdot,2)==size(struct.sym.x0,2)]))
    error('Sizes of x0, xdot and x do not agree!')
end

% complete optional fields
if(~isfield(struct.sym,'u'))
    struct.sym.u = sym.empty(0,0);
end
if(~isfield(struct.sym,'k'))
    struct.sym.k = sym.empty(0,0);
end
if(~isfield(struct.sym,'root'))
    struct.sym.root = sym.empty(0,1);
end
if(~isfield(struct,'t0'))
    struct.t0 = 0;
end

if(any(ismember(struct.sym.k,struct.sym.p)))
    error(['Invalid Model: ' char(struct.sym.k(find(ismember(struct.sym.k,struct.sym.p),1))) ' is contained in both p and k!'])
end

fprintf('reordering | ')
% optimal reordering (aproximate minimum degree jacobian)
J=jacobian(struct.sym.xdot,struct.sym.x);
M = double(logical(J~=sym(zeros(size(J)))));

if(~struct.noreorder)
    opts.dense = 100;
    opts.agressive = 1;
    r = amd(M,opts);
else
    r = 1:length(struct.sym.x);
end

[ubw,lbw] = cw_bandwidth(M(r,r));

struct.ubw = ubw;
struct.lbw = lbw;
struct.r = r;
% inverse permutation
struct.rt(r) = 1:length(r);
struct.nnz = length(find(M(:)));

struct.sym.x = struct.sym.x(r);
struct.sym.xdot = struct.sym.xdot(r);
struct.sym.x0 = struct.sym.x0(r);

nx = length(struct.sym.x);
np = length(struct.sym.p);
nk = length(struct.sym.k);
ny = length(struct.sym.y);
%remove zero-roots
ir = 1;
while ir <= length(struct.sym.root)
    if(isequaln(struct.sym.root(ir),0))
        struct.sym.root(ir) = [];
    else 
        ir = ir + 1;
    end
end
nr = length(struct.sym.root);

struct.nx = nx;
struct.ny = ny;
struct.nr = nr;

% short strings
xs = cell(nx,1);
vs = cell(nx,1);
ps = cell(np,1);
ks = cell(nk,1);
sx = cell(nx,np);
xBs = cell(nx,ny);

for j=1:nx
    xs{j} = sprintf('x[%i]',j-1);
    vs{j} = sprintf('v[%i]',j-1);
end
for j=1:np
    ps{j} = sprintf('p[%i]',j-1);
end
for j=1:nk
    ks{j} = sprintf('k[%i]',j-1);
end
for j = 1:nx
    for i = 1:np
        sx{j,i} = sprintf('sx[%i]', j-1);
    end
end

for j = 1:nx
    for i = 1:ny
        xBs{j,i} = sprintf('xB[%i]', j-1 + (i-1)*nx);
    end
end

% transform into syms
strsym.xs = sym(xs);
strsym.vs = sym(vs);
strsym.ps = sym(ps);
strsym.ks = sym(ks);
strsym.sx = sym(sx);
strsym.xBs = sym(xBs);

% replace states by short strings
struct.sym.xdot = mysubs(struct.sym.xdot,struct.sym.x,strsym.xs);
struct.sym.y = mysubs(struct.sym.y,struct.sym.x,strsym.xs);
struct.sym.root = mysubs(struct.sym.root,struct.sym.x,strsym.xs);

struct.sym.xdot = mysubs(struct.sym.xdot,struct.sym.p,strsym.ps);
struct.sym.y = mysubs(struct.sym.y,struct.sym.p,strsym.ps);
struct.sym.x0 = mysubs(struct.sym.x0,struct.sym.p,strsym.ps);
struct.sym.root = mysubs(struct.sym.root,struct.sym.p,strsym.ps);

struct.sym.xdot = mysubs(struct.sym.xdot,struct.sym.k,strsym.ks);
struct.sym.y = mysubs(struct.sym.y,struct.sym.k,strsym.ks);
struct.sym.x0 = mysubs(struct.sym.x0,struct.sym.k,strsym.ks);
struct.sym.root = mysubs(struct.sym.root,struct.sym.k,strsym.ks);

% compute derivatives
fprintf('J | ')
struct.sym.J=jacobian(struct.sym.xdot,strsym.xs);
fprintf('JB | ')
struct.sym.JB=-transpose(struct.sym.J);

fprintf('dxdotdp | ')
struct.sym.dxdotdp=jacobian(struct.sym.xdot,strsym.ps);
if(~struct.noforward)
    fprintf('sx | ')
    struct.sym.sx=struct.sym.J*strsym.sx+struct.sym.dxdotdp;
    fprintf('sx0 | ')
    struct.sym.sx0=jacobian(struct.sym.x0,strsym.ps);
else
    struct.sym.sx=sym(zeros(nx,np));
    struct.sym.sx0=sym(zeros(nx,np));
end

fprintf('dydx | ')
struct.sym.dydx=jacobian(struct.sym.y,strsym.xs);
if(~struct.noforward)
    fprintf('dydp | ')
    struct.sym.dydp=jacobian(struct.sym.y,strsym.ps);
    fprintf('sy | ')
    struct.sym.sy=struct.sym.dydp + struct.sym.dydx*strsym.sx ;
else
    struct.sym.dydp=sym(zeros(ny,np));
    struct.sym.sy=sym(zeros(ny,np));
end

if(nr>0 && not(struct.noforward))
    fprintf('drootdx | ')
    struct.sym.drootdx = jacobian(struct.sym.root,strsym.xs);
    fprintf('drootdt | ')
    struct.sym.drootdt = diff(struct.sym.root,sym('t')) + struct.sym.drootdx*struct.sym.xdot;
    fprintf('drootdp | ')
    struct.sym.drootdp = jacobian(struct.sym.root,strsym.ps) + struct.sym.drootdx*strsym.sx;
    fprintf('sroot | ')
    for ir = 1:nr
        struct.sym.sroot(ir,:) = -(struct.sym.drootdt(ir,:))\(struct.sym.drootdp(ir,:));
    end
    struct.sym.s2root =  sym(zeros(nr,np,np));
    struct.sym.srootval = struct.sym.drootdp;
    struct.sym.s2rootval =  sym(zeros(nr,np,np));
else
    struct.sym.drootdx = sym.empty(0,nx);
    struct.sym.drootdt = sym.empty(0,1);
    struct.sym.drootdp = sym.empty(0,np);
    struct.sym.sroot = sym.empty(0,np);
    struct.sym.s2root =  sym.empty(0,np,np);
    struct.sym.s2rootval = sym.empty(0,np,np);
end

% spils solvers
fprintf('Jv | ')
struct.sym.Jv = struct.sym.J*strsym.vs;
fprintf('JvB | ')
struct.sym.JvB = -transpose(struct.sym.J)*strsym.vs;

% adjoint sensitivities
if(~struct.noadjoint)
    fprintf('xBdot | ')
    struct.sym.xBdot = -transpose(struct.sym.J)*strsym.xBs;
    
    fprintf('int | ')
    struct.sym.int = -transpose(strsym.xBs)*struct.sym.dxdotdp;
    
    fprintf('mu0 | ')
    struct.sym.mu0 = transpose(struct.sym.dydx);
else
    struct.sym.xBdot = sym(zeros(nx,ny));
    
    struct.sym.int = sym(zeros(ny,np));
    
    struct.sym.mu0 = sym(zeros(nx,ny));
end


if(any([strfind(char(struct.sym.xdot),'spline'),strfind(char(struct.sym.root),'spline')]))
    struct.splineflag = true;
else
    struct.splineflag = false;
end

% sparse stuff

J = M(r,r);
fprintf('spls | ')
struct.sparseidx = find(J);
I = arrayfun(@(x) find(J(:,x))-1,1:nx,'UniformOutput',false);
struct.rowvals = [];
for ix = 1:nx
    struct.colptrs(ix) = length(struct.rowvals);
    struct.rowvals = [struct.rowvals; I{ix}];
end
struct.colptrs(ix+1) = length(struct.rowvals);

JB = transpose(M(r,r));
struct.sparseidxB = find(JB);
I = arrayfun(@(x) find(JB(:,x))-1,1:nx,'UniformOutput',false);
struct.rowvalsB = [];
for ix = 1:nx
    struct.colptrsB(ix) = length(struct.rowvalsB);
    struct.rowvalsB = [struct.rowvalsB; I{ix}];
end
struct.colptrsB(ix+1) = length(struct.rowvalsB);

fprintf('discontinuities | ')
% discontinuities
ndisc = 0;
for ix = 1:length(struct.sym.xdot)
    tmp_str = char(struct.sym.xdot(ix));
    dl = zeros([length(tmp_str),1]);
    idx_start = [strfind(tmp_str,'heaviside') + 9,strfind(tmp_str,'dirac') + 5]; % find
    if(~isempty(idx_start))
        for iocc = 1:length(idx_start) % loop over occurances
            if(dl(idx_start(iocc))==0) % ensure we do not treat nested discontinuities
                brl = 1; % init bracket level
                str_idx = idx_start(iocc); % init string index
                % this code should be improved at some point
                while brl >= 1 % break as soon as initial bracket is closed
                    str_idx = str_idx + 1;
                    if(strcmp(tmp_str(str_idx),'(')) % increase bracket level
                        brl = brl + 1;
                    end
                    if(strcmp(tmp_str(str_idx),')')) % decrease bracket level
                        brl = brl - 1;
                    end
                end
                idx_end = str_idx;
                arg = tmp_str((idx_start(iocc)+1):(idx_end-1));
                if(ndisc>0)
                    if(strcmp(tmp_str(max((idx_start(iocc)-9),1):(idx_start(iocc)-1)),'heaviside'))
                        if(~any(isequaln(abs(struct.sym.rdisc),abs(sym(arg)))))
                            ndisc = ndisc + 1;
                            struct.sym.rdisc(ndisc) = sym(arg);
                            dl((idx_start(iocc)+1):(idx_end-1))=1;
                        else
                        end
                    elseif(strcmp(tmp_str(max((idx_start(iocc)-5),1):(idx_start(iocc)-1)),'dirac'))
                        if(~any(isequaln(abs(struct.sym.rdisc),abs(sym(arg)))))
                            if(isempty(strfind(arg,','))) % no higher order derivatives
                                ndisc = ndisc + 1;
                                struct.sym.rdisc(ndisc) = sym(arg);
                                dl((idx_start(iocc)+1):(idx_end-1))=1;
                            end
                        else
                        end
                    end
                else
                    if(isempty(strfind(arg,',')))
                        ndisc = ndisc + 1;
                        struct.sym.rdisc(ndisc) = sym(arg);
                        dl((idx_start(iocc)+1):(idx_end-1))=1;
                    end
                end
            end
        end
    end
end

struct.sym.deltadisc = sym(zeros(nx,ndisc));
struct.sym.f = sym(zeros(nx,ndisc));
struct.sym.v = sym(zeros(nx,ndisc));
struct.sym.w = sym(zeros(nx,ndisc));
struct.sym.z = sym(zeros(nx,ndisc));
struct.sym.sdeltadisc = sym(zeros(nx,np,ndisc));
struct.sym.sf = sym(zeros(nx,np,ndisc));
struct.sym.sv = sym(zeros(nx,np,ndisc));
struct.sym.sw = sym(zeros(nx,np,ndisc));
struct.sym.sz = sym(zeros(nx,np,ndisc));


% convert deltas in xdot
syms foo
if(ndisc>0)
    for ix = 1:nx
        summand_ignore = [];
        tmp_str_d = char(struct.sym.xdot(ix));
        idx_start_d = strfind(tmp_str_d,'dirac' ) + 5 ;
        if(~isempty(idx_start_d))
            for iocc_d = 1:length(idx_start_d)
                brl = 1; % init bracket level
                str_idx_d = idx_start_d(iocc_d); % init string index
                % this code should be improved at some point
                while brl >= 1 % break as soon as initial bracket is closed
                    str_idx_d = str_idx_d + 1;
                    if(strcmp(tmp_str_d(str_idx_d),'(')) % increase bracket level
                        brl = brl + 1;
                    end
                    if(strcmp(tmp_str_d(str_idx_d),')')) % decrease bracket level
                        brl = brl - 1;
                    end
                end
                idx_end_d = str_idx_d;
                arg_d = tmp_str_d((idx_start_d(iocc_d)+1):(idx_end_d-1));
                if(strfind(arg_d,',')) % higher order derivatives
                    if(regexp(tmp_str_d((idx_start_d(iocc_d)):(idx_start_d(iocc_d)+2)),'\([1-2],'))
                        arg_d = tmp_str_d((idx_start_d(iocc_d)+3):(idx_end_d-1));
                    elseif(regexp(tmp_str_d((idx_end_d-3):(idx_end_d)),', [1-2])'))
                        arg_d = tmp_str_d((idx_start_d(iocc_d)+1):(idx_end_d-4));
                    end
                end
                str_arg_d = ['dirac(' char(sym(arg_d)) ')' ]; % the char(sym(...)) hopefully ensures consistent ordering of the argument
                str_arg_d1 = ['dirac(1, ' char(sym(arg_d)) ')' ]; % the char(sym(...)) hopefully ensures consistent ordering of the argument
                str_arg_1d = ['dirac(' char(sym(arg_d)) ', 1)' ]; % the char(sym(...)) hopefully ensures consistent ordering of the argument
                str_arg_d2 = ['dirac(2, ' char(sym(arg_d)) ')' ]; % the char(sym(...)) hopefully ensures consistent ordering of the argument
                str_arg_2d = ['dirac(' char(sym(arg_d)) ', 2)' ]; % the char(sym(...)) hopefully ensures consistent ordering of the argument
                dotarg_d = diff(sym(arg_d),'t') + jacobian(sym(arg_d),strsym.xs)*struct.sym.xdot;
                % find the corresponding discontinuity
                for idisc = 1:ndisc;
                    if(isequaln(abs(sym(arg_d)),abs(struct.sym.rdisc(idisc))))
                        if(length(children(struct.sym.xdot(ix)+foo))==2) % if the term is not a sum
                            summands = struct.sym.xdot(ix);
                        else
                            summands = children(struct.sym.xdot(ix));
                        end
                        for is = 1:length(summands)
                            if(isempty(find(is==summand_ignore,1))) % check if we already added that term
                                str_summand = char(summands(is));
                                if(strfind(str_summand,str_arg_d)) % v_i
                                    struct.sym.v(ix,idisc) = struct.sym.v(ix,idisc) + sym(strrep(str_summand,str_arg_d,char(1/abs(dotarg_d))));
                                    summand_ignore = [summand_ignore is];
                                elseif(strfind(str_summand,str_arg_d1)) % w_i
                                    struct.sym.w(ix,idisc) = struct.sym.w(ix,idisc) ...
                                        + sym(strrep(str_summand,str_arg_d1,char(dotarg_d/abs(dotarg_d))));
                                    summand_ignore = [summand_ignore is];
                                elseif(strfind(str_summand,str_arg_1d)) % w_i
                                    struct.sym.w(ix,idisc) = struct.sym.w(ix,idisc) ...
                                        + sym(strrep(str_summand,str_arg_1d,char(dotarg_d/abs(dotarg_d))));
                                    summand_ignore = [summand_ignore is];
                                elseif(strfind(str_summand,str_arg_d2)) % z_i
                                    struct.sym.z(ix,idisc) = struct.sym.z(ix,idisc) ...
                                        + sym(strrep(str_summand,str_arg_d2,char(dotarg_d^2/abs(dotarg_d))));
                                    summand_ignore = [summand_ignore is];
                                elseif(strfind(str_summand,str_arg_2d)) % z_i
                                    struct.sym.z(ix,idisc) = struct.sym.z(ix,idisc) ...
                                        + sym(strrep(str_summand,str_arg_2d,char(dotarg_d^2/abs(dotarg_d))));
                                    summand_ignore = [summand_ignore is];
                                else % f
                                    struct.sym.f(ix,idisc) = struct.sym.f(ix,idisc) + summands(is);
                                    summand_ignore = [summand_ignore is];
                                end
                            end
                        end
                    else
                        struct.sym.f(ix,idisc) = struct.sym.xdot(ix);
                    end
                end
            end
        else
            struct.sym.f(ix,:) = repmat(struct.sym.xdot(ix),[1,ndisc]);
        end
    end
    % compute deltadisc
    for idisc = 1:ndisc
        nonzero_z = any(struct.sym.z(:,idisc)~=0);
        if(nonzero_z)
            M = jacobian(struct.sym.f(:,idisc),strsym.xs)*struct.sym.z(:,idisc) ...
                + struct.sym.w(:,idisc) ...
                - diff(struct.sym.z(:,idisc),'t') - jacobian(struct.sym.z(:,idisc),strsym.xs)*struct.sym.xdot;
        else
            M = struct.sym.w(:,idisc);
        end
        nonzero_M = any(M~=0);
        for ix = 1:nx
            % delta_x_i =
            % sum_j(df_idx_j*(sum_k(df_idx_k*z_k) +
            % w_j - dot{z}_j) + v_i - dot{w}_i +
            % ddot{z}_i
               
            if(struct.sym.v(ix,idisc)~=0)
                struct.sym.deltadisc(ix,idisc) = struct.sym.deltadisc(ix,idisc) ...
                    +struct.sym.v(ix,idisc);
            end
            
            if(struct.sym.w(ix,idisc)~=0)
                struct.sym.deltadisc(ix,idisc) = struct.sym.deltadisc(ix,idisc) ...
                    - diff(struct.sym.w(ix,idisc),'t') - jacobian(struct.sym.w(ix,idisc),strsym.xs)*struct.sym.xdot;
            end
            
            if(struct.sym.z(ix,idisc)~=0)
                struct.sym.deltadisc(ix,idisc) = struct.sym.deltadisc(ix,idisc) ...
                    + diff(struct.sym.z(ix,idisc),'t',2) + diff(jacobian(struct.sym.z(ix,idisc),strsym.xs)*struct.sym.xdot,'t') + ...
                jacobian(jacobian(struct.sym.z(ix,idisc),strsym.xs)*struct.sym.xdot,strsym.xs)*struct.sym.xdot;
            end
            
            if(nonzero_M)
                struct.sym.deltadisc(ix,idisc) = struct.sym.deltadisc(ix,idisc) ...
                    +jacobian(struct.sym.f(ix,idisc),strsym.xs)*M;
            end
            
        end
    end
end

% convert deltas in sx

if(ndisc>0)
    for ix = 1:nx
        for ip = 1:np
            tmp_str_d = char(struct.sym.sx(ix,ip));
            idx_start_d = strfind(tmp_str_d,'dirac' ) + 5 ;
            summand_ignore = [];
            
            if(~isempty(idx_start_d))
                for iocc_d = 1:length(idx_start_d)
                    brl = 1; % init bracket level
                    str_idx_d = idx_start_d(iocc_d); % init string index
                    % this code should be improved at some point
                    while brl >= 1 % break as soon as initial bracket is closed
                        str_idx_d = str_idx_d + 1;
                        if(strcmp(tmp_str_d(str_idx_d),'(')) % increase bracket level
                            brl = brl + 1;
                        end
                        if(strcmp(tmp_str_d(str_idx_d),')')) % decrease bracket level
                            brl = brl - 1;
                        end
                    end
                    idx_end_d = str_idx_d;
                    arg_d = tmp_str_d((idx_start_d(iocc_d)+1):(idx_end_d-1));
                    if(strfind(arg_d,',')) % higher order derivatives
                        if(regexp(tmp_str_d((idx_start_d(iocc_d)):(idx_start_d(iocc_d)+2)),'\([1-2],'))
                            arg_d = tmp_str_d((idx_start_d(iocc_d)+3):(idx_end_d-1));
                        elseif(regexp(tmp_str_d((idx_end_d-3):(idx_end_d)),', [1-2])'))
                            arg_d = tmp_str_d((idx_start_d(iocc_d)+1):(idx_end_d-4));
                        end
                    end
                    str_arg_d = ['dirac(' char(sym(arg_d)) ')' ]; % the char(sym(...)) hopefully ensures consistent ordering of the argument
                    str_arg_d1 = ['dirac(1, ' char(sym(arg_d)) ')' ]; % the char(sym(...)) hopefully ensures consistent ordering of the argument
                    str_arg_1d = ['dirac(' char(sym(arg_d)) ', 1)' ]; % the char(sym(...)) hopefully ensures consistent ordering of the argument
                    str_arg_d2 = ['dirac(2, ' char(sym(arg_d)) ')' ]; % the char(sym(...)) hopefully ensures consistent ordering of the argument
                    str_arg_2d = ['dirac(' char(sym(arg_d)) ', 2)' ]; % the char(sym(...)) hopefully ensures consistent ordering of the argument
                    dotarg_d = diff(sym(arg_d),'t') + jacobian(sym(arg_d),strsym.xs)*struct.sym.xdot;
                    % find the corresponding discontinuity
                    % xdot_i = f_i(x,t) + v_i*delta(t-t0) +
                    % w_i*delta(1)(t-t0) + z_i*delta(t-t0);
                    for idisc = 1:ndisc;
                        if(isequaln(abs(sym(arg_d)),abs(struct.sym.rdisc(idisc))))
                            if(length(children(struct.sym.sx(ix,ip)+foo))==2) % if the term is not a sum
                                summands = struct.sym.sx(ix,ip);
                            else
                                summands = children(struct.sym.sx(ix,ip));
                            end
                            for is = 1:length(summands)
                                if(isempty(find(is==summand_ignore,1))) % check if we already added that term
                                    str_summand = char(summands(is));
                                    if(strfind(str_summand,str_arg_d)) % v_i
                                        struct.sym.sv(ix,ip,idisc) = struct.sym.sv(ix,ip,idisc) + sym(strrep(str_summand,str_arg_d,char(1/abs(dotarg_d))));
                                        summand_ignore = [summand_ignore is];
                                    elseif(strfind(str_summand,str_arg_d1)) % w_i
                                        struct.sym.sw(ix,ip,idisc) = struct.sym.sw(ix,ip,idisc) ...
                                            + sym(strrep(str_summand,str_arg_d1,char(dotarg_d/abs(dotarg_d))));
                                        summand_ignore = [summand_ignore is];
                                    elseif(strfind(str_summand,str_arg_1d)) % w_i
                                        struct.sym.sw(ix,ip,idisc) = struct.sym.sw(ix,ip,idisc) ...
                                            + sym(strrep(str_summand,str_arg_1d,char(dotarg_d/abs(dotarg_d))));
                                        summand_ignore = [summand_ignore is];
                                    elseif(strfind(str_summand,str_arg_d2)) % z_i
                                        struct.sym.sz(ix,ip,idisc) = struct.sym.sz(ix,ip,idisc) ...
                                            + sym(strrep(str_summand,str_arg_d2,char(dotarg_d^2/abs(dotarg_d))));
                                        summand_ignore = [summand_ignore is];
                                    elseif(strfind(str_summand,str_arg_2d)) % z_i
                                        struct.sym.sz(ix,ip,idisc) = struct.sym.sz(ix,ip,idisc) ...
                                            + sym(strrep(str_summand,str_arg_2d,char(dotarg_d^2/abs(dotarg_d))));
                                        summand_ignore = [summand_ignore is];
                                    else % f
                                        struct.sym.sf(ix,ip,idisc) = struct.sym.sf(ix,ip,idisc) + summands(is);
                                        summand_ignore = [summand_ignore is];
                                    end
                                end
                            end
                        else
                            struct.sym.sf(ix,ip,idisc) = struct.sym.sx(ix,ip);
                        end
                    end
                end
            else
                struct.sym.sf(ix,ip,:) = repmat(struct.sym.sx(ix,ip),[1,1,ndisc]);
            end
        end
    end
    % compute sdeltadisc

    nonzero_z = any(struct.sym.z(:,idisc)~=0);
    if(nonzero_z)
        M = jacobian(struct.sym.f(:,idisc),strsym.xs)*struct.sym.z(:,idisc) ...
            + struct.sym.w(:,idisc) ...
            - diff(struct.sym.z(:,idisc),'t') - jacobian(struct.sym.z(:,idisc),strsym.xs)*struct.sym.xdot;
    else
        M = struct.sym.w(:,idisc);
    end
    nonzero_M = any(M~=0);
    for idisc = 1:ndisc
        for ip = 1:np
            nonzero_sz = any(struct.sym.sz(:,ip,idisc)~=0);
            if(nonzero_z)
                if(nonzero_sz)
                    N = jacobian(struct.sym.sf(:,ip,idisc),strsym.xs)*struct.sym.z(:,idisc) ...
                        + jacobian(struct.sym.sf(:,ip,idisc),strsym.sx(:,ip))*struct.sym.sz(:,ip,idisc) ...
                        + struct.sym.sw(:,ip,idisc) ...
                        - diff(struct.sym.sz(:,ip,idisc),'t') - jacobian(struct.sym.sz(:,ip,idisc),strsym.xs)*struct.sym.xdot ...
                        - jacobian(struct.sym.sz(:,ip,idisc),strsym.sx(:,ip))*struct.sym.sx(:,ip);
                else
                    N = jacobian(struct.sym.sf(:,ip,idisc),strsym.xs)*struct.sym.z(:,idisc) ...
                        + struct.sym.sw(:,ip,idisc);
                end
            else
                if(nonzero_sz)
                    N = + jacobian(struct.sym.sf(:,ip,idisc),strsym.sx(:,ip))*struct.sym.sz(:,ip,idisc) ...
                        + struct.sym.sw(:,ip,idisc) ...
                        - diff(struct.sym.sz(:,ip,idisc),'t') - jacobian(struct.sym.sz(:,ip,idisc),strsym.xs)*struct.sym.xdot ...
                        - jacobian(struct.sym.sz(:,ip,idisc),strsym.sx(:,ip))*struct.sym.sx(:,ip);
                else
                    N = struct.sym.sw(:,ip,idisc);
                end
            end
            nonzero_N = any(N~=0);
            for ix = 1:nx
                % delta_x_i =
                % sum_j(df_idx_j*(sum_k(df_idx_k*z_k) +
                % w_j - dot{z}_j) + v_i - dot{w}_i +
                % ddot{z}_i
                if(struct.sym.sz(ix,ip,idisc)~=0)
                    dszdx = jacobian(struct.sym.sz(ix,ip,idisc),strsym.xs);
                    dszdsx = jacobian(struct.sym.sz(ix,ip,idisc),strsym.sx(:,ip));
                    
                    struct.sym.sdeltadisc(ix,ip,idisc) = struct.sym.sdeltadisc(ix,ip,idisc) ...
                        + diff(struct.sym.sz(ix,ip,idisc),'t',2) ...
                        + diff(dszdx*struct.sym.xdot,'t') ...
                        + diff(dszdsx*struct.sym.sx(:,ip),'t') ...
                        + jacobian(dszdx*struct.sym.xdot,strsym.xs)*struct.sym.xdot ...
                        + jacobian(dszdx*struct.sym.xdot,strsym.sx(:,ip))*struct.sym.sx(:,ip) ...
                        + jacobian(dszdsx*struct.sym.sx(:,ip),strsym.xs)*struct.sym.xdot ...
                        + jacobian(dszdsx*struct.sym.sx(:,ip),strsym.sx(:,ip))*struct.sym.sx(:,ip);
                end
                
                if(struct.sym.sv(ix,ip,idisc)~=0)
                    struct.sym.sdeltadisc(ix,ip,idisc) = struct.sym.sdeltadisc(ix,ip,idisc) ...
                        + struct.sym.sv(ix,ip,idisc);
                end
                
                if(struct.sym.sw(ix,ip,idisc)~=0)
                    struct.sym.sdeltadisc(ix,ip,idisc) = struct.sym.sdeltadisc(ix,ip,idisc) ...
                        - diff(struct.sym.sw(ix,ip,idisc),'t') - jacobian(struct.sym.sw(ix,ip,idisc),strsym.xs)*struct.sym.xdot ...
                        - jacobian(struct.sym.sw(ix,ip,idisc),strsym.sx(:,ip))*struct.sym.sx(:,ip);
                        
                end
                
                if(nonzero_M)
                    struct.sym.sdeltadisc(ix,ip,idisc) = struct.sym.sdeltadisc(ix,ip,idisc) ...
                        + jacobian(struct.sym.sf(ix,ip,idisc),strsym.xs)*M;
                end
                if(nonzero_N)
                    struct.sym.sdeltadisc(ix,ip,idisc) = struct.sym.sdeltadisc(ix,ip,idisc) ...
                        + jacobian(struct.sym.sf(ix,ip,idisc),strsym.sx(:,ip))*N;
                end
            end
        end
    end
end

struct.ndisc = ndisc;
if(ndisc>0)
    if(isempty(struct.sym.root));
        struct.sym.root = struct.sym.rdisc;
    else
        struct.sym.root = [struct.sym.root;transpose(struct.sym.rdisc)];
    end
end

fprintf('\r')

end

% better subs
function out = mysubs(in, old, new)
if(~isnumeric(in) && ~isempty(old) && ~isempty(symvar(in)))
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

% old matlabs do not have built-in bandwidth ...
function [ubw,lbw] = cw_bandwidth(M)
[i,j] = find(M);
ubw = max(j-i);
lbw = max(i-j);
end