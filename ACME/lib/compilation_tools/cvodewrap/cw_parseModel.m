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
    fprintf('sx | ')
    struct.sym.sx=struct.sym.J*strsym.sx+struct.sym.dxdotdp;

    fprintf('sx0 | ')
    struct.sym.sx0=jacobian(struct.sym.x0,strsym.ps);

    fprintf('dydx | ')
    struct.sym.dydx=jacobian(struct.sym.y,strsym.xs);
    fprintf('dydp | ')
    struct.sym.dydp=jacobian(struct.sym.y,strsym.ps);
    fprintf('sy | ')
    struct.sym.sy=struct.sym.dydp + struct.sym.dydx*strsym.sx ;
    

    if(nr>0)
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
    if(struct.noadjoint)
        fprintf('xBdot | ')
        struct.sym.xBdot = -transpose(struct.sym.J)*strsym.xBs;
        
        fprintf('int | ')
        struct.sym.int = -transpose(strsym.xBs)*struct.sym.dxdotdp;
        
        fprintf('mu0 | ')
        struct.sym.mu0 = transpose(struct.sym.dydx);
    else
        fprintf('xBdot | ')
        struct.sym.xBdot = sym.zeros(nx,ny);
        
        fprintf('int | ')
        struct.sym.int = sym.zeros(ny,np);
        
        fprintf('mu0 | ')
        struct.sym.mu0 = sym.zeros(nx,ny);
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
        idx_start = strfind(tmp_str,'heaviside') + 9; % find 
        if(~isempty(idx_start))
            for iocc = 1:length(idx_start) % loop over occurances
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
                if(~any(isequaln(struct.sym.rdisc,sym(arg))))
                    ndisc = ndisc + 1;
                    struct.sym.rdisc(ndisc) = sym(arg);
                else
                end
                else
                    ndisc = ndisc + 1;
                    struct.sym.rdisc(ndisc) = sym(arg);
                end
            end
        end
    end
    struct.ndisc = ndisc;
    if(ndisc>0)
        if(isempty(struct.sym.root));
            struct.sym.root = struct.sym.rdisc;
        else
            struct.sym.root = [struct.sym.root;struct.sym.rdisc];
        end
    end
    
    
    fprintf('\r')
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

% old matlabs do not have built-in bandwidth ...
function [ubw,lbw] = cw_bandwidth(M)
[i,j] = find(M);
ubw = max(j-i);
lbw = max(i-j);
end