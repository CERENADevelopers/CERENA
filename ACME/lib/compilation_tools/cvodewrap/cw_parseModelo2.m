% parseModelo2 parses the model definition struct which contains the symbolic description of the ODE system and
% symbolically derives all expressions required by CVODES for the computation of second order sensitivies
%
% USAGE:
% ======
% [struct] = parseModel(modelstruct)
% 
% INPUTS:
% =======
% modelstruct ... is the struct which contains the symbolic description of the ODE system
% 
% OUTPUS:
% =======
% struct ... is the struct which can be subsequently used by generateC and generateM

function [struct] = cw_parseModelo2(struct)
    
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
        error('Model struct is missing the value for the absolute integration tolerance (.atol)!')
    end
    if(~isfield(struct,'rtol'))
        error('Model struct is missing the value for the relative integration tolerance (.rtol)!')
    end
    if(~isfield(struct,'maxsteps'))
        error('Model struct is missing the maximum number of allowed steps (.maxsteps)!')
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
        struct.sym.root = sym.empty(0,0);
    end
    if(~isfield(struct,'t0'))
        struct.t0 = 0;
    end
    
    if(any(ismember(struct.sym.k,struct.sym.p)))
        error(['Invalid Model: ' char(struct.sym.k(find(ismember(struct.sym.k,struct.sym.p),1))) ' is contained in both p and k!'])
    end
    
    nx = length(struct.sym.x);
    np = length(struct.sym.p);
    nk = length(struct.sym.k);
    ny = length(struct.sym.y);
    nr = length(struct.sym.root);
    
    syms Sx Sdot Sy S0 
    
    Sx = sym(zeros(nx,np));
    
    for j = 1:nx
        for k = 1:np
            eval(['syms S' num2str(j) '_' num2str(k)]);
            eval(['Sx(j,k) = S' num2str(j) '_' num2str(k) ';']);
        end
    end 
    
    % augment system
    fprintf('Sdot | ')
    Sdot = jacobian(struct.sym.xdot,struct.sym.x)*Sx+jacobian(struct.sym.xdot,struct.sym.p);
    
    Sy = sym(zeros(ny,np));
    
    if(size(struct.sym.y,2)>size(struct.sym.y,1))
        struct.sym.y = transpose(struct.sym.y);
    end
    
    % dy/dp = sum(dy/dx*dx/dp) + dy/dp
    for j = 1:nx
        Sy = Sy + repmat(diff(struct.sym.y,struct.sym.x(j)),1,np).*repmat(Sx(j,:),ny,1);
    end
    Sy = Sy + jacobian(struct.sym.y,struct.sym.p);
    
    if(size(struct.sym.x0,2)>size(struct.sym.x0,1))
        struct.sym.x0 = transpose(struct.sym.x0);
    end
    
    S0 = jacobian(struct.sym.x0,struct.sym.p);
    
    struct.sym.x = [struct.sym.x;reshape(Sx,[numel(Sx),1])];
    struct.sym.xdot = [struct.sym.xdot;reshape(Sdot,[numel(Sdot),1])];
    struct.sym.y = [struct.sym.y;reshape(Sy,[numel(Sy),1])];
    struct.sym.x0 = [struct.sym.x0;reshape(S0,[numel(S0),1])];
    
    
    fprintf('reordering | ')
    % optimal reordering (lowest bandwith jacobian)
    J=jacobian(struct.sym.xdot,struct.sym.x);
    M = double(logical(J~=sym(zeros(size(J)))));
    try
        opts.dense = 100;
        opts.agressive = 1;
        r = colamd(M,opts);
    catch
        r = symrcm(M+M');
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
    
    
    nxtrue = nx;
    nytrue = ny;
    
    nx = length(struct.sym.x);
    ny = length(struct.sym.y);
    
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
    
    % event detection
    
    if(nr>0)
       
        dxdp = reshape(strsym.xs(struct.rt((nxtrue+1):end)),nxtrue,np);
        
        % dgdx
        fprintf('dgdx | ')
        struct.sym.drootdx = jacobian(struct.sym.root,strsym.xs(struct.rt(1:nxtrue)));
        
        % ddgdxdx
        fprintf('ddgdxdx | ')
        struct.sym.ddrootdxdx = sym(zeros(nr,nxtrue,nxtrue));
        for ir = 1:nr
            struct.sym.ddrootdxdx(ir,:,:) = hessian(struct.sym.root(ir),strsym.xs(struct.rt(1:nxtrue)));
        end
        
        fprintf('ddgdxdt | ')
        % ddgdxdt = ddgdxdt + ddgdxdx*dxdt
        for ir = 1:nr
            struct.sym.ddrootdxdt(ir,:) = diff(struct.sym.drootdx(ir,:),sym('t')) ...
                + jacobian(jacobian(struct.sym.root,strsym.xs(struct.rt(1:nxtrue)))*struct.sym.xdot(struct.rt(1:nxtrue)),strsym.xs(struct.rt(1:nxtrue))) ... 
                + transpose(squeeze(struct.sym.ddrootdxdx(ir,:,:))*struct.sym.xdot(struct.rt(1:nxtrue)));
        end
        
        fprintf('dgdt | ')
        % dgdt = dgdt + pdgpdx*dxdt
        struct.sym.drootdt = diff(struct.sym.root,sym('t')) ... 
            + struct.sym.drootdx*struct.sym.xdot(struct.rt(1:nxtrue));
        
        fprintf('ddgdtdt | ')
        % ddgdtdt = pddgdtpdt + ddgdxdt*dxdt
        struct.sym.ddrootdtdt = diff(struct.sym.drootdt,sym('t')) ...
            + struct.sym.ddrootdxdt*struct.sym.xdot(struct.rt(1:nxtrue));
        
        fprintf('ddxdtdp | ')
        % ddxdtdp = pddxdtpdp + ddxdtdx*dxdp
        struct.sym.ddxdtdp = jacobian(struct.sym.xdot(struct.rt(1:nxtrue)),strsym.ps) ...
            + jacobian(struct.sym.xdot(struct.rt(1:nxtrue)),strsym.xs(struct.rt(1:nxtrue)))*dxdp;
        
        fprintf('dgdp | ')
        % dgdp = pdgpdp + dgdx*dxdp
        struct.sym.drootdp = jacobian(struct.sym.root,strsym.ps) ...
            + struct.sym.drootdx*dxdp;
        
        fprintf('ddgdpdx | ')
        % ddgdpdx =
        struct.sym.ddrootdpdx = sym(zeros(nr,np,nxtrue));
        for ir = 1:nr
            struct.sym.ddrootdpdx(ir,:,:) = jacobian(struct.sym.drootdp(ir,:),strsym.xs(struct.rt(1:nxtrue)));
        end
        
        fprintf('ddgdpdt | ')
        % ddgdpdt =  pddgdppdt + ddgdpdx*dxdt
        struct.sym.ddrootdpdt = sym(zeros(nr,np,1));
        for ir = 1:nr
            struct.sym.ddrootdpdt(ir,:,:) = diff(struct.sym.drootdp,sym('t')) ...
                + transpose(permute(struct.sym.ddrootdpdx(ir,:,:),[2,3,1])*struct.sym.xdot(struct.rt(1:nxtrue)));
        end
        
        fprintf('ddgdpdt | ')
        % ddgdpdt =  pddgdtpdp + ddgdtdx*dxdp
        struct.sym.ddrootdtdp = sym(zeros(nr,1,np));
        for ir = 1:nr
            struct.sym.ddrootdtdp(ir,:,:) = jacobian(struct.sym.drootdt,strsym.ps) ...
                + struct.sym.ddrootdxdt*dxdp;
        end
        
        fprintf('ddgdpdp | ')
        % dgdx*ddxdpdp
        drootdx_ddxdpdp = sym(zeros(nr,np,np));
        ddxdpdp = reshape(strsym.sx(struct.rt((nxtrue+1):end),:),nxtrue,np,np);
        for ir = 1:nr
            for ix = 1:nxtrue
                drootdx_ddxdpdp(ir,:,:) = drootdx_ddxdpdp(ir,:,:) + struct.sym.drootdx(ir,ix)*ddxdpdp(ix,:,:);
            end
        end
        
        % ddgdpdp = pddgdppdp + ddgdpdx*dxdp + dgdx*dxdp
        struct.sym.ddrootdpdp = sym(zeros(nr,np,np));
        for ir = 1:nr
            struct.sym.ddrootdpdp(ir,:,:) = jacobian(struct.sym.drootdp,strsym.ps) ...
                + permute(struct.sym.ddrootdpdx(ir,:,:),[2,3,1])*dxdp ...
                + squeeze(drootdx_ddxdpdp(ir,:,:));
        end
        
        struct.sym.srootval = struct.sym.drootdp;
        
        struct.sym.s2rootval = struct.sym.ddrootdpdp;

        fprintf('dtaudp | ')
        % dtaudp = -dgdt\dgdp
        struct.sym.sroot=-1/(struct.sym.drootdt)*(struct.sym.drootdp);
        
        
        fprintf('dtaudpdp | ')
        % ddtaudpdp = -dgdt\(dtaudp*ddgdtdt*dtaudp + 2*ddgdpdt*dtaudp
        % + ddgdpdp)
        struct.sym.s2root = sym(zeros(nr,np,np));
        for ir = 1:nr
            struct.sym.s2root(ir,:,:) = -1/(struct.sym.drootdt(ir))*(...
                transpose(struct.sym.sroot(ir,:))*struct.sym.ddrootdtdt(ir)*struct.sym.sroot(ir,:) ...
                + squeeze(struct.sym.ddrootdtdp(ir,:,:))*struct.sym.sroot(ir,:) ...
                + transpose(squeeze(struct.sym.ddrootdtdp(ir,:,:))*struct.sym.sroot(ir,:)) ...
                + squeeze(struct.sym.s2rootval(ir,:,:)));
        end
    else
        struct.sym.sroot = sym.empty(0,np);
        struct.sym.s2root = sym.empty(0,np,np);
        struct.sym.srootval = sym.empty(0,np);
        struct.sym.s2rootval = sym.empty(0,np,np);
    end
    
    % spils solvers
    fprintf('Jv | ')
    struct.sym.Jv = struct.sym.J*strsym.vs;
    
    
    % adjoint sensitivities
    
    fprintf('JvB | ')
    struct.sym.JvB = -transpose(struct.sym.J)*strsym.vs;
    
    % disable this for the time being   
    %struct.sym.xBdot = -transpose(struct.sym.J)*strsym.xBs;
    fprintf('xBdot | ')
    struct.sym.xBdot = sym(zeros(nx,ny));
    
    % disable this for the time being   
    % struct.sym.int = -transpose(strsym.xBs)*struct.sym.dxdotdp;
    fprintf('int | ')
    struct.sym.int = sym(zeros(ny,np));
    
    fprintf('mu0 | ')
    struct.sym.mu0 = transpose(struct.sym.dydx);
    
    if(any([strfind(char(struct.sym.xdot),'spline'),strfind(char(struct.sym.root),'spline')]))
        struct.splineflag = true;
    else
        struct.splineflag = false;
    end
    
    % sparse stuff
    
    fprintf('sparse | ')
    J = M(r,r);
    struct.sparseidx = find(J);
    I = arrayfun(@(x) find(J(:,x))-1,1:nx,'UniformOutput',false);
    struct.rowvals = [];
    for ix = 1:nx
        struct.colptrs(ix) = length(struct.rowvals);
        struct.rowvals = [struct.rowvals; I{ix}];
    end
    struct.colptrs(ix+1) = length(struct.rowvals);
    
    fprintf('sparseB | ')
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