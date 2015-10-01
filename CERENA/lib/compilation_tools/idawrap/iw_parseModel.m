% parseModel parses the model definition struct which contains the symbolic description of the ODE system and
% symbolically derives all expressions required by CVODES
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

function [struct] = iw_parseModel(struct)
    
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
    if(~isfield(struct.sym,'f'))
        error('Model struct is missing the definition of the right hand side f (.sym.f)!')
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
    if(~all([size(struct.sym.x,2)==size(struct.sym.f,2),size(struct.sym.f,2)==size(struct.sym.x0,2)]))
        error('Sizes of x0, xdot and x do not agree!')
    end   
     
    % complete optional fields
    if(~isfield(struct.sym,'u'))
        struct.sym.u = sym.empty(0,0);
    end
    if(~isfield(struct.sym,'k'))
        struct.sym.k = sym.empty(0,0);
    end
    if(~isfield(struct,'t0'))
        struct.t0 = 0;
    end

    syms cj;
    
    nx = length(struct.sym.x);
    nu = length(struct.sym.u);
    nv = length(struct.sym.f);
    np = length(struct.sym.p);
    nk = length(struct.sym.k);
    ny = length(struct.sym.y);
    
    % simplify
    struct.sym.f = simplify(struct.sym.f);
    struct.sym.x0 = simplify(struct.sym.x0);
    struct.sym.dx0 = simplify(struct.sym.dx0);
    struct.sym.M = simplify(struct.sym.M);
    struct.sym.y = simplify(struct.sym.y);
    
    % short strings
    xs = cell(nx,1);
    dxs = cell(nx,1);
    vs = cell(nv,1);
    us = cell(nu,1);
    ps = cell(np,1);
    ks = cell(nk,1);

    
    for j=1:nx
        xs{j} = sprintf('x[%i]',j);
        dxs{j} = sprintf('dx[%i]',j);
    end
    if(nu>0)
        for j=nu
            us{j} = sprintf('u[%i]',j);
        end
    end
    for j=1:nv
        vs{j} = sprintf('v[%i]',j);
    end
    for j=1:np
        ps{j} = sprintf('p[%i]',j);
    end
    for j=1:nk
        ks{j} = sprintf('k[%i]',j);
    end  
    
    % transform into syms
    strsym.xs = sym(xs);
    strsym.dxs = sym(dxs);
    strsym.us = sym(us);
    strsym.vs = sym(vs);
    strsym.ps = sym(ps);
    strsym.ks = sym(ks);
    
    % replace syms by strsyms
    struct.sym.f = mysubs(struct.sym.f,struct.sym.x,strsym.xs);
    struct.sym.y = mysubs(struct.sym.y,struct.sym.x,strsym.xs);
    struct.sym.M = mysubs(struct.sym.M,struct.sym.x,strsym.xs);
    struct.sym.dx0 = mysubs(struct.sym.dx0,struct.sym.x,strsym.xs);

    struct.sym.f = mysubs(struct.sym.f,struct.sym.p,strsym.ps);
    struct.sym.y = mysubs(struct.sym.y,struct.sym.p,strsym.ps);
    struct.sym.x0 = mysubs(struct.sym.x0,struct.sym.p,strsym.ps);
    struct.sym.dx0 = mysubs(struct.sym.dx0,struct.sym.p,strsym.ps);
    struct.sym.M = mysubs(struct.sym.M,struct.sym.p,strsym.ps);
    struct.sym.u = mysubs(struct.sym.u,struct.sym.p,strsym.ps);

    struct.sym.f = mysubs(struct.sym.f,struct.sym.k,strsym.ks);
    struct.sym.y = mysubs(struct.sym.y,struct.sym.k,strsym.ks);
    struct.sym.x0 = mysubs(struct.sym.x0,struct.sym.k,strsym.ks);
    struct.sym.dx0 = mysubs(struct.sym.dx0,struct.sym.k,strsym.ks);
    struct.sym.M = mysubs(struct.sym.M,struct.sym.k,strsym.ks);
    struct.sym.u = mysubs(struct.sym.u,struct.sym.k,strsym.ks);

    struct.sym.f = mysubs(struct.sym.f,struct.sym.u,strsym.us);
    struct.sym.y = mysubs(struct.sym.y,struct.sym.u,strsym.us);
    struct.sym.x0 = mysubs(struct.sym.x0,struct.sym.u,strsym.us);
    struct.sym.dx0 = mysubs(struct.sym.dx0,struct.sym.u,strsym.us);
    struct.sym.M = mysubs(struct.sym.M,struct.sym.u,strsym.us);



    % compute rhs
    if(size(struct.sym.f,2)>size(struct.sym.f,1))
        struct.sym.f = simplify(-transpose(struct.sym.M*strsym.dxs)+struct.sym.f);
    else
        struct.sym.f = simplify(-struct.sym.M*strsym.dxs+struct.sym.f);
    end
    
    % compute derivatives
    struct.sym.dfvdx=simplify(jacobian(struct.sym.f,strsym.xs));
    struct.sym.dfvddx=simplify(jacobian(struct.sym.f,strsym.dxs));
    if(nu>0)
        struct.sym.dfvdu=simplify(jacobian(struct.sym.f,strsym.us));
    else
        struct.sym.dfvdu=sym(ones(nv,0));
    end 
    struct.sym.dfvdp=simplify(jacobian(struct.sym.f,strsym.ps));
    
    % find non-zero elements
    
    fu_nonzero = logical(struct.sym.u ~= 0);
    dvdx_nonzero = logical(struct.sym.dfvdx ~= 0);
    dvddx_nonzero = logical(struct.sym.dfvddx ~= 0);
    if(nu>0)
        dvdu_nonzero = logical(struct.sym.dfvdu ~= 0);
    end
    dvdp_nonzero = logical(struct.sym.dfvdp ~= 0);

    dvdx = cell(nv,nx);
    dvddx = cell(nv,nx);
    dvdu = cell(nv,nu);
    dvdp = cell(nv,np);
    
    for j = 1:nv
        for i = 1:nx
            if(dvdx_nonzero(j,i))
                dvdx{j,i} = sprintf('dvdx[%i]', j + (i-1)*nv);
            else
                dvdx{j,i} = '0';
            end
        end
    end
    
    for j = 1:nv
        for i = 1:nx
            if(dvddx_nonzero(j,i))
                dvddx{j,i} = sprintf('dvddx[%i]', j + (i-1)*nv);
            else
                dvddx{j,i} = '0';
            end
        end
    end
    
    for j = 1:nv
        for i = 1:nu
            if(dvdu_nonzero(j,i))
                dvdu{j,i} = sprintf('dvdu[%i]', j + (i-1)*nv);
            else
                dvdu{j,i} = '0';
            end
        end
    end
    
    for j = 1:nv
        for i = 1:np
            if(dvdp_nonzero(j,i))
                dvdp{j,i} = sprintf('dvdp[%i]', j);
            else
                dvdp{j,i} = '0';
            end
        end
    end
    
    sx = cell(nx,np);
    sdx = cell(nx,np);
    seqsx = cell(nx,np);
    su = cell(nu,np);
    seqsu = cell(nu,np);
    sv = cell(nv,np);
    mus = cell(nx,ny);
    
    for j = 1:nx
        for i = 1:np
            sx{j,i} = sprintf('sx[%i]', j);
            sdx{j,i} = sprintf('sdx[%i]', j);
            seqsx{j,i} = sprintf('sx[%i]', j);
        end
    end

    for j = 1:nu
        for i = 1:np
            if(fu_nonzero(j))
                su{j,i} = sprintf('su[%i]', j);
                seqsu{j,i} = sprintf('sx[%i]', j + (i-1)*nu);
            else
                su{j,i} = '0';
                seqsu{j,i} = '0';
            end
        end
    end
    
    for j = 1:nx
        for i = 1:ny
            mus{j,i} = sprintf('mu[%i]', j + (i-1)*nx);
        end
    end
    
    % transform into syms
    strsym.dvdx = sym(dvdx);
    strsym.dvddx = sym(dvddx);
    strsym.dvdu = sym(dvdu);
    strsym.dvdp = sym(dvdp);
    strsym.sx = sym(sx);
    strsym.sdx = sym(sdx);
    strsym.seqsx = sym(seqsx);
    strsym.su = sym(su);
    strsym.seqsu = sym(seqsu);
    strsym.mus = sym(mus);
    % remove zero values
    
    struct.sym.fsv=strsym.dvdx*strsym.sx+strsym.dvddx*strsym.sdx+strsym.dvdu*strsym.su+strsym.dvdp;
    fsv_nonzero = logical(struct.sym.fsv ~= 0);
    
    % compute symbolics
    struct.sym.fv=struct.sym.f;
    
    struct.sym.fx=strsym.vs;
    struct.sym.fpx0=struct.sym.x0;
    struct.sym.fpdx0=struct.sym.dx0;
    struct.sym.dfxdx= strsym.dvdx;
    struct.sym.dfxddx= strsym.dvddx;
    struct.sym.J = struct.sym.dfxdx + cj*struct.sym.dfxddx;
    
    struct.sym.fsx0=simplify(jacobian(struct.sym.fpx0,strsym.ps));
    struct.sym.fsdx0=simplify(jacobian(struct.sym.fpdx0,strsym.ps));
    
    struct.sym.fu=struct.sym.u;
    if(nu>0)
        struct.sym.dfudp=simplify(jacobian(struct.sym.fu,strsym.ps));
    else
        struct.sym.dfudp=strsym(ones(0,np));
    end
    struct.sym.fy=struct.sym.y;
    struct.sym.dfydx=simplify(jacobian(struct.sym.fy,strsym.xs));
    if(nu>0)
        struct.sym.dfydu=simplify(jacobian(struct.sym.fy,strsym.us));
    else
        struct.sym.dfydu=strsym(ones(ny,0));
    end
    struct.sym.dfydp=simplify(jacobian(struct.sym.fy,strsym.ps));
    if(nu>0)
        struct.sym.fsy=struct.sym.dfydp + struct.sym.dfydx*strsym.seqsx + struct.sym.dfydu*strsym.seqsu ;
    else
        struct.sym.fsy=struct.sym.dfydp + struct.sym.dfydx*strsym.seqsx;
    end
    if(nu>0)
        struct.sym.dfxdp=strsym.dvdp + strsym.dvdu * struct.sym.dfudp;
    else
        struct.sym.dfxdp=strsym.dvdp ;
    end
    
    struct.id = sum(struct.sym.M,2)~=0;
    
    % adjoint sensitivities
    
    struct.sym.mudot = -transpose(struct.sym.dfvdx)*strsym.mus;
    
    struct.sym.int = -transpose(strsym.mus)*struct.sym.dfxdp;
    
    struct.sym.mu0 = transpose(struct.sym.dfydx);
    
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