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

function [struct] = lw_parseModel(struct)
    
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
        struct.sym.root = sym.empty(0,1);
    end
    if(~isfield(struct,'t0'))
        struct.t0 = 0;
    end
    
    if(any(ismember(struct.sym.k,struct.sym.p)))
        error(['Invalid Model: ' char(struct.sym.k(find(ismember(struct.sym.k,struct.sym.p),1))) ' is contained in both p and k!'])
    end
    
    % optimal reordering (lowest bandwith jacobian)
    J=simplify(jacobian(struct.sym.xdot,struct.sym.x));
    M = double(logical(J~=sym(zeros(size(J)))));
    try
        %r = hsl_mc73_order(sparse(M+M'),3);
        r = amd(M);
    catch
        r = symrcm(M+M');
    end
    % disable for the time being
    %r = 1:length(struct.sym.xdot);
    [ubw,lbw] = cw_bandwidth(double(M(r,r)));

    struct.ubw = ubw;
    struct.lbw = lbw;
    struct.r = r;
    % inverse permutation
    struct.rt(r) = 1:length(r);
    
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
    xBs = cell(nx,1);

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
        xBs{j} = sprintf('xB[%i]', j-1);
    end

    % transform into syms
    strsym.xs = sym(xs);
    strsym.vs = sym(vs);
    strsym.ps = sym(ps);
    strsym.ks = sym(ks);
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
    struct.sym.J=simplify(jacobian(struct.sym.xdot,strsym.xs));
    struct.sym.JB=-transpose(struct.sym.J);

    struct.sym.dxdotdp=simplify(jacobian(struct.sym.xdot,strsym.ps));

    struct.sym.sx0=simplify(jacobian(struct.sym.x0,strsym.ps));

    struct.sym.dydx=simplify(jacobian(struct.sym.y,strsym.xs));
    struct.sym.dydp=simplify(jacobian(struct.sym.y,strsym.ps));
    
    % events
    struct.sym.drootdx=simplify(jacobian(struct.sym.root,strsym.xs));
    struct.sym.drootdt=simplify(diff(struct.sym.root,sym('t')));
    struct.sym.drootdp=simplify(jacobian(struct.sym.root,strsym.ps));
    
    struct.sym.dtdp = simplify(-1/(struct.sym.drootdx*struct.sym.xdot+struct.sym.drootdt)*(struct.sym.drootdp));
    struct.sym.dtdx = simplify(-1/(struct.sym.drootdx*struct.sym.xdot+struct.sym.drootdt)*(struct.sym.drootdx));
    
    % spils solvers
    struct.sym.Jv = struct.sym.J*strsym.vs;
    struct.sym.JvB = -transpose(struct.sym.J)*strsym.vs;
    
    % adjoint sensitivities
    struct.sym.xBdot = -transpose(struct.sym.J)*strsym.xBs;
    
    struct.sym.int = -transpose(strsym.xBs)*struct.sym.dxdotdp;
    
    struct.sym.mu0 = transpose(struct.sym.dydx);
    
    if(any([strfind(char(struct.sym.xdot),'spline'),strfind(char(struct.sym.root),'spline')]))
        struct.splineflag = true;
    else
        struct.splineflag = false;
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

% old matlabs do not have built-in bandwidth ...
function [ubw,lbw] = cw_bandwidth(M)
[i,j] = find(M);
ubw = max(j-i);
lbw = max(i-j);
end