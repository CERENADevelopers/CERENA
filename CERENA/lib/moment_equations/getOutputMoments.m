% function [y,H,My,Iy,n_Iy,M,I] = getOutputMoments(System,I,M,X,XM,uncentM,c,c_kappa,options)
function varargout = getOutputMoments(varargin)
if nargin >= 1
    System = varargin{1};
    n_s = System.state.number;
else
    error('At least one input argument is required!')
end

options.moment_order = 2;
options.moment_order_output = 2;
options.output.calculate = 0;
if nargin >=9
    options = setdefault(varargin{9},options);
end

if nargin >= 2 && ~isempty(varargin{2})
    I = varargin{2};
else
    [I ,~] = getMomentIndexSet(n_s ,options.moment_order);
end
alpha   = convertI2alpha(  I,n_s);
if nargin >= 3 && ~isempty(varargin{3})
    M = varargin{3};
else
    mu_ind = I(find(sum(I>=1,2) == 1),:);
    C_ind  = I(find((2 <= sum(I>=1,2)).*(sum(I>=1,2) <= options.moment_order)),:);
    % Mean
    mu  = getMu(mu_ind);
    % Moments with order >= 2
    C  = getC(C_ind);
    % All moments
    M  = [mu ;C ];
end
n_I = size(I,1);

if nargin >= 4 && ~isempty(varargin{4})
    X = varargin{4};
else
    X = sym('X',[n_s,1]);
end
if nargin >= 5 && ~isempty(varargin{5})
    XM = varargin{5};
else
    XM   = prod(repmat(transpose(X),  n_I,1).^alpha  ,2);
end
if nargin >= 6 && ~isempty(varargin{6})
    uncentM = varargin{6};
else
    mu = M(1:n_s);
    uncentM = convertUncent2Cent(I,alpha,mu,n_s);
end

if nargin >= 7 && ~isempty(varargin{7})
    c = varargin{7};
else
    n_c = length(System.parameter.variable);
    c = sym(zeros(n_c,1));
    for i = 1:n_c
        c(i) = sym(['c' num2str(i,'%d')]);
    end
end
if nargin >= 8 && ~isempty(varargin{8})
    c_kappa = varargin{8};
else
    if isfield(System,'kappa')
        n_kappa = length(System.kappa.variable);
        c_kappa = sym(zeros(n_kappa,1));
        for i = 1:n_kappa
            c_kappa(i) = sym(['c_kappa' num2str(i,'%d')]);
        end
    else
        c_kappa = [];
    end
end

n_sy = System.output.number;
[Iy,n_Iy] = getMomentIndexSet(n_sy,options.moment_order_output);
% Mean
muy = getMu(Iy(1:n_sy,:),'y');
% Moments with order >= 2
Cy = getC(Iy(n_sy+1:end,:),'y');
% All moments
My = [muy;Cy];

%% Output functions
for o = 1:n_sy
    h{o} = System.output.function(o);
    h{o} = mysubs(h{o},System.parameter.variable,c);
    if isfield(System,'kappa')
        h{o} = mysubs(h{o},System.kappa.variable,c_kappa);
    end
    exph{o} = mysubs(h{o},System.state.variable,X);
end

%% Output Map
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
        expH = expand(exph{J});
        H(i) = mysubs(expH,XM(end:-1:1),uncentM(end:-1:1));
        muH(J) = H(i);
        %% MOMENTS WITH ORDER >= 2
    else
        expH = [];
        for k=1:n_sy
            expH = [expH; exph{k} - muH(k)];
        end
        expH = prod(expH(Iyi));
        expH = expand(expH);
        H(i) = mysubs(expH,XM(end:-1:1),uncentM(end:-1:1));
    end
end
H = simplify(H);
%% Calculation of outputs
if isfield(options.output,'calculate') && options.output.calculate == 1
    %% Assembling the output function
    str_y = ['@(x,theta,kappa) ['];
    for i = 1:length(H)-1
        str_y = [str_y,char(H(i)),','];
    end
    str_y = [str_y,char(H(end)),']'];
    for i=length(M):-1:1
        str_y = strrep(str_y,char(M(i)),['x(:,',num2str(i),')']);
    end
    for i=length(c):-1:1
        str_y = strrep(str_y,char(c(i)),['theta(',num2str(i),')']);
    end
    for i=length(c_kappa):-1:1
        str_y = strrep(str_y,char(c_kappa(i)),['kappa(',num2str(i),')']);
    end
    fun_y = eval(str_y);
    %% Assembling the output arguments
    varargout{1} = fun_y;
    varargout{2} = H;
    varargout{3} = My;
    varargout{4} = Iy;
    varargout{5} = n_Iy;
    if nargout > 5
        varargout{6} = M;
        varargout{7} = I;
    end
else
    varargout{1} = H;
    varargout{2} = My;
    varargout{3} = Iy;
    varargout{4} = n_Iy;
    if nargout > 4
        varargout{5} = M;
        varargout{6} = I;
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
