% gethistogram - compute histogram
%
% H = gethistorgam(X) bins the elements of X into 10 equally spaced
%    containers and returns the relative amount of elements of X in
%    one bin (discreticed histogram).
%
% H = gethistorgam(X,W) determines a weighted historgam there
%    each element is weighted with W(i).
%
% H = gethistorgam(X,W,G) uses the grid G for binning.
%
% H = gethistorgam(X,W,G,options) can be used to access additional
%    options. In particular on can switch from 'bins' to 'hats'
%
% [H,G] = gethistorgam(...) returns also the gird G.
%
% X: matrix in which each line represents one element.
%    The columms contain different properties according to
%    which binning is performed:
%
%          properties
%          ---------- 
%         |          |
%     X = |          | elements
%         |          |
%          ---------- 
%
% W: vector containing the  weights of the entries of the matrix X.
%    If no 'W' is provided all weights are set to one.
% G: Grid used for bining and hat-function weight computation.
%    formate: G(i).vec = [...], there index i indicates the property
%             of X for which this grid is used.
%
% 27/08/09 - Jan Hasenauer

function [H,G] = gethistogram(varargin)
    
%% CHECK INPUTS AND ASSIGN DEFAULT VALUES
% Data
if nargin >= 1
    X = varargin{1};
    [m,n] = size(X);
else
    error('Not enought inputs.');
end
% Weights
W = ones(m,1);
if nargin >= 2
    W = setdefault(varargin{2},W);
end
% Grid generation
for j = 1:n
    xjmin = min(X(:,j));
    xjmax = max(X(:,j));
    G(j).vec = xjmin:(xjmax-xjmin)/10:xjmax;
end
if nargin >= 3
    G = setdefault(varargin{3},G);
end
% Options
options.type = 'bin'; % type: 'bin' or 'hat'
options.doubleboundary = 1;
options.warning = 'on'; % type: 'off' or 'on'
options.normalize = 0;
if nargin == 4
    options = setdefault(varargin{4},options);
end
nig = 0; % boolean indication whether all x(i,j) are in G

%% ASSIGNMENT OF DIMENSION VECTOR
dim = ones(1,max(n,2));
for j = 1:n
    dim(j) = length(G(j).vec);
end
% Initalization of histogram
if strcmp(options.type,'bin');
    dim(find(dim>1)) = dim(find(dim>1)) - 1;
end
H = zeros(dim);

%% CONSTRUCTION OF VECTOR USED FOR SELECTION OF HISTOGRAM ELEMENTS
% which has to be updated: mp = [1 dim(1) dim(1)*dim(2) dim(1)*dim(2)*dim(3) ...]
mp = ones(1,n);
for j = 2:n
    mp(j) = mp(j-1)*dim(j-1);
end

%% DEFINITION OF UPDATE FUNCTIONS
switch options.type
    case 'bin'
        fpi = @(pi,r,bl,br) 1;
        fni = @(ni,nj,mpj) ni + mpj*(nj-1);
    case 'hat'
        fpi = @(pi,r,bl,br) [(1+bl*options.doubleboundary)*pi*(1-r);(1+br*options.doubleboundary)*pi*r];
        fni = @(ni,nj,mpj) [ni + mpj*(nj-1); ni + mpj*nj];
    otherwise
        error('Unknown option.');
end
        

%% CALCULATION OF HISTOGRAM
for i = 1:m
    ig = 1;
    ni = 1;
    pi = 1;
    for j = 1:n % loop over properties
        v = G(j).vec; % for speed-up
        if ig && (v(1)<=X(i,j)) && (X(i,j)<=v(end)) % element contained in grid?
            nj = sum(X(i,j)>=v(1:end-1));
            pi = fpi(pi,(X(i,j)-v(nj))/(v(nj+1)-v(nj)),nj==1,(nj+1)==dim(j));
            ni = fni(ni,nj,mp(j));
        else  % ni is the element of D which has to be updated (in vector notation)
            ig = 0;
        end
    end
    % updates histogram
    if ig == 1 % within grid?
        H(ni) = H(ni) + W(i)*pi; 
    else
        nig = 1;
    end

end
% Normalization
if options.normalize == 1
    s = sum(mat2vec(H));
    if s > 0
        H = 1/s * H;
    end
end

%% Warning
if (nig == 1) && strcmp(options.warning,'on')
    warning('Not all elements of X are contained in provided grid G!');
end

    