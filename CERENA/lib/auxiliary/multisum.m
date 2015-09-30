% multisum - sum of components
%    The function multisum computes the sum of the
%    component of the matrix along the dimensions dim.
%
% USAGE:
% ======
% X = multisum(X) sum of all elements
% X = multisum(X,dim) sum along the dimensions dim
%
% X: matrix which is summate up along
% dim: the dimensions along which is minimized
%
% Output Arguments:
% =================
% X = marix of sums with singleton dimensions removed.
%
% 23/08/09 - Jan Hasenauer

function X = multisum(X,dim)

if nargin < 2 
    dim = 1:length(size(X));
end

for i = 1:length(dim)
    X = sum(X,dim(i));
end
X = squeeze(X);

