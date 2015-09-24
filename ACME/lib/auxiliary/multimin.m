% multimin - smallest components
%    The function multimin computes the smallest
%    component of the matrix X along the dimesions
%    specified by dim.
%
% USAGE:
% ======
% X = multimin(X) minimum of all elements
% X = multimin(X,dim) minimum along the dimensions dim
%
% X: matrix there the smallest components are extracted along
% dim: the dimensions along which is minimized
%
% Output Arguments:
% =================
% X = marix of smallest components with singleton dimensions removed.
%
% 23/08/09 - Jan Hasenauer

function X = multimin(X,dim)

if nargin < 2 
    dim = 1:length(size(X));
end    

for i = 1:length(dim)
    [X,ind] = min(X,[],dim(i));
end
X = squeeze(X);
