% multimax - largest components
%    The function multimin computes the largest
%    component of the matrix X along the dimesions
%    specified by dim.
%
% USAGE:
% ======
% X = multimax(X) maximum of all elements
% X = multimax(X,dim) maximum along the dimensions dim
%
% X: matrix there the largest components are extracted along
% dim: the dimensions along which is minimized
%
% Output Arguments:
% =================
% X = marix of largest components with singleton dimensions removed.
%
% 23/08/09 - Jan Hasenauer

function X = multimax(X,dim)

if nargin < 2 
    dim = 1:length(size(X));
end

for i = 1:length(dim)
    X = max(X,[],dim(i));
end
X = squeeze(X);
