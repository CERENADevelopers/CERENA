% mat2vec performs a matrix -> vector conversion.
%   Given an n1 x n2 x ... x xq matrix X, this function returns
%   the vector x of length that contains the 
%   columns of the matrix X, stacked below each other.
%
% X: matrix
% x: vector
% n: dimensions of X
%
% 19/08/09 - Jan Hasenauer

function [x,n] = mat2vec(X)

n = size(X);
x = reshape(X,prod(n),1);
