% vec2mat performs a vector -> matrix conversion.
%   Given a vector x a matrix X, with dimensions n
%   is cunstructed.
%
% Usage:
% ======
% vec2mat(a,n)
%
% x: vector
% n: dimensions of X
% X: matrix
%
% 19/08/09 - Jan Hasenauer

function X = vec2mat(x,n)

% check dimensions
if prod(n) ~= length(x)
    error('Dimension of x and n do not agree.');
end

% reshape vector
X = reshape(x,n);
