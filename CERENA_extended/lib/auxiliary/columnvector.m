% columnvector.m converts a vector to a column vector

function vec = columnvector(vec)

if size(vec,1) == 1
    vec = vec.';
end
