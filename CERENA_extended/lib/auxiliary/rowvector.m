% rowvector.m converts a vector to a row vector

function vec = rowvector(vec)

if size(vec,2) == 1
    vec = vec.';
end
