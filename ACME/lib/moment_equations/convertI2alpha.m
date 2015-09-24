function alpha = convertI2alpha(I,n)

alpha = zeros(size(I,1),n);
for i = 1:size(I,1)
    for j = 1:n
        alpha(i,j) = sum(I(i,:)==j);
    end    
end


