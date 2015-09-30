function I_alpha = getAlphaIndex(I,n)

I_alpha = zeros(size(I,1),n);
for k = 1:n
    I_alpha(:,k) = sum(I == k,2);
end


