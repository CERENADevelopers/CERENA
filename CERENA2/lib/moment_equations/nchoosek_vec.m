function c = nchoosek_vec(n,k)

c = 1;
for i = 1:length(n)
    if n(i) >= k(i)
        c = c*nchoosek(n(i),k(i));
    else
        c = 0;
    end
end
