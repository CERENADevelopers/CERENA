function y=binomial(n,k)
if(k>0.5*n)
    k=n-k;
end
y=n-k+1;
for i=2:k
    y=y*((n-k+i)/i);
end