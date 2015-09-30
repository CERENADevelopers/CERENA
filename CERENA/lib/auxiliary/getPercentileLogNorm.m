function x=getPercentileLogNorm(sigma, mu, width, p)
j=1;
P = lognpdf(width,mu,sigma);
while(sum(P(1:j))<p)
    j=j+1;
end
x=width(j);

