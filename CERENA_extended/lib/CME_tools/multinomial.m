function y=multinomial(n,k)
y=1;
for i=2:length(k)
   y=y*binomial(sum(k(1:i)),sum(k(1:i-1)));
end