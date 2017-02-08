function x=getProbInt2(model, t, theta, mu, sigma, n)
if(nargin<6)
    n=8;
end
helper=mod(n,2);
if(helper~=0)
    n=n+1;
end
k=linspace(0.001, 0.5, n/2);
j=logspace(log10(0.501),log10(0.999),n/2);
l=unique([j,k]);
l=sort(l);
y=icdf('logn', l, mu, sigma);
mode=model;
mode.A=model.A(y(1),theta);    
P{1}=simulateFSPext(mode,t)*lognpdf(y(1),mu,sigma);
mode.A=model.A(y(2),theta);
P{2}=simulateFSPext(mode,t)*lognpdf(y(2),mu,sigma);
x=(P{1}+P{2})*0.5*(y(2)-y(1));
helper2=sum(x);
for i=1:size(x,2)
x(:,i)=x(:,i)/helper2(i);
end
for(i=3:length(l))   
P{1}=P{2};
mode.A=model.A(y(i),theta);
P{2}=simulateFSPext(mode,t)*lognpdf(y(i),mu,sigma);
x=x+((P{1}+P{2})*0.5*(y(i)-y(i-1)));
helper2=sum(x);
for o=1:size(x,2)
x(:,o)=x(:,o)/helper2(o);
end
end
