function PopDens=getProbInt(model, t, theta, mu, sigma, n)
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
phi_i=icdf('logn', l, mu, sigma);
mode=model;
mode.A=model.A(phi_i(1),theta);    
P{1}=simulateFSPext(mode,t)*lognpdf(phi_i(1),mu,sigma);
mode.A=model.A(phi_i(2),theta);
P{2}=simulateFSPext(mode,t)*lognpdf(phi_i(2),mu,sigma);
PopDens=(P{1}+P{2})*0.5*(phi_i(2)-phi_i(1));

for(i=3:length(l))   
P{1}=P{2};
mode.A=model.A(phi_i(i),theta);
P{2}=simulateFSPext(mode,t)*lognpdf(phi_i(i),mu,sigma);
PopDens=PopDens+((P{1}+P{2})*0.5*(phi_i(i)-phi_i(i-1)));
end
helper2=sum(PopDens);
for o=1:size(PopDens,2)
PopDens(:,o)=PopDens(:,o)/helper2(o);
end
