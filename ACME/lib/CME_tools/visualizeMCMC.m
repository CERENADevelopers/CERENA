function p=visualizeMCMC(CHAIN,SSCHAIN, RESULTS)
if(size(CHAIN,2)==2)
scatter((CHAIN(:,1)), 10.^(CHAIN(:,2)), 1, 1:length(CHAIN))
p=figure;
parallelcoords(CHAIN,'labels', RESULTS.names);
else 
if(size(CHAIN,2)==1)
plot((CHAIN),-0.5*SSCHAIN,'.')
xlabel(RESULTS.names{1})
ylabel('loglikelihood')
p=0;
else
plotmatrix((CHAIN))
p=figure;
parallelcoords(CHAIN,'labels', RESULTS.names);
end
end
