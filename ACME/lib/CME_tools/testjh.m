logt = linspace(log10(0.5)-2e-2,log10(0.5)+2e-2,15);
[logt1,logt2] = meshgrid(logt,logt);
LL = zeros(length(logt));

for i = 1:length(logt)
    for j = 1:length(logt)
        LL(i,j) = MODEL.ssfun([logt1(i,j),logt2(i,j)],0);
    end
end

figure;
surf(LL)
