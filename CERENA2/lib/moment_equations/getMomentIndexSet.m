function [M,n_M] = getMomentIndexSet(n_species,order)

% Complete possible set
v = [0:n_species]';
M = v;
for i = 1:order-1
    M = [kron(M,ones(length(v),1)),kron(ones(size(M,1),1),v)];
end
% Sorting and elimination of redundant moments
M = unique(sort(M(2:end,:),2),'rows');
n_M = size(M,1);
