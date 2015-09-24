function M = getMomentsFSP_uncentered(P,index,order)

%% INITIALIZATION
n_t = size(P,2);
n_s = size(index,2);

%% CONSTRUCT MOMENT MATRIX
% Complete possible set
v = [0:n_s]';
M_ind = v;
for i = 1:order-1
    M_ind = [kron(M_ind,ones(length(v),1)),kron(ones(size(M_ind,1),1),v)];
end
% Sorting and elimination of redundant moments
M_ind = unique(sort(M_ind(2:end,:),2),'rows');
n_M = size(M_ind,1);

%% MOMENT COMPUTATION
M = zeros(n_t,n_M);
% Loop: moments
for i = 1:n_M
    M(:,i) =  P'*prod(index(:,M_ind(i,find(M_ind(i,:)~=0))),2);
end
