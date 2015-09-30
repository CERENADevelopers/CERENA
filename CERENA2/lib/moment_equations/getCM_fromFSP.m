function M = getCM_fromFSP(P,index,z_index,z_ind,order)

%% INITIALIZATION
n_z = length(z_ind);
n_z_index = size(z_index,1);
n_y = size(index,2)-n_z;
y_ind = setdiff(1:size(index,2),z_ind);

%% CONSTRUCT MOMENT MATRIX
[E_ind,n_E] = getMomentIndexSet(n_y,order);

%% MOMENT COMPUTATION
M = zeros(size(P,2),n_z_index+n_E*n_z);
% Loop: Stochastic states
for z = 1:n_z_index
    % Subset of states
    J = find(sum(bsxfun(@eq,index(:,z_ind),z_index(z,:)),2) == n_z);
    % State probability
    M(:,z) = sum(P(J,:),1)';
    % Loop: moments
    for i = 1:n_E
        I = y_ind(E_ind(i,find(E_ind(i,:)~=0)));
        M(:,n_z_index+n_E*(z-1)+i) = P(J,:)'*prod(index(J,I),2);
        M(:,n_z_index+n_E*(z-1)+i) = M(:,n_z_index+n_E*(z-1)+i)./(sum(P(J,:),1)'+1e-12);
    end
end
