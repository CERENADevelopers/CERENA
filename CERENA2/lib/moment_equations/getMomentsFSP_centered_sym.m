function [M,M_ind] = getMomentsFSP_centered_sym(P,index,order)

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
M = sym(zeros(n_t,n_M));
% Loop: moments
for i = 1:n_M
    I = M_ind(i,find(M_ind(i,:)~=0));
    if length(I) == 1
        % Mean
        M(:,i) =  transpose(P)*prod(index(:,I),2);
    else
        % Central moments with order >= 2
        M(:,i) =  transpose(P)*prod(index(:,I),2) - prod(M(:,I),2);
    end
end
