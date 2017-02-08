function [barM,barM_ind] = getMomentsBarMu_sym(y_index,y_ind,z_ind,M_ind_I,p_cMM,cmu_cMM,cC_cMM,moment_order)

% p_cMM{iy}   = cM(:,iy);
% cmu_cMM{iy} = cM(:,(n_z+n_C)*(iy-1)+n_y_index    +[1:n_z]);
% cC_cMM{iy}  = cM(:,(n_z+n_C)*(iy-1)+n_y_index+n_z+[1:n_C]);

%cM
%cMFSP

n_y_index = length(y_index);
n_sy = length(y_ind);
n_sz = length(z_ind);
n_s  = n_sy + n_sz;

M_ind_alpha = getAlphaIndex(M_ind_I,n_sz);

%% GENERATION OF HIGHER-ORDER MOMENTS
[barM_ind,n_M] = getMomentIndexSet(n_s,moment_order);

%% MEAN
muy = zeros(1,n_sy);
muz = zeros(1,n_sz);
% States y
for iy = 1:n_y_index
    muy = muy + (p_cMM(iy)*ones(1,length(y_index(iy,:)))).*y_index(iy,:);
end
% States z
for iy = 1:n_y_index
    muz = muz + (p_cMM(iy)*ones(1,length(cmu_cMM{iy}))).*transpose(cmu_cMM{iy});
end
% Combination
mu = [muy,muz];

%% HIGHER-ORDER MOMENTS
I = barM_ind(n_s+1:end,:);
C = sym(zeros(1,size(I,1)));

% Loop: Moments
for i = 1:(n_M-n_s);
    % Conversion to alpha-index
    I_I  = I(i,find(I(i,:)~=0));
    Iy_I = I_I(find(I_I <= n_sy));
    Iz_I = I_I(find(I_I >  n_sy))-n_sy;
    Iy_alpha = getAlphaIndex(Iy_I,n_sy);
    Iz_alpha = getAlphaIndex(Iz_I,n_sz);
    
    % Loop: Discrete states
    for iy = 1:n_y_index
        % Y components
        CIy = prod((y_index(iy,:) - muy).^(Iy_alpha));
        
        % Z Components
        CIz = zeros(1,1);
        % Generation of combinations
        [K_alpha,K_I] = getLset(Iz_I,n_sz,'full');
        % Loop: combination
        for k = 1:size(K_alpha,1)
            k_I = K_I(k,find(K_I(k,:)~=0));
            k_alpha = K_alpha(k,:);
            if sum(k_alpha) == 0
                CIz = CIz + nchoosek_vec(Iz_alpha,k_alpha)...
                        .*prod((transpose(cmu_cMM{iy})-muz).^(Iz_alpha-k_alpha),2);
            elseif sum(k_alpha) >= 2
                ind = findM(M_ind_alpha,k_alpha);
                if ~isempty(ind)
                    CIz = CIz + nchoosek_vec(Iz_alpha,k_alpha)...
                            .*prod((transpose(cmu_cMM{iy})-muz).^(Iz_alpha-k_alpha),2)...
                            .*cC_cMM{iy}(ind);
                end
            end
        end
        % Moment
        C(:,i) = C(:,i) + CIy.*CIz.*p_cMM(iy); 
    end
end

%% ASIGNMENT Of OUTPUT
barM = [mu,C];