function [p,m,C] = getMofFSPsol(varargin)

P = varargin{1};
index = varargin{2};
y_ind = varargin{3};
y_index = varargin{4};
z_ind = varargin{5};
I_hoMz = varargin{6};
if nargin == 7
    type = varargin{7};
else
    type = 'central';
end
if nargin == 8
    eps = varargin{8};
else
    eps = 10^-35;
end
%% INITIALIZATION
n_y = length(y_ind);
p = {};
m = {};
C = {};

for iy = 1: size(y_index,1)
    % Subset of states
    J = find(sum(bsxfun(@eq,index(:,y_ind),y_index(iy,:)),2) == n_y);
    % Marginal probabilities
    p{iy} = sum(P(:,J),2);
    
    % Conditional mean
    for iz = 1:length(z_ind)
        m{iy}(:,iz) = (P(:,J)*index(J,z_ind(iz)))./(p{iy}(:)+eps);
    end
    
    % Conditional central moments
    for iC = 1:size(I_hoMz,1)
        I = I_hoMz(iC,find(I_hoMz(iC,:)~=0));
        for k = 1:size(P,1)
            C{iy}(k,iC) = P(k,J)*prod(bsxfun(@minus,index(J,z_ind(I)),m{iy}(k,I)),2)/(p{iy}(k)+eps);
        end
    end
end

