% function uncentMom = convertUncent2Cent(I_I,I_alpha,mu,n_z,y_index_str,iy_hoC)
function uncentMom = convertUncent2Cent(varargin)
I_I = varargin{1};
I_alpha = varargin{2};
mu = varargin{3};
n_z = varargin{4};
if nargin>4
    y_index_str = varargin{5};
    iy_hoC = varargin{6};
end

uncentMom = sym(zeros(size(I_I,1),1));
for i = 1:size(I_I,1)
    [K_alpha,K_I] = getLset(I_I(i,:),n_z,'full');
    s = 0;
    for j = 1:size(K_alpha,1)
        if nargin == 4 % moments
            s = s + nchoosek_vec(I_alpha(i,:),K_alpha(j,:)) * prod(mu.^((I_alpha(i,:)-K_alpha(j,:))')) * getC(K_I(j,:));
        elseif nargin == 6 % conditional moments
            s = s + nchoosek_vec(I_alpha(i,:),K_alpha(j,:)) * prod(mu{iy_hoC}.^((I_alpha(i,:)-K_alpha(j,:))')) * getC(K_I(j,:),y_index_str{iy_hoC});
        end
    end
    uncentMom(i) = s;
end

