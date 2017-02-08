% function [L_alpha,L_I] = getLset(I,n_s,option)
function [L_alpha,L_I] = getLset(varargin)

I = varargin{1};
n_s = varargin{2};

if nargin == 2
    option = 'not full';
else
    option = varargin{3};
end

%% COSTRUCTION OF ALPHA-INDEX
for i = 1:n_s
    alpha(i) = sum(I == i);
end
d = cumsum(alpha);

%% COSTRUCTION OF ALPHA-INDEX AND I-INDEX MATRIX
% Initialization
L_alpha = [0:alpha(1)]';
for j = 0:alpha(1)
    L_I(j+1,1:j) = 1*ones(1,j);
end
% Loop: dimensions
for i = 2:length(alpha)
    n_L = size(L_alpha,1);
    if alpha(i) >= 1
        L_alpha = [repmat(L_alpha,alpha(i)+1,1),...
                   reshape(repmat(0:alpha(i),n_L,1),n_L*(alpha(i)+1),1)];
        L_I = [repmat(L_I,alpha(i)+1,1),zeros(n_L*(alpha(i)+1),alpha(i))];
        for j = 1:alpha(i)
            L_I(n_L*j+1:n_L*(j+1),d(i-1)+[1:j]) = i*ones(n_L,j);
        end
    else
        L_alpha = [L_alpha,zeros(n_L,1)];
    end
end
% Delete zero row
if strcmp(option,'not full')
    L_alpha([end],:) = [];
    L_I([end],:) = [];
end
% Reorder
[~,ind] = sort(sum(L_alpha,2));
L_alpha = L_alpha(ind,:);
L_I = sort(L_I(ind,:),2);


