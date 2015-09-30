% function M = getMu(I,str)
function M = getMu(varargin)

I = varargin{1};
if nargin == 2
    str = varargin{2};
else
    str = '';
end

n_M = size(I,1);
M = sym(zeros(n_M,1));
for i = 1:n_M
    Ii = I(i,find(I(i,:)~=0));
    if length(Ii) == 1
        M(i) = sym(['mu_' strrep(num2str(Ii,'%d'),' ','') str]);
    else
        error();
    end
end
