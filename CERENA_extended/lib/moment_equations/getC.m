function M = getC(varargin)

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
    if ~isempty(Ii)
        if length(Ii) == 1
            M(i) = 0;
        else
%             M(i) = sym(['C_' strrep(num2str(Ii,'%d'),' ','') str]);
            M(i) = sym(['C_' strrep(num2str(Ii),'  ','_') str]);
            
        end
    else
        M(i) = 1;
    end
end
