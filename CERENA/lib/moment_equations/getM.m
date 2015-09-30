function M = getM(I)

n_M = size(I,1);
M = sym(zeros(n_M,1));
for i = 1:n_M
    Ii = I(i,find(I(i,:)~=0));
    if ~isempty(Ii)
%         M(i) = sym(['M_' strrep(num2str(Ii,'%d'),' ','')]);
        M(i) = sym(['M_' strrep(num2str(Ii),'  ','_')]);
    else
        M(i) = 1;
    end
end
