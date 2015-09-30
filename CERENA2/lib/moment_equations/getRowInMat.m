function row = getRowInMat(A,a)

row = ones(size(A,1),1);
for i = 1:size(A,2)
    row = row.*(A(:,i)==a(i));
end
row = find(row);