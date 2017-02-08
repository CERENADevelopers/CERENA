% getindstr determines the positions
%    the elements of string S2 have in string S1,
%    and returns this indices as a vector, ind.

function ind = getindstr(S1,S2)

ind = zeros(length(S2),1);
for i = 1:length(S2)
    j = find(strcmp(S1,S2(i)));
    if ~isempty(j)
        ind(i) = j;
    else
        ind(i) = NaN;
    end
end