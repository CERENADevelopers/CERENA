% setdefault sets the not defined fields of the 
%    structure object to the default values specified
%    in the structure defaultobject.
%
% [object] = setdefault(object,defaultobject)

function defaultobject = cw_setdefault(object,defaultobject)

fieldlist = fieldnames(object);
for i = 1:length(fieldlist)
    defaultobject.(fieldlist{i}) = object.(fieldlist{i});
end

end
