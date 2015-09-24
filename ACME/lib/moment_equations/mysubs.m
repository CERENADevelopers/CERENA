function in=mysubs(in,old,new,repeat,mpnb)

% mysubs   Symbolic substitution.
%   mysubs(S,old,new) replaces old with new in the symbolic expression S.
%   old is a symbolic variable, a string representing a variable name, or
%   a string (quoted) expression. new is a symbolic or numeric variable
%   or expression.
%
%   old and new are vectors or arrays of the same size, each element
%   of old is replaced by the corresponding element of new.  
%
%   when repeat is nonzero, mysubs appies the substitutions as many
%   times as needed until a fixed point is reached
%
% Copyright (C) 2006  Joao Hespanha

% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License along
% with this program; if not, write to the Free Software Foundation, Inc.,
% 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
%
% Modification history:
%
% Created May 19, 2009


if isempty(old)
    return
end

if nargin<4
    repeat=0;
end

if nargin<5
    mpnb=symengine;
end

old=toflatstrlist(old);
new=toflatstrlist(new);
oldin=in;
in=char(sym(in));
if repeat
    while 1
        cmd=['subsex(',char(sym(in)),',zip([',old,'],[',new,'],_equal))'];
        in=evalin(mpnb,cmd);
        if oldin==in
            break;
        end
        oldin=in;
    end
else
    cmd=['subsex(',char(sym(in)),',zip([',old,'],[',new,'],_equal))'];
    in=evalin(mpnb,cmd);
end

if nargin<5
    %delete(mpnb)  % should be used, but apparently create errors
end



return

% old code
str='';
if ischar(old)
    if isnumeric(new)
        str=[old,'=',num2str(new)];
    else
        str=[old,'=',char(new)];
    end
else
    for i=1:prod(size(old))
        if iscell(old)
            o=old{i};
        else
            o=old(i);
        end
        if iscell(new)
            n=new{i};
        else
            n=new(i);
        end
        if isnumeric(n)
            str=[str,char(o),'=',num2str(n),','];
        else
            str=[str,char(o),'=',char(n),','];
        end
    end
    str=str(1:end-1);
end
if ~repeat
    str=['{',str,'}'];
end

% preform substitution
if isa(in,'cell')
    for j=1:numel(in)
        if ischar(in{j})
            in{j}=feval(symengine,'subs',in{j},str);
        else
            for i=1:numel(in{j})
                in{j}(i)=feval(symengine,'subs',char(in{j}(i)),str);
            end
        end
    end
else
    if ischar(in)
        in=feval(symengine,'subs',in,str);
    else
        for i=1:numel(in)
            in(i)=feval(symengine,'subs',char(in(i)),str);
        end
    end
end


function y=toflatstrlist(x)

% no change if already char
if ischar(x)
    y=x;
    return
end

% convert cell numeric to symbolic
if iscell(x) 
    k=find(cellfun(@isnumeric,x));
    if ~isempty(k)
        x{k}=sym(x{k});
    end
end

y=char(sym(x(:).'));
y=strrep(y,']','');
y=strrep(y,'[','');
y=regexprep(y,'^matrix\((.*)\)$','$1');

if nargin<5
    %delete(mpnb)  % should be used, but apparently create errors
end
