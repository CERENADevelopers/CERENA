function [fun,exprstr] = sym2fun(expr,varargin)
% fun = sym2fun(expr,argsym1,argsym2,...)
%  generate anonymous function from symbolic expression
%  Both expr and argsym* should be symbolic scalars, vectors, or matrices (max. two-dimensional).
% 2008-09-23, sw

if nargin==1
  error('sym2fun needs at least two input arguments.');
end;

argstr = 'x1';
for k=2:nargin-1
  argstr = [argstr ', x' num2str(k)];
end;

%% COMPLEXE FUNCTIONAL EXPRESSION
if isempty(str2num(char(expr)))

for k=1:nargin-1
    x{k} = sym([]);
    argsym = varargin{k};
    if isempty(argsym)
        x{k} = [];
    end;
    for i=1:size(argsym,1)
        for j=1:size(argsym,2)
            x{k}(i,j) = sym(['x' num2str(k) '_' num2str(i) '_' num2str(j)]);
            if ~isreal(argsym(i,j))
                expr = subs(expr,argsym(i,j),x{k}(i,j));
            end;
        end;
    end;
end;

exprstr = '@(?0) [';
for d1=1:size(expr,1)
	for d2=1:size(expr,2)
		next = char(expr(d1,d2));
		for k=1:nargin-1
                  for i=size(x{k},1):-1:1
			for j=size(x{k},2):-1:1
				next = strrep(next,char(x{k}(i,j)),['?' num2str(k) '(' num2str(i) ',' num2str(j) ')']);
			end;
                  end;
                end;
		exprstr = [exprstr next];
		if d2<size(expr,2)
			exprstr = [exprstr ','];
		elseif d1<size(expr,1)
			exprstr = [exprstr ';'];
		end;
	end;
end;
exprstr = [exprstr '];'];
for k=1:nargin-1
  exprstr = strrep(exprstr,['?' num2str(k)],['x' num2str(k)]);
end;
exprstr = strrep(exprstr,'?0',argstr);

fun = eval(exprstr);

%% ONLY A CONSTANT
else 
    exprstr = char(expr);
    if ~isempty(str2num(char(expr)))
        c = str2num(char(expr));
    end
    exprstr = ['@(' argstr ') [' num2str(c) '];'];
    fun = eval(exprstr);
end



