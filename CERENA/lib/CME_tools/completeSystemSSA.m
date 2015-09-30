function system = completeSystemSSA(system)

for j = 1:system.state.number
    x{j,1} = sprintf('x(%i)',j);
end
for j = 1:length(system.parameter.variable)
    Theta{j,1} = sprintf('Theta(%i)',j);
end
x = sym(x);
Theta = sym(Theta);
v = [];
for j = 1:length(system.reaction)
    v = [v,system.reaction(j).propensity];
end
if isfield(system,'input')
    v = mysubs(v,system.input.variable,system.input.function);
end
v = mysubs(v,system.state.variable,x);
v = mysubs(v,system.parameter.variable,Theta);
v = mysubs(v,system.time,'t');
v = subs(v,'Theta','theta');
str_v = char(v);
str_v = strrep(str_v,',',';');
v = eval(['@(t,x,theta)',str_v(9:end-2)]);

S = system.stoichiometry;

if isfield(system,'output')
    H = transpose(system.output.function);
    if isfield(system,'input')
        H = mysubs(H,system.input.variable,system.input.function);
    end
    H = mysubs(H,system.state.variable,x);
    H = mysubs(H,system.parameter.variable,Theta);
    H = mysubs(H,system.time,'t');
    H = subs(H,'Theta','theta');
    str_H = char(H);
    str_H = strrep(str_H,'x(','x(:,');
    str_H = strrep(str_H,'*','.*');
    str_H = strrep(str_H,'/','./');
    str_H = strrep(str_H,'^','.^');
else
    H = [];
    str_H = '[]';
end
if length(H) <= 1
    H = eval(['@(t,x,theta)',str_H]);
else
    H = eval(['@(t,x,theta)',str_H(11:end-2)]);
end

system.v = v;
system.S = S;
system.H = H;

system.name = 'SSA';
end


% better subs
function out = mysubs(in, old, new)
% if(~isnumeric(in) && ~isempty(old) && ~isempty(findsym(in)))
if(~isnumeric(in) && ~isempty(old) && ~isempty(symvar(in)))
    matVer = ver('MATLAB');
    if(str2double(matVer.Version)>=8.1)
        out = subs(in, old(:), new(:));
    else
        out = subs(in, old(:), new(:), 0);
    end
else
    out = in;
end
end