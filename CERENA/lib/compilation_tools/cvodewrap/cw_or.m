function fun = cw_or(a,b)
% syms x y
% f = symfun(sym('cw_or(x,y)'),[x y]);
% fun = f(a,b);
fun = 1-cw_and(a,b);
end