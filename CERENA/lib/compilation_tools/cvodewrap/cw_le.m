function fun = cw_le(a,b)
% syms x y
% f = symfun(sym('cw_le(x,y)'),[x y]);
% fun = f(a,b);
fun = heaviside(b-a);
end