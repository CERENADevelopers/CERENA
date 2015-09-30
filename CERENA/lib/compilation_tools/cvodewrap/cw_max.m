function fun = cw_max(a, b)
syms x y
f = symfun(sym('cw_max(x,y)'),[x y]);
fun = f(a,b);
% fun = b + heaviside(a-b)*(a-b);
end