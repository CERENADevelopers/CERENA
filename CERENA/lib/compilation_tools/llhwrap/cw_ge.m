function fun = cw_ge(a,b)
syms x y
f = symfun(sym('cw_ge(x,y)'),[x y]);
fun = f(a,b);
end