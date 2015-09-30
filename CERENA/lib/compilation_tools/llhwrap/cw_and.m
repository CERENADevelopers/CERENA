function fun = cw_and(a,b)
syms x y
f = symfun(sym('cw_and(x,y)'),[x y]);
fun = f(a,b);
end