function fun = cw_gt(a,b)
syms x y
f = symfun(sym('cw_gt(x,y)'),[x y]);
fun = f(a,b);
end