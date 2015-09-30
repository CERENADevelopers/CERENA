function fun = cv_max(a, b)
syms x y
f = symfun(sym('cw_max(x,y)'),[x y]);
fun = f(a,b);
end