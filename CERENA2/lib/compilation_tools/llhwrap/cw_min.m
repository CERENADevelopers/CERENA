function fun = cv_min(a, b)
syms x y
f = symfun(sym('cw_min(x,y)'),[x y]);
fun = f(a,b);
end