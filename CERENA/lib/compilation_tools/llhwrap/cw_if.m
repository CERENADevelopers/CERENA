function fun = cw_if(condition, truepart, falsepart)
syms x y z
f = symfun(sym('cw_if(x,y,z)'),[x y z]);
fun = f(condition, truepart, falsepart);
end