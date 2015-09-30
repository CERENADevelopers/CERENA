function tw = getTaylorExpansion(w,X,M,order)
if isempty(order)
    order = 2;
end
n_s = length(X);
taylor1 = subs(w,X,M(1:n_s));
taylor2 = 0;
for k=1:n_s
    dw = diff(w,X(k));
    dw = subs(dw,X,M(1:n_s));
    taylor2 = taylor2 + dw * (X(k) - M(k));
end
if order > 1
    taylor3 = 0;
    for k=1:n_s
        for l=1:n_s
            dw = diff(w,X(k));
            d2w = diff(dw,X(l));
            d2w = subs(d2w,X,M(1:n_s));
            taylor3 = taylor3 + 0.5 * d2w * (X(k) - M(k)) * (X(l) - M(l));
        end
    end
end
if order > 2
    taylor4 = 0;
    for k = 1:n_s
        for l = 1:n_s
            for m = 1:n_s
                dw = diff(w,X(k));
                d2w = diff(dw,X(l));
                d3w = diff(d2w,X(m));
                d3w = subs(d3w,X,M(1:n_s));
                taylor4 = taylor4 + 1/6 * d3w * (X(k) - M(k)) * (X(l) - M(l)) * (X(m) - M(m));
            end
        end
    end
end
if order > 3
    taylor5 = 0;
    for k = 1:n_s
        for l = 1:n_s
            for m = 1:n_s
                for n = 1:n_s
                    dw = diff(w,X(k));
                    d2w = diff(dw,X(l));
                    d3w = diff(d2w,X(m));
                    d4w = diff(d3w,X(n));
                    d4w = subs(d4w,X,M(1:n_s));
                    taylor5 = taylor5 + 1/24 * d4w * (X(k) - M(k)) * (X(l) - M(l)) * (X(m) - M(m)) * (X(n) - M(n));
                end
            end
        end
    end
end
if order > 4
    taylor6 = 0;
    for k = 1:n_s
        for l = 1:n_s
            for m = 1:n_s
                for n = 1:n_s
                    for kk=1:n_s
                        dw = diff(w,X(k));
                        d2w = diff(dw,X(l));
                        d3w = diff(d2w,X(m));
                        d4w = diff(d3w,X(n));
                        d5w = diff(d4w,X(kk));
                        d5w = subs(d5w,X,M(1:n_s));
                        taylor6 = taylor6 + 1/120 * d5w * (X(k) - M(k)) * (X(l) - M(l)) * (X(m) - M(m)) * (X(n) - M(n))* (X(kk) - M(kk));
                    end
                end
            end
        end
    end
end
if order > 5
    taylor7 = 0;
    for k = 1:n_s
        for l = 1:n_s
            for m = 1:n_s
                for n = 1:n_s
                    for kk=1:n_s
                        for ll=1:n_s
                            dw = diff(w,X(k));
                            d2w = diff(dw,X(l));
                            d3w = diff(d2w,X(m));
                            d4w = diff(d3w,X(n));
                            d5w = diff(d4w,X(kk));
                            d6w = diff(d5w,X(ll));
                            d6w = subs(d6w,X,M(1:n_s));
                            taylor7 = taylor7 + 1/720 * d6w * (X(k) - M(k)) * (X(l) - M(l)) * (X(m) - M(m)) * (X(n) - M(n))* (X(kk) - M(kk))* (X(ll) - M(ll));
                        end
                    end
                end
            end
        end
    end
end
if order == 1
    tw = taylor1 + taylor2;
elseif order == 2
    tw = taylor1 + taylor2 + taylor3;
elseif order == 3
    tw = taylor1 + taylor2 + taylor3 + taylor4;
elseif order == 4
    tw = taylor1 + taylor2 + taylor3 + taylor4 + taylor5;
elseif order == 5
    tw = taylor1 + taylor2 + taylor3 + taylor4 + taylor5 + taylor6;
elseif order == 6
    tw = taylor1 + taylor2 + taylor3 + taylor4 + taylor5 + taylor6 + taylor7;
end