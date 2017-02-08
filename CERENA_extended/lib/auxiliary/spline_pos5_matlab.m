function yy = spline_pos5_matlab(t,t1,p1,t2,p2,t3,p3,t4,p4,t5,p5,ss,dudt)
    
    tt = [t1,t2,t3,t4,t5];
    xx = log([p1,p2,p3,p4,p5]);
    [b,c,d] = splined2d(tt,xx);
    
    for it = 1:length(t)
        yy(it) = exp(sevald2d(t(it),tt,xx,b,c,d));
    end
    
end

function [b,c,d] = splined2d(x,y)
    n = length(x);
    nm1 = n-1;
    
    d(1) = x(2) - x(1);
    c(2) = (y(2) - y(1)) / d(1);
    for i = 2:nm1
        d(i)   = x(i+1) - x(i);
        b(i)   = 2.0 * (d(i-1) + d(i));
        c(i+1) = (y(i+1) - y(i)) / d(i);
        c(i)   = c(i+1) - c(i);
    end
    
    b(1)   = -d(1);
    b(n) = -d(n-1);
    c(1)   = 0.0;
    c(n) = 0.0;
    c(1)   = c(3) / (x(4) - x(2)) - c(2) / (x(3) - x(1));
    c(n) = c(n-1) / (x(n) - x(n-2)) - c(n-2) / (x(n-1) - x(n-3));
    c(1)   = c(1) * d(1) * d(1) / (x(4) - x(1));
    c(n) = -c(n) * d(n-1) * d(n-1) / (x(n) - x(n-3));
    for i=2:n
        t    = d(i-1) / b(i-1);
        b(i) = b(i) - t * d(i-1);
        c(i) = c(i) - t * c(i-1);
    end
    
    c(n) = c(n) / b(n);
    for ib = 1:nm1
        i    = n - ib ;
        c(i) = (c(i) - d(i) * c(i+1)) / b(i);
    end
    
    b(n) = (y(n) - y(n-1)) / d(n-1) + d(n-1) * (c(n-1) + 2.0 * c(n));
    for i = 1:nm1
        b(i) = (y(i+1) - y(i)) / d(i) - d(i) * (c(i+1) + 2.0 * c(i));
        d(i) = (c(i+1) - c(i)) / d(i);
        c(i) = 3.0 * c(i);
        
    end
    c(n) = 3.0 * c(n);
    d(n) = d(n-1);
end

function yy = sevald2d(t,tt,xx,b,c,d)
    
    n = length(tt);
    i = 1;
    if(tt(i) > t || tt(i+1) < t)
        i = 1;
        j = n;
        while (j > i+1)
            k = floor((i+j)/2);
            if(t<tt(k))
                j = k;
            else
                i = k;
            end
        end
    end
    
    w = t - tt(i);
    yy = xx(i) + w*(b(i) + w*(c(i)+w*d(i)));
    
end