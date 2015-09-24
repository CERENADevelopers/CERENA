function gamma = getDerMatchExponent(m_bar,m)

k = size(m,1);
%gamma = sym('g',[k,1]);
A=zeros(k);
B= zeros(k,1);

for s = 1:k
    ms = m(s,:);
    for p = 1:k
        mp = m(p,:);
        if (all(ge(mp,ms)))
            A(s,p) = prod(factorial(mp)./(factorial(mp-ms).*factorial(ms)));
        end
    end
    if (all(ge(m_bar,ms)))
        B(s,1) = prod(factorial(m_bar)./(factorial(m_bar-ms).*factorial(ms)));
    end
end

gamma = linsolve(A,B);