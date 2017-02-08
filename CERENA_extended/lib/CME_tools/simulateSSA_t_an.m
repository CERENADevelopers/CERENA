function X = simulateSSA_t_an(intv,intvSol,intv_nonlin,S,Tsim,x0)

%% SIMULATION
% Initialization
nr = size(S,2);
t_ = 0;
x = x0;
X = zeros(length(x0),length(Tsim));
X(:,1) = x0;
% Propensities
Tk = zeros(nr,1);
Pk = log(1./rand(nr,1));
counter = 2;
iter = 0;
while t_ < Tsim(end)
    tau = nan(nr,1);
    tau(2:end) = intvSol(t_,x,Pk-Tk);
    tau(isinf(tau)) = NaN;
    d = min(tau);
    if ~isnan(d)
        A_nonlin = intv_nonlin(t_,d,x);
        if A_nonlin>=(Pk(1)-Tk(1))
            tau(1) = fzero(@(tau) (intv_nonlin(t_,tau,x) - (Pk(1)-Tk(1))),0);
        else
            tau(1) = NaN;
        end
    else
        tau(1) = fzero(@(tau) (intv_nonlin(t_,tau,x) - (Pk(1)-Tk(1))),0);
    end
    [d,Ak,mu] = nextReac(intv,t_,tau,x,Tsim(end)-t_);
    
    if ~isempty(mu)
        t_ = t_ + d;
        % Store
        while ((t_>Tsim(counter)) && (counter<=length(Tsim)))
            X(:,counter) = x;
            counter = counter + 1;
        end
        x = x + S(:,mu);
        Tk = Tk + Ak;
        Pk(mu) = Pk(mu) + log(1/rand);
    else
        t_ = Tsim(end);
        X(:,counter:end) = repmat(x,1,size(X,2)-counter+1);
    end
    iter = iter + 1;
end

end

function [d,Ak,mu] = nextReac(intv,t,tau,x,Tend)
[d,mu] = min(tau);
if (d > Tend)
    d = [];
    mu = [];
    Ak = [];
else
    Ak = intv(t,d,x);
end
end