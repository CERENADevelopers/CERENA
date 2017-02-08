% function [T,X] = simulateSSA_t(v,S,Tsim,x0,options)
function X = simulateSSA_t(v,S,Tsim,x0,options)

%% SIMULATION
% Initialization
nr = size(S,2);
t_ = 0;
x = x0;
T = 0;
X = zeros(length(x0),length(Tsim));
X(:,1) = x0;
% Propensities
Tk = zeros(nr,1);
%ak = v(0,x);
Pk = log(1./rand(nr,1));
counter = 2;
iter = 0;
while t_ < Tsim(end)
%     dtk = (Pk-Tk)./ak;
%     [d,mu] = min(dtk);
%     t = t + d;
%     x = x + S(:,mu);
%     Tk = Tk + ak*d;
%     Pk(mu) = Pk(mu) + log(1/rand);
%     ak = v(0,x);

    op_ode = odeset('RelTol', 1e-10,'AbsTol',1e-10,'events',@(t,A) events(t,A,Pk-Tk));
    [~,~,d,Ak,mu] = ode15s(@(tau,A) v(t_+tau,x),[0,Tsim(end)-t_],zeros(nr,1),op_ode);
    if ~isempty(mu)
        t_ = t_ + d(1);
        % Store
        while ((t_>Tsim(counter)) && (counter<=length(Tsim)))
            X(:,counter) = x;
            counter = counter + 1;
        end
        x = x + S(:,mu(1));
        Tk = Tk + Ak(1,:)';
        Pk(mu(1)) = Pk(mu(1)) + log(1/rand);
    else
        t_ = Tsim(end);
        X(:,counter:end) = repmat(x,1,size(X,2)-counter+1);
    end
    iter = iter + 1;
    % Store
%     if (T(end) == t_)
%        T(end) = t_;
%        X(:,end) = x;
%     else
%         X = [X,x];
%         T = [T,t_];
%     end
    
   
end

end
% % Calculation of trajectory
% while t < Tsim(end)
%     % Evaluation of propensities
%     w_i = v(t,x);
%     % w_i = SSA_TF(x); 
%     cs_w_i = cumsum(w_i);
%     % Sampling of reaction time
%     t = 1/cs_w_i(end)*log(1/rand) + t;
%     % Check whether state has to be saved
%     % Sampling of reaction index
%     I = max(find([0;cs_w_i(1:end-1)] <= rand*cs_w_i(end)));
%     % Update of state
%     x = x + S(:,I);
% end


% % Calculation of trajectory
% for i = 2:length(T)
%     if t >= T(i) 
%         % Save time point
%         X(:,i) = X(:,i-1);       
%     else
%         % Simulate till next time point
%         while t < T(i)
%             % Evaluation of propensities
%             w_i = v(t,x);
%             % w_i = SSA_TF(x); 
%             cs_w_i = cumsum(w_i);
%             % Sampling of reaction time
%             t = 1/cs_w_i(end)*log(1/rand) + t;
%             % Check whether state has to be saved
%             if t >= T(i)
%                 % Save time point
%                 X(:,i) = x;       
%             end
%             % Sampling of reaction index
%             I = max(find([0;cs_w_i(1:end-1)] <= rand*cs_w_i(end)));
%             % Update of state
%             x = x + S(:,I);
%         end
%     end
% end

function [value,isterminal,direction] = events(t,A,A0)
    value = A0 - A;
    isterminal = ones(size(A));
    direction = 0;
end
