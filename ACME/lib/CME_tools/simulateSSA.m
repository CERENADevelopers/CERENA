function X = simulateSSA(v,S,T,x0)

%% SIMULATION
% Initialization
t = T(1);
X = zeros(length(x0),length(T));
X(:,1) = x0;
x = x0;
% Calculation of trajectory
for i = 2:length(T)
    if t >= T(i) 
        % Save time point
        X(:,i) = X(:,i-1);       
    else
        % Simulate till next time point
        while t < T(i)
            % Evaluation of propensities
            w_i = v(x);
            % w_i = SSA_TF(x); 
            cs_w_i = cumsum(w_i);
            % Sampling of reaction time
            if cs_w_i(end) > 0
                t = 1/cs_w_i(end)*log(1/rand) + t;
            else
                t = inf;
            end
            % Check whether state has to be saved
            if t >= T(i)
                % Save time point
                X(:,i) = x;       
            end
            % Sampling of reaction index
            I = max(find([0;cs_w_i(1:end-1)] <= rand*cs_w_i(end)));
            % Update of state
            x = x + S(:,I);
        end
    end
end
