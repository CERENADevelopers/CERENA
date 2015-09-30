function [p_cMM,cmu_cMM,cC_cMM] = sim_cMM_dm(cMM,t,theta,x0,dx0) 

data.params = theta;

% Determine dimension
n_y = size(cMM.state.stochatic.FSP_index,1);
n_z = length(cMM.state.expectation.state_index);
n_C = size(cMM.state.expectation.C_index,1);
% Set solver options
options = IDASetOptions('RelTol',1e-10,...
'AbsTol',1e-10,...
'VariableTypes',ones(n_y*(1+n_z+n_C),1),...
'suppressAlgVars','on',...
'MaxNumSteps', 10^4,...
'LinearSolver','Dense');
% Initiallize IDAS and solve ODE
IDAInit(@(t,x,dx) cMM_dm(t,x,dx,data),0,x0,dx0,options);
[status,ty,y] = IDASolve(t(2:end),'Normal');
cM = [x0';y'];
% Free memory
IDAFree;

% Reorder result
for iy = 1:n_y
       p_cMM{iy}   = cM(:,iy);
       cmu_cMM{iy} = cM(:,(n_z+n_C)*(iy-1)+n_y    +[1:n_z]);
       cC_cMM{iy}  = cM(:,(n_z+n_C)*(iy-1)+n_y+n_z+[1:n_C]);
end
