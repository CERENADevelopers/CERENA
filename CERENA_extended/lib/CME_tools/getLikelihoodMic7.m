%getLikelyhoodMic.m
%
%Calculates the log-likelyhood for  microscopy-measurements, given 
%the measurement data and a FSP-model for the measured 
%system
%
%Usage:
%L=getLikelyhoodMic(model, measurement)
%...
function [LL,b]=getLikelihoodMic7(model, measurement, NumTraj)
%cPr = @(Ym) 1/(2*pi)*exp(-0.5*sum(bsxfun(@minus,Y',Ym).^2,1));A=sparse(model.A);
%options=odeset('Jacobian', A, 'Vectorized', 'on');
Y = model.index;
data=model.A;
function cP = cPr(Ym,x) 
        for j=1:length(Ym)
                cP(:,j)=ones(size(x,1),1)*Ym(j);
                cP(:,j)=poisspdf(cP(:,j),x(:,j));
        end
        cP=(prod(cP,2));
    end 
%cPr = @(Ym) 1/(2*pi)*exp(-0.5*sum(bsxfun(@minus,Y',Ym).^2,1));
%cPr = @(Ym,x) min(Ym == x); 
% Transition probability matrix
dt=measurement.data.time(2)-measurement.data.time(1);
for(num=1:NumTraj)
Ym=measurement.data.values{num}';

% Loop: Time points
pj = model.p0.*cPr(Ym(:,1),Y);
logL = log(sum(pj));
pj = pj/sum(pj);
    options = CVodeSetOptions('JacobianFn', @Jac,...
                          'RelTol',1.e-7,...
                          'AbsTol',1.e-7,...
                          'LinearSolver','GMRES',...
                          'UserData',data,...
                          'MaxNumSteps', 100000);
 CVodeInit(@rhsfn, 'BDF', 'Newton', 0, pj, options);                     
for k = 2:length(measurement.data.time)
CVodeReInit(0, pj, options);
[flag, dummyt, b]=CVode(dt, 'Normal');
for o=1:length(b)
if b(o)<0
Pj(o)=0;
else
    Pj(o)=b(o);
end
end
    pj = (Pj'.*cPr(Ym(:,k),Y));
logL = logL + log(sum(pj));
pj = pj/sum(pj);

end
CVodeFree;
LL(:,num)=logL;
end
LL=sum(LL);
    function [J,flag,new_data]=Jac(~,~,~,y,data)
    J=data*y;
    flag=0;
    new_data=[];
    end
    function [yd, flag, new_data]= rhsfn(~,y,data)
        yd=data*y;
        flag=0;
        new_data=[];
    end
end
