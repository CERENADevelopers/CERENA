%getLikelyhoodMic.m
%
%Calculates the log-likelyhood for  microscopy-measurements, given 
%the measurement data and a FSP-model for the measured 
%system
%
%Usage:
%L=getLikelyhoodMic(model, measurement)
%...
function [LL,W]=getLikelihoodMicMitRauschen(model, measurement, NumTraj)
    function cP = cPr(Ym,x) 
        cP=ones(length(x),1);
        for i=1:length(x(1,:))
            cP=cp.*poisspdf(Ym,x);
        end
    end
%cPr = @(Ym,x) min(Ym == x); 
% Transition probability matrix
dt=measurement.data.time(2)-measurement.data.time(1);
PT = expm(model.A*dt)';
PT = sparse(PT);
for(num=1:NumTraj)
Ym=measurement.data.values{num}';
% Construction of conditional probability
for k = 1:length(measurement.data.time)
% Initialize matrix
W{k} = sparse(zeros(size(model.index,1)));
% Loop: States
for i = 1:size(model.index,1)
% Assign of conditional probability
W{k}(i,i) = cPr(Ym(:,k),model.index(i,:)');
%W{k}(i,i) = cPr(Ym(:,k),model.index(i,:)');
end
end


% Loop: Time points
pj = model.p0'*W{1}
logL = log(sum(pj));
pj = pj/sum(pj);
for k = 2:length(measurement.data.time)
pj = pj*PT*W{k};
logL = logL + log(sum(pj));
pj = pj/sum(pj);
end
LL(:,num)=logL;
end
end
