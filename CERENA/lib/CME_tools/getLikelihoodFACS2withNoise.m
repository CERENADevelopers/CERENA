%getLikelihoodFACS.m
%
%Calculates the log-likelyhood for  FACS-measurements, given 
%the measurement data and a probability distribution for the measured 
%system
%
%
%
%Usage:
%L=getLikelihoodFACS(P,  histogram)
%...
function L=getLikelihoodFACS2withNoise(P, histogram, index)
%% Calculate probability
% This routine calculates the probability Ppdf of measuring a certain state
% while being in any possible state. 
for i=1:length(index)
    Ppdf(:,i)=ones(length(index),1);
    %Ppdf is a column matrix constructed as 
    %Ppdf(i)=(state(i,1)^(state(j,1)/state(j,1)!)*exp(-state(i,1))*
    %(state(i,2)^(state(j,2)/state(j,2)!)*exp(-state(i,2))*...* 
    %(state(i,n)^(state(j,n)/state(j,n)!)*exp(-state(i,n))
    %where n is the size of one state
for q=1:length(index(1,:))
    Ppdf(:,i)=Ppdf(:,i).*poisspdf(index(:,q), repmat(index(i,q),length(index(:,q)),1));
end
end
for j=1:length(histogram.time)
    %calculate the probability of measuring a certain distribution:
    %P(measuring state(i) at time t)=P(state(i) at time t)*
    %P(measuring state(i) while being in state(1))+...+
    %P(state(i) at time t)*P(measuring state(i) while being in state(n))
    P(:,j)=Ppdf*P(:,j);  
    %% prepare data
    for i=1:length(histogram.data.values{1}(:,j))
        if(histogram.data.values{1}(i,j)~=0)
        %calculate log(factorial)
        lh(i,j)=sum(log(1:histogram.data.values{1}(i,j)));
        else
            lh(i,j)=0;
        end
        if (P(i,j)~=0)
            lP(i,j)=log(P(i,j));
        else lP(i,j)=0;
        end
    end
    %% Calculate the Log-Likelihood-function

    L(j)=(sum(log(1:histogram.data.cellsMeasured{j}))-sum(lh(:,j))+sum(lP(:,j).*histogram.data.values{1}(:,j)));
end
       L=sum(L);
  
    
