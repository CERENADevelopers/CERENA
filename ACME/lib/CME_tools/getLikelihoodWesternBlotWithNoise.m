%getLilelihoodWesterBlot.m
%Computes the loglikelihood for Westernblot-data
%
%USAGE
%======
%L=getLikelihoodWesternBlot(expectation, measurement, noise, numberMeas)
%
%INPUTS
%======
%expectation ... expectation for probability distribution
%measurement ... Westernblot-data
%noise ... information about the noise which corrupted the measurement
%   .sigma ... standard deviation
%   .mu ... mean
%numberMeas ... number of different experiments
%...
function L=getLikelihoodWesternBlotWithNoise(expectation, measurement, noise, numberMeas)
n=length(measurement.data.time);
for i=1:numberMeas
L(:,i)=-0.5*n*log(2*pi)+sum(-0.5*n*log(noise.sigma.^2)-(1./(2*noise.sigma.^2)).*...
    (sum(abs(measurement.data.values{i}-expectation),2)-noise.mu).^2);
    
end
L=sum(L);
end