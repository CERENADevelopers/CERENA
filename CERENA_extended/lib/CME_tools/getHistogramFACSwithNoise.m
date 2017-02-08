%getHistogramFACS.m
%Produces a histogram out of FACS-data given the measurement and an index
%mapping between the states and the bins of the histogram
%
%Usage:
%histogram=getHistogramFACS(measurement, index, P)
%
%Input: 
% measurement ... information about measurement:
%   .time ... vector of time instances then measurements are
%           performed.
%   .data ... sturcture containing the experimental data:
%       .measurands ... names of measurands.
%       .values{i} ... measurement data at time instance i:
%
% index ... index mapping between system state and histogram bins
%           
%
%Output:
% histogram ... information about measurement:
%   .time ... vector of time instances then measurements are
%           performed.
%   .data ... sturcture containing the histogram data:
%       .values{i} ... measurement data at time instance i:
%       .measurands
%

function h=getHistogramFACS(measurement, bins)
h.data.values=zeros(length(bins),length(measurement.time));
h.time=measurement.time;
h.data.measurands=measurement.data.measurands;
for j=1:length(h.time)
    h.data.cellsMeasured(j)=length(measurement.data.values{j});
    for k=1:h.data.cellsMeasured(j)
        mes=abs(repmat(measurement.data.values{j}(k,:), length(bins), 1)-bins);
        y=0;
        a=0;
        while(a~=1)
            y=y+1;
            a=(mes(y,:)==min(mes));
        end
        ind(k,j)=y;
        h.data.values(ind(k,j),j)=h.data.values(ind(k,j),j)+1;
    end
end