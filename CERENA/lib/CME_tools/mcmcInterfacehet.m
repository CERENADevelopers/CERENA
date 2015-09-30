%mcmcInterface.m
%Interface for Likelihood-Functions to the mcmc-Toolbox
%USAGE: LL=mcmcInterface(par,system,model,data)
%should be given to mcmc-Toolbox-Functions as anonymous funtion of
%par and t: @(par, ~) -2*mcmcInterface(...)
%INPUT:
% model ... struct containing necessary information about the model
% 	.A ... Function handle for a matrix of FSP (dp/dt = A p, p(0) = p0).
% 	.p0 ... initial condition of FSP.
% 	.index ... index mapping of FSP.
%   .species ... measureed species
% theta ... String of a rowvector containing the known parameters and 
%           par(i) for the ith unknown parameter e.g, '[5 par(1) 5 par(2)]'  
% measurementFACS ... FACS measurement data
% measurementMic ... Microscopy measurement data
% measurementWest ... Westernblot measurement data
% data ... struct of necessary data
%   .noise ... information about noise in case of westernblot-data
%   .trajM ... number of trajectories in case of Mic-data
%   .trajW ... number of trajectories in case of Westernblot-data
%   .indexM ... microscopy index
%
function LL=mcmcInterfacehet(par, model, data, measurementFACS,n)
%function LL=mcmcInterfaceFACS(par,data)
%eval(['model.A=system.A(' theta ');']);
log10(par)
    if(~isempty(measurementFACS))
            P=getProbInt(model, measurementFACS.time,par(1:(end-2)),par(end-1),par(end), n);
[P,index]=getMarginalization(P,model.index, getindstr(model.species,measurementFACS.data.measurands));
          LL=getLikelihoodFACS2withNoise(P, measurementFACS, index)

    end
end
        


