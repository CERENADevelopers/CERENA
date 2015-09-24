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
function LL=mcmcInterface(par, model, data, measurementMic, measurementFACS, measurementWest)
%function LL=mcmcInterfaceFACS(par,data)
%eval(['model.A=system.A(' theta ');']);
log10(par)
mod=struct('A',model.A(par), 'p0', model.p0, 'index', model.index);
    if(~isempty(measurementWest))
  P=simulateFSPext(mod,measurementWest.data.time);
[P,index]=getMarginalization(P,mod.index, getindstr(model.species,measurementWest.measurands));
 expectation=getExpectationFSP(P, index);
        LLw=getLikelihoodWesternBlotWithNoise(expectation, measurementWest, data.noise, data.trajW);
    else
        LLw=0;
    end

    if(~isempty(measurementMic))
modm=mod;
modm.index=data.indexM;
modm.p0=zeros(length(modm.index),1);
modm.p0(1)=1;
    LLm=getLikelihoodMic7(modm, measurementMic, data.trajM);
    else
        LLm=0;
    end
    if(~isempty(measurementFACS))
            P=simulateFSPext(mod,measurementFACS.time);
[P,index]=getMarginalization(P,mod.index, getindstr(model.species,measurementFACS.data.measurands));
          LLf=getLikelihoodFACS2withNoise(P, measurementFACS, index)

    else
        LLf=0;
    end
    LL=LLm+LLf+LLw;
end
        


