% genSimFile generates simulation .mex files from the model definition and using the specifed simluation method.
%
% USAGE:
% ======
% genSimFile(modelName,modelDefFileName,simMethod)
%
% INPUTS:
% =======
% modelName ... specifies the name of the model which will be used to name
% the simulation files (string)
% modelDefFileName ... the name of the definition file of the model (string)
% simMethod ... the simulation method to be used. Available options are:
%   'RRE' ... Reaction Rate Equations. First order System Size Expansion of the mean
%   'EMRE' ... Effective Mesoscopic Rate Equations. Second order System Size Expansion of the mean
%   'LNA' ... Linear Noise Approximation. First order System Size Expansion of the mean and variance
%   'IOS' ... Inverse Omega Square Method. Second order System Size Expansion of the mean and variance
%   'MEC_X_CC_Y_SC' ... Centered Moment Equations for moments of order 1
%   to X with sclosure CC and output order Y and scale SC
%                  The options for SC are 
%                  c ... concentration
%                  a ... absolute 
%                  The options for the closure CC are
%                  MF ... Mean field
%                  LD ... Low Dispersion
%                  ZC ... Zero Cumulants
%                  DM ... Derivative Matching
%   'MEUC_X_CC_Y' ... Unentered Moment Equations for moments of order 1 to X with closure CC and output order Y
%                  The options for the closure CC are
%                  MF ... Mean field
%                  LD ... Low Dispersion
%                  ZC ... Zero Cumulants
%                  DM ... Derivative Matching
% 
% [10.12.2014] Atefeh Kazeroonian


function genSimFileo2(modelName,modelDefFileName,simMethod)
    genmexp(modelName,modelDefFileName,simMethod)
    cvodewrapo2(modelName,[simMethod,'_',modelName,'_syms'])
end


