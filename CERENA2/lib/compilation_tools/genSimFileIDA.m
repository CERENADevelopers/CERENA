% genSimFileIDA generates simulation .mex files from the model definition and using the specifed simluation method.
%
% USAGE:
% ======
% genSimFileIDA(modelName,modelDefFileName,simMethod)
%
% INPUTS:
% =======
% modelName ... specifies the name of the model which will be used to name
% the simulation files (string)
% modelDefFileName ... the name of the definition file of the model (string)
% simMethod ... the simulation method to be used. Available options are:
%   'CMEC_X_CC_Y_SC' ... Centered Conditional Moment Equations for moments of order 1
%   to X with sclosure CC and output order Y and scale SC
%                  The options for SC are 
%                  c ... concentration
%                  a ... absolute 
%                  The options for the closure CC are
%                  MF ... Mean field
%                  LD ... Low Dispersion
%                  ZC ... Zero Cumulants
%                  DM ... Derivative Matching

% 
% [10.12.2014] Atefeh Kazeroonian


function System = genSimFileIDA(modelName,modelDefFileName,simMethod)
    System = genmexp(modelName,modelDefFileName,simMethod);
    idawrap(modelName,[simMethod,'_',modelName,'_syms'])
end
