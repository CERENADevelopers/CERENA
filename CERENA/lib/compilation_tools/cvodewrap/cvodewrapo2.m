% cvodewrap generates c mex files for the simulation of systems of ordinary differential equations via CVODES.
% for a detailed description of the required input files and generated simulation files please look at one of 
% the example models or the help file cvodewrap.txt. 
%
% USAGE:
% ======
% cvodewrap( modelname, symfun)
% 
% INPUTS:
% =======
% modelname ... specifies the name of the model which will be later used for the naming of the simualation file
% symfun ... specifies a function which returns a struct containing the symbolic definition of the set of differential
%            equations

function cvodewrapo2( modelname,symfun )

    % generate modelstruct
    disp('Generating model struct ...')
    modelstruct = eval(symfun);
    
    % do symbolic computations of modelstruct
    disp('Parsing model struct ...')
    modelstructo2 = cw_parseModelo2(modelstruct);
    modelstruct = cw_parseModel(modelstruct);
    
    % generate C code out of symbolic computiations
    disp('Generating C code ...')
    cw_generateC(modelname,modelstruct);
    cw_generateC([modelname '_o2'],modelstructo2);
    
    % compile the previously generated C code
    disp('Compiling mex file ...')
    cw_compileC(modelname,modelstructo2);
    cw_compileC([modelname '_o2'],modelstructo2);
      
    % generate the matlab wrapper
    cw_generateMo2(modelname,modelstruct,modelstructo2);
    
end


