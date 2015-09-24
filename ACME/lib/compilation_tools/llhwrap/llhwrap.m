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

function llhwrap( modelname,symfun )
    
    warning('off','symbolic:mupadmex:MuPADTextWarning')
    % extract odewrap directory
    if(~ispc)
        llhwrap_path = strrep(which('lw_compileC.m'),'lw_compileC.m','');
    else
        llhwrap_path = strrep(which('lw_compileC.m'),'lw_compileC.m','');
    end
    
    % save current directory
    pwdtmp = pwd;
    
    % change to llhwrap directory to avoid path issues
    eval(['cd ' llhwrap_path])

    % generate modelstruct
    disp('Generating model struct ...')
    modelstruct = eval(symfun);
    
    % do symbolic computations of modelstruct
    disp('Parsing model struct ...')
    modelstruct = lw_parseModel(modelstruct);
    
    % generate C code out of symbolic computations
    disp('Generating C code ...')
    lw_generateC(modelname,modelstruct);
    
    % compile the previously generated C code
    disp('Compiling mex file ...')
    lw_compileC(modelname);
      
    % generate the matlab wrapper
    lw_generateM(modelname,modelstruct);
    
    % change to original directory
    eval(['cd(''' pwdtmp ''')'])
    warning('on','symbolic:mupadmex:MuPADTextWarning')
end

