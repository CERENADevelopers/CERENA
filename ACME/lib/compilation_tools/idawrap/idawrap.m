% idawrap generates c mex files for the simulation of systems of differential algebraic equations via IDAS.
% for a detailed description of the required input files and generated simulation files please look at one of 
% the example models or the help file idawrap.txt.
%
% USAGE:
% ======
% idawrap( modelname, symfun)
% idawrap( modelname, symfun, IC)
% 
% INPUTS:
% =======
% modelname ... specifies the name of the model which will be later used for the naming of the simualation file
% symfun ... specifies a function which returns a struct containing the symbolic definition of the set of differential
%            algebraic equations
% IC   ... User-provided numeric or symbolic inital conditions
%     (default is the prespecified initial conditions in model definition.) 

function idawrap( modelname,symfun )
    
    % extract odewrap directory
    if(~ispc)
        daewrap_path = strrep(which('iw_compileC.m'),'iw_compileC.m','');
    else
        daewrap_path = strrep(which('iw_compileC.m'),'iw_compileC.m','');
    end
    
    
    % save current directory
    pwdtmp = pwd;
    
    % change to odewrap directory to avoid path issues
    eval(['cd ' daewrap_path])
   
    if nargin > 2
        IC = varargin{3};
    else
        IC = [];
    end
    
    % generate modelstruct
    disp('Generating model struct ...')
    if isempty(IC)
        modelstruct = eval(symfun);
    else
        modelstruct = eval([symfun,'(IC)']);
    end
    
    % do symbolic computations of modelstruct
    disp('Parsing model struct ...')
    modelstruct = iw_parseModel(modelstruct);
    
    % generate C code out of symbolic computations
    disp('Generating C code ...')
    iw_generateC(modelname,modelstruct);
    
    % compile the previously generated C code
    disp('Compiling mex file ...')
    iw_compileC(modelname);
      
    % generate the matlab wrapper
    iw_generateM(modelname,modelstruct);
    
    % add the generated files to he MATLAB path
    addpath([daewrap_path 'models/' modelname '/'])
    % change to original directory
    eval(['cd ' pwdtmp])
end

