% cvodewrap generates c mex files for the simulation of systems of ordinary differential equations via CVODES.
% for a detailed description of the required input files and generated simulation files please look at one of
% the example models or the help file cvodewrap.txt.
%
% USAGE:
% ======
% llhwrap( modelname, symfun)
%
% INPUTS:
% =======
% modelname ... specifies the name of the model which will be later used for the naming of the simualation file
% symfun ... specifies a function which returns a struct containing the symbolic definition of the set of differential
%            equations

function llhwrap( varargin )

if(nargin<2)
    error('Must provide modelname and symfun.')
end
modelname = varargin{1};
if(~ischar(modelname))
    error(' modelname must be a string')
end
symfun = varargin{2};
if(~ischar(symfun))
    error(' symfun must be a string')
end
if(exist(symfun,'file') ~= 2)
    error(' symfun must be the name of a matlab function in the matlab path')
end
if nargin > 2
    tdir = varargin{3};
else
    tdir = [];
end

if(isempty(mex.getCompilerConfigurations('fortran')))
    error('No Fortran compiler setup. Please install and configure with MATLAB')
end

if(isempty(mex.getCompilerConfigurations('C')))
    error('No C compiler setup. Please install and configure with MATLAB')
end

warningreset = warning;
warning('off','symbolic:mupadmex:MuPADTextWarning')
warning('off','symbolic:mldivide:RankDeficientSystem')
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
lw_compileC(modelname,modelstruct);

% generate the matlab wrapper
lw_generateM(modelname,modelstruct);
[odewrap_path,~,~]=fileparts(which('lw_compileC.m'));
if(~isempty(tdir))
    %         [odewrap_path,~,~]=fileparts(which('lw_compileC.m'));
    movefile(fullfile(odewrap_path,'models',modelname,['llh_' modelname '.m']),fullfile(tdir,['llh_' modelname '.m']))
    movefile(fullfile(odewrap_path,'models',modelname,[ 'lw_' modelname '.' mexext]),fullfile(tdir,['lw_' modelname '.' mexext]))
    addpath(tdir)
else
    addpath([odewrap_path '/models/' modelname '/'])
end

warning(warningreset)
end

