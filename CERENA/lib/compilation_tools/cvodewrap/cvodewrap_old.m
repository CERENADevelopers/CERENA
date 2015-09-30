% cvodewrap generates c mex files for the simulation of systems of ordinary differential equations via CVODES.
% for a detailed description of the required input files and generated simulation files please look at one of
% the example models or the help file cvodewrap.txt.
%
% USAGE:
% ======
% cvodewrap( modelname, symfun)
% cvodewrap( modelname, symfun, tdir)
% cvodewrap( modelname, symfun, tdir, IC)
%
% INPUTS:
% =======
% modelname ... specifies the name of the model which will be later used
%     for the naming of the simualation file [string]
% symfun ... specifies a function which returns a struct containing the
%     symbolic definition of the set of differential equations. [string]
%     the struct returned by symfun should have the following fields:
%     .sym.x ... vector of states [symbolic]
%     .sym.p ... vector of parameters (for these sensitivities are computed
%     [symbolic]
%     .sym.k (optional) ... vector of constants (for these no sensitivities
%     are computed [symbolic]
%     .sym.xdot ... right-hand-side of the differential equation [symbolic]
%     .sym.x0 ... initial state of the differential equation [symbolic]
%     .sym.y ... output of the differential equation system [symbolic]
%     .sym.root (optional) ... event vector. should contain expressions
%          that have value equal to zero iff the event happens [symbolic]
%     .param (optional ... parametetrisation of the problem
%     .atol (optional) ... default absolute tolerance for simulation
%         (default = 1e-8) [double]
%     .rtol (optional) ... default relative tolerance for simulation
%         (default = 1e-8) [double]
%     .maxsteps (optional) ... default maximum number of steps
%         (default = 1e4) [double]
%     .noadjoint (optional) ... flag that allows the disabling of adjoint
%         sensitivities (default = false) [boolean]
%     .noreorder (optional) ... flag that allows the disabling of state
%         reordering (default = false) [boolean]
% tdir ... target directory where the simulation file should be placed
%     (default is the /models/modelname/ folder in which
% IC   ... User-provided numeric or symbolic inital conditions
%     (default is the prespecified initial conditions in model definition.)


function cvodewrap_old( varargin )

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
if nargin > 3
    IC = varargin{4};
else
    IC = [];
end

if(isempty(mex.getCompilerConfigurations('C')))
    error('No C compiler setup. Please install and configure with MATLAB')
end
if(~isempty(tdir))
    if(exist(tdir,'file') ~= 7)
        error('provided tdir is not a valid path')
    end
end

warningreset = warning;
warning('off','symbolic:mupadmex:MuPADTextWarning')

% generate modelstruct
disp('Generating model struct ...')
if isempty(IC)
    modelstruct = eval(symfun);
else
    modelstruct = eval([symfun,'(IC)']);
end

% do symbolic computations of modelstruct
disp('Parsing model struct ...')
modelstruct = cw_parseModel(modelstruct);

% generate C code out of symbolic computations
disp('Generating C code ...')
cw_generateC(modelname,modelstruct);

% compile the previously generated C code
disp('Compiling mex file ...')
cw_compileC(modelname);

% generate the matlab wrapper
cw_generateM(modelname,modelstruct);

[odewrap_path,~,~]=fileparts(which('cw_compileC.m'));
if(~isempty(tdir))
    %         [odewrap_path,~,~]=fileparts(which('cw_compileC.m'));
    movefile(fullfile(odewrap_path,'models',modelname,['simulate_' modelname '.m']),fullfile(tdir,['simulate_' modelname '.m']))
    movefile(fullfile(odewrap_path,'models',modelname,[ 'cw_' modelname '.' mexext]),fullfile(tdir,['cw_' modelname '.' mexext]))
    addpath(tdir)
else
    addpath([odewrap_path '/models/' modelname '/'])
end
warning(warningreset);
end

