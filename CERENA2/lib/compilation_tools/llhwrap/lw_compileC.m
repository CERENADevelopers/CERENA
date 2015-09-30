% compileC compiles the mex simulation file. 
%
% USAGE:
% ======
% compileC( modelname)
% 
% INPUTS:
% =======
% modelname ... specifies the name of the model which will be later used for the naming of the simualation file. this
%               should be the same which was used when calling either cvodewrap or generateC.
% modelstruct ... is the struct which contains the symbolic description of the ODE system

function lw_compileC(filename,modelstruct)

if isempty(modelstruct)
    modelstruct.coptim = '-O3';
    modelstruct.debug = false;
else
    if(~isfield(modelstruct,'coptim'))
        modelstruct.coptim = '-O3';
    end
    if(~isfield(modelstruct,'debug'))
        modelstruct.debug = false;
    end
end
    
[odewrap_path,~,~] = fileparts(which('lw_compileC.m'));

sundials_path = fullfile(odewrap_path,'sundials-2.6.1');
sundials_ver = '2.6.1';

ssparse_path = fullfile(odewrap_path,'SuiteSparse');
ssparse_ver = '4.4.4';

lapack_path = fullfile(odewrap_path,'lapack-3.5.0'); % currently not used, lapack implementation still needs to be done
lapack_ver = '3.5.0';

includesstr = '';
includesstr = strcat(includesstr,' -I"', fullfile(sundials_path, 'include'), '"');
includesstr = strcat(includesstr,' -I"', fullfile(sundials_path, 'src','cvodes'), '"');
includesstr = strcat(includesstr,' -I"', fullfile(odewrap_path, 'models', filename ), '"');
includesstr = strcat(includesstr,' -I"', fullfile(odewrap_path), '"');
includesstr = strcat(includesstr,' -I"', fullfile(ssparse_path, 'KLU','Include'), '"');
includesstr = strcat(includesstr,' -I"', fullfile(ssparse_path, 'AMD','Include'), '"');
includesstr = strcat(includesstr,' -I"', fullfile(ssparse_path, 'COLAMD','Include'), '"');
includesstr = strcat(includesstr,' -I"', fullfile(ssparse_path, 'BTF','Include'), '"');
includesstr = strcat(includesstr,' -I"', fullfile(ssparse_path, 'SuiteSparse_config'), '"');

% compile directory
if(~exist(fullfile(odewrap_path,'models',mexext), 'dir'))
    mkdir(fullfile(odewrap_path,'models',mexext))
end

% read version number from versions.txt and decide whether object have to
% be regenerated
if(~exist(fullfile(odewrap_path,'models',mexext,'versions.txt'),'file'))
    del_sundials = true;
    del_ssparse = true;
    del_lapack = true;
    fid = fopen(fullfile(odewrap_path,'models',mexext,'versions.txt'),'w');
    fprintf(fid,[sundials_ver '\r']);
    fprintf(fid,[ssparse_ver '\r']);
    fprintf(fid,[lapack_ver '\r']);
    fclose(fid);
else
    fid = fopen(fullfile(odewrap_path,'models',mexext,'versions.txt'),'r');
    sundials_objver = fgetl(fid);
    ssparse_objver = fgetl(fid);
    lapack_objver = fgetl(fid);
    fclose(fid);
    del_sundials = isnewer(sundials_ver,sundials_objver);
    if(del_sundials)
        display('Newer version of Sundials! Regenerating object files.')
    end
    del_ssparse = isnewer(ssparse_ver,ssparse_objver);
    if(del_ssparse)
        display('Newer version of SuiteSparse! Regenerating object files.')
    end
    del_lapack = isnewer(lapack_ver,lapack_objver);
    if(del_lapack)
        display('Newer version of Lapack! Regenerating object files.')
    end
    fid = fopen(fullfile(odewrap_path,'models',mexext,'versions.txt'),'w');
    fprintf(fid,[sundials_ver '\r']);
    fprintf(fid,[ssparse_ver '\r']);
    fprintf(fid,[lapack_ver '\r']);
    fclose(fid);
end

sources_sundials = {
    %    'src/cvodes/cvodes_lapack.c';
    fullfile('src','cvodes','cvodes_band.c');
    fullfile('src','cvodes','cvodes_bandpre.c');
    fullfile('src','cvodes','cvodes_bbdpre.c');
    fullfile('src','cvodes','cvodes_direct.c');
    fullfile('src','cvodes','cvodes_dense.c');
    fullfile('src','cvodes','cvodes_sparse.c');
    fullfile('src','cvodes','cvodes_diag.c');
    fullfile('src','cvodes','cvodea.c');
    fullfile('src','cvodes','cvodes.c');
    fullfile('src','cvodes','cvodes_io.c');
    fullfile('src','cvodes','cvodea_io.c');
    fullfile('src','cvodes','cvodes_spils.c');
    fullfile('src','cvodes','cvodes_spbcgs.c');
    fullfile('src','cvodes','cvodes_spgmr.c');
    fullfile('src','cvodes','cvodes_sptfqmr.c');
    fullfile('src','cvodes','cvodes_klu.c');
    fullfile('src','sundials','sundials_band.c');
    fullfile('src','sundials','sundials_dense.c');
    fullfile('src','sundials','sundials_sparse.c');
    fullfile('src','sundials','sundials_iterative.c');
    fullfile('src','sundials','sundials_nvector.c');
    fullfile('src','sundials','sundials_direct.c');
    fullfile('src','sundials','sundials_spbcgs.c');
    fullfile('src','sundials','sundials_spgmr.c');
    fullfile('src','sundials','sundials_sptfqmr.c');
    fullfile('src','sundials','sundials_math.c');
    fullfile('src','nvec_ser','nvector_serial.c');
    };
sourcesstr_sundials = '';
for j=1:length(sources_sundials)
    sourcesstr_sundials = strcat(sourcesstr_sundials,' "',fullfile(sundials_path,sources_sundials{j}),'"');
end

sources_ssparse = {
    fullfile('KLU','Source','klu_analyze_given.c');
    fullfile('KLU','Source','klu_analyze.c');
    fullfile('KLU','Source','klu_defaults.c');
    fullfile('KLU','Source','klu_diagnostics.c');
    fullfile('KLU','Source','klu_dump.c');
    fullfile('KLU','Source','klu_extract.c');
    fullfile('KLU','Source','klu_factor.c');
    fullfile('KLU','Source','klu_free_numeric.c');
    fullfile('KLU','Source','klu_free_symbolic.c');
    fullfile('KLU','Source','klu_kernel.c');
    fullfile('KLU','Source','klu_memory.c');
    fullfile('KLU','Source','klu_refactor.c');
    fullfile('KLU','Source','klu_scale.c');
    fullfile('KLU','Source','klu_sort.c');
    fullfile('KLU','Source','klu_solve.c');
    fullfile('KLU','Source','klu_tsolve.c');
    fullfile('KLU','Source','klu.c');
    fullfile('AMD','Source','amd_1.c');
    fullfile('AMD','Source','amd_2.c');
    fullfile('AMD','Source','amd_aat.c');
    fullfile('AMD','Source','amd_control.c');
    fullfile('AMD','Source','amd_defaults.c');
    fullfile('AMD','Source','amd_dump.c');
    fullfile('AMD','Source','amd_global.c');
    fullfile('AMD','Source','amd_info.c');
    fullfile('AMD','Source','amd_order.c');
    fullfile('AMD','Source','amd_post_tree.c');
    fullfile('AMD','Source','amd_postorder.c');
    fullfile('AMD','Source','amd_preprocess.c');
    fullfile('AMD','Source','amd_valid.c');
    %fullfile('AMD','Source','amd.f');
    %fullfile('AMD','Source','amdbar.f');
    fullfile('COLAMD','Source','colamd_global.c');
    fullfile('COLAMD','Source','colamd.c');
    fullfile('BTF','Source','btf_maxtrans.c');
    fullfile('BTF','Source','btf_order.c');
    fullfile('BTF','Source','btf_strongcomp.c');
    fullfile('SuiteSparse_config','SuiteSparse_config.c');
    };
sourcesstr_ssparse = '';
for j=1:length(sources_ssparse)
    sourcesstr_ssparse = strcat(sourcesstr_ssparse,' "', fullfile(ssparse_path,sources_ssparse{j}),'"');
end


% objects
objects_sundials = {
%    'cvodes_lapack.o';
    'cvodes_band.o';
    'cvodes_bandpre.o';
    'cvodes_bbdpre.o';
    'cvodes_direct.o';
    'cvodes_dense.o';
    'cvodes_sparse.o';
    'cvodes_diag.o';
    'cvodea.o';
    'cvodes.o';
    'cvodes_io.o';
    'cvodea_io.o';
    'cvodes_spils.o';
    'cvodes_spbcgs.o';
    'cvodes_spgmr.o';
    'cvodes_sptfqmr.o';
    'cvodes_klu.o';
    'sundials_band.o';
    'sundials_dense.o';
    'sundials_sparse.o';
    'sundials_iterative.o';
    'sundials_nvector.o';
    'sundials_direct.o';
    'sundials_spbcgs.o';
    'sundials_spgmr.o';
    'sundials_sptfqmr.o';
    'sundials_math.o';
    'nvector_serial.o';
    };
if(ispc)
    objects_sundials = strrep(objects_sundials, '.o', '.obj');
end

objects_ssparse = {
    'klu_analyze_given.o';
    'klu_analyze.o';
    'klu_defaults.o';
    'klu_diagnostics.o';
    'klu_dump.o';
    'klu_extract.o';
    'klu_factor.o';
    'klu_free_numeric.o';
    'klu_free_symbolic.o';
    'klu_kernel.o';
    'klu_memory.o';
    'klu_refactor.o';
    'klu_scale.o';
    'klu_sort.o';
    'klu_solve.o';
    'klu_tsolve.o';
    'klu.o';
    'amd_1.o';
    'amd_2.o';
    'amd_aat.o';
    'amd_control.o';
    'amd_defaults.o';
    'amd_dump.o';
    'amd_global.o';
    'amd_info.o';
    'amd_order.o';
    'amd_post_tree.o';
    'amd_postorder.o';
    'amd_preprocess.o';
    'amd_valid.o';
    %'amd.o';
    %'amdbar.o';
    'colamd_global.o';
    'colamd.o';
    'btf_maxtrans.o';
    'btf_order.o';
    'btf_strongcomp.o'
    'SuiteSparse_config.o';
    };
if(ispc)
    objects_ssparse = strrep(objects_ssparse, '.o', '.obj');
end

objectsstr = '';
for j=1:length(objects_sundials)
    objectsstr = strcat(objectsstr,' "',fullfile(odewrap_path,'models',mexext,objects_sundials{j}),'"');
end

for j=1:length(objects_ssparse)
    objectsstr = strcat(objectsstr,' "',fullfile(odewrap_path,'models',mexext,objects_ssparse{j}),'"');
end

for j=1:length(sources_sundials)
    if(~exist(fullfile(odewrap_path,'models',mexext,objects_sundials{j}), 'file') || del_sundials)
        eval(['mex COPTIMFLAGS=''-O3 -DNDEBUG'' -c -outdir ' fullfile(odewrap_path,'models',mexext) includesstr ' ' fullfile(sundials_path,sources_sundials{j})]);
    end
end

for j=1:length(sources_ssparse)
    if(~exist(fullfile(odewrap_path,'models',mexext,objects_ssparse{j}), 'file') || del_ssparse)
        eval(['mex COPTIMFLAGS=''-O3 -DNDEBUG'' -c -outdir ' fullfile(odewrap_path,'models',mexext) includesstr ' ' fullfile(ssparse_path,sources_ssparse{j})]);
    end
end

if(isunix)
    if(~ismac)
        CLIBS = 'CLIBS=''\$CLIBS -lrt''';
    else
        CLIBS = [];
    end
else
    CLIBS = [];
end

COPT = ['COPTIMFLAGS=''' modelstruct.coptim ' -DNDEBUG'''];

try
    if(modelstruct.debug)
        eval(['mex -g ' CLIBS ' -output ' fullfile(odewrap_path,'models',filename,['lw_' filename]) ' ' fullfile(odewrap_path,'llhwrap.c')  includesstr objectsstr ])
    else
        eval(['mex ' COPT ' ' CLIBS ' -output ' fullfile(odewrap_path,'models',filename,['lw_' filename]) ' ' fullfile(odewrap_path,'llhwrap.c')  includesstr objectsstr ])
    end
catch
    CLIBS = []; % some systems need lrt libraries, some don't. I do not know how to check for this from matlab, so we'll just try both and see what works.
    if(modelstruct.debug)
        eval(['mex -g ' CLIBS ' -output ' fullfile(odewrap_path,'models',filename,['lw_' filename]) ' ' fullfile(odewrap_path,'llhwrap.c')  includesstr objectsstr ])
    else
        eval(['mex ' COPT ' ' CLIBS ' -output ' fullfile(odewrap_path,'models',filename,['lw_' filename]) ' ' fullfile(odewrap_path,'llhwrap.c')  includesstr objectsstr ])
    end
end

% $CFLAGS=''$CFLAGS -framework Accelerate''

end

function result = isnewer(ver1str,ver2str)

ver1Parts = getParts(ver1str);
ver2Parts = getParts(ver2str);
if ver2Parts(1) ~= ver1Parts(1)     % major version
    result = ver2Parts(1) < ver1Parts(1);
elseif ver2Parts(2) ~= ver1Parts(2) % minor version
    result = ver2Parts(2) < ver1Parts(2);
else                                  % revision version
    result = ver2Parts(3) < ver1Parts(3);
end
end

function parts = getParts(V)
parts = sscanf(V, '%d.%d.%d')';
if length(parts) < 3
    parts(3) = 0; % zero-fills to 3 elements
end
end