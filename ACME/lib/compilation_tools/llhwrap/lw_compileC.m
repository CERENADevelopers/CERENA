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


function lw_compileC(filename)
    
if(~ispc)
    llhwrap_path = strrep(which('lw_compileC.m'),'lw_compileC.m','');
    sundials_path = [strrep(which('lw_compileC.m'),'/lw_compileC.m','') '/sundials-2.5.0/']; % sundials 2.5.0
else
    llhwrap_path = strrep(which('lw_compileC.m'),'lw_compileC.m','');
    sundials_path = [strrep(which('lw_compileC.m'),'\lw_compileC.m','') '\sundials-2.5.0\']; % sundials 2.5.0
end

includesstr = '';
includesstr = strcat(includesstr, [' -I"' sundials_path 'include"']);
includesstr = strcat(includesstr, [' -I"' sundials_path 'src/cvodes"']);
includesstr = strcat(includesstr, [' -I"' llhwrap_path  'models/' filename '"']);
includesstr = strcat(includesstr, [' -I"' llhwrap_path  '"']);

sources = {
    %'src/cvodes/cvodes_lapack.c';
    'src/cvodes/cvodes_band.c';
    'src/cvodes/cvodes_bandpre.c';
    'src/cvodes/cvodes_bbdpre.c';
    'src/cvodes/cvodes_direct.c';
    'src/cvodes/cvodes_dense.c';
    'src/cvodes/cvodes_diag.c';
    'src/cvodes/cvodea.c';
    'src/cvodes/cvodes.c';
    'src/cvodes/cvodes_io.c';
    'src/cvodes/cvodea_io.c';
    'src/cvodes/cvodes_spils.c';
    'src/cvodes/cvodes_spbcgs.c';
    'src/cvodes/cvodes_spgmr.c';
    'src/cvodes/cvodes_sptfqmr.c';
    'src/sundials/sundials_band.c';
    'src/sundials/sundials_dense.c';
    'src/sundials/sundials_iterative.c';
    'src/sundials/sundials_nvector.c';
    'src/sundials/sundials_direct.c';
    'src/sundials/sundials_spbcgs.c';
    'src/sundials/sundials_spgmr.c';
    'src/sundials/sundials_sptfqmr.c';
    'src/sundials/sundials_math.c';
    'src/nvec_ser/nvector_serial.c';
    };
sourcesstr = '';
for j=1:length(sources)
    sourcesstr = strcat(sourcesstr, [' "' sundials_path sources{j} '"']);
end

% objects
objects = {
    %'cvodes_lapack.o';
    'cvodes_band.o';
    'cvodes_bandpre.o';
    'cvodes_bbdpre.o';
    'cvodes_direct.o';
    'cvodes_dense.o';
    'cvodes_diag.o';
    'cvodea.o';
    'cvodes.o';
    'cvodes_io.o';
    'cvodea_io.o';
    'cvodes_spils.o';
    'cvodes_spbcgs.o';
    'cvodes_spgmr.o';
    'cvodes_sptfqmr.o';
    'sundials_band.o';
    'sundials_dense.o';
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
    objects = strrep(objects, '.o', '.obj');
end

% compile directory
if(~exist([cd '/models/' mexext], 'dir'))
    mkdir([cd '/models/' mexext])
end

objectsstr = '';
for j=1:length(objects)
    objectsstr = strcat(objectsstr, [' ./models/' mexext '/' objects{j}]);
end

for j=1:length(sources)
     if(~exist(['./models/' mexext '/' objects{j}], 'file'))
        eval(['mex -g -c -outdir ./models/' mexext '/' includesstr ' ' sundials_path sources{j}]);
     end
end

eval(['mex -g -output ' llhwrap_path 'models/' filename '/' filename ' ' llhwrap_path 'llhwrap.c '  includesstr objectsstr ])