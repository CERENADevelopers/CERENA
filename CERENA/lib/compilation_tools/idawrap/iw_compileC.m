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


function iw_compileC(filename)
    
if(~ispc)
    odewrap_path = strrep(which('iw_compileC.m'),'iw_compileC.m','');
    sundials_path = [strrep(which('iw_compileC.m'),'/iw_compileC.m','') '/sundials-2.5.0/']; % sundials 2.5.0
else
    odewrap_path = strrep(which('iw_compileC.m'),'iw_compileC.m','');
    sundials_path = [strrep(which('iw_compileC.m'),'\jw_compileC.m','') '\sundials-2.5.0\']; % sundials 2.5.0
end

includesstr = '';
includesstr = strcat(includesstr, [' -I"' sundials_path 'include"']);
includesstr = strcat(includesstr, [' -I"' sundials_path 'src/cvodes"']);
includesstr = strcat(includesstr, [' -I"' odewrap_path  'models/' filename '"']);
includesstr = strcat(includesstr, [' -I"' odewrap_path  '"']);

sources = {
    'src/idas/idaa_io.c';
    'src/idas/idaa.c';
    'src/idas/idas_band.c';
    'src/idas/idas_bbdpre.c';
    'src/idas/idas_dense.c';
    'src/idas/idas_direct.c';
    'src/idas/idas_ic.c';
    'src/idas/idas_io.c';
    'src/idas/idas_spbcgs.c';
    'src/idas/idas_spgmr.c';
    'src/idas/idas_spils.c';
    'src/idas/idas_sptfqmr.c';
    'src/idas/idas.c';
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
    'idaa_io.o';
    'idaa.o';
    'idas_band.o';
    'idas_bbdpre.o';
    'idas_dense.o';
    'idas_direct.o';
    'idas_ic.o';
    'idas_io.o';
    'idas_spbcgs.o';
    'idas_spgmr.o';
    'idas_spils.o';
    'idas_sptfqmr.o';
    'idas.o';
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
        eval(['mex -c -outdir ./models/' mexext '/' includesstr ' ' sundials_path sources{j}]);
     end
end

if(~exist(['models/' mexext '/' objects{j}], 'file'))
    eval(['mex -c -outdir models/' mexext '/' includesstr ' ' filename '/' filename  '.c']);
end

eval(['mex -output ' odewrap_path 'models/' filename '/' filename ' ' odewrap_path 'idawrap.c '  includesstr objectsstr ])
%eval(['mex -output ' odewrap_path 'models/' filename '/' filename '_ASA ' odewrap_path 'idawrap_ASA.c '  includesstr objectsstr ])