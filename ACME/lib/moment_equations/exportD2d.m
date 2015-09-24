function exportD2d(t,sys,options)

clear([options.filename '.def']);
fid = fopen([options.filename '.def'],'w');
str_description = ['DESCRIPTION \n' ...
    '"' options.modelname '" \n\n'...
    'PREDICTOR \n'...
    't               T   min         time   '	num2str(t(1)) '   '	num2str(t(2)) '\n\n'...
    'COMPARTMENTS \n'];
if ~isempty(sys.system.compartments)
    str_compartments = [char(sys.system.compartments) repmat('           V   "pl"      "vol."   ',length(sys.system.compartments),1) ...
        num2str(sys.system.volumes) repmat(' \n',length(sys.system.compartments),1)];
else
    str_compartments = ['\n'];
end
char_states = char(transpose(sys.MM.sym.state.moments));
char_states = char_states(10:end-3);
char_states = textscan(char_states,'%s','delimiter',',');
char_states = char(char_states{1});
str_states = [char_states repmat('        C   "nmol/l"      "conc."',length(sys.MM.sym.state.moments),1) ...
    repmat(' \n',length(sys.MM.sym.state.moments),1)];

str_input = ['\nINPUTS \n'...
    char(sys.system.input.name) '        C   "au"      "conc."      "' char(sys.system.input.function) '" \n\n'];
str_input = strrep(str_input,'time','t');

char_odes = char(transpose(sys.MM.sym.state.derivative));
char_odes = char_odes(10:end-3);
char_odes = ['"' char_odes '"'];
char_odes = strrep(char_odes,',','",');
char_odes = strrep(char_odes,', ',', "');
char_odes = textscan(char_odes,'%s','delimiter',',');
char_odes = char(char_odes{1});
str_odes = [char_odes  repmat(' \n',length(sys.MM.sym.state.moments),1)];

char_derived_name = sym([]);
for i = 1:length(sys.MM.sym.output.function)
    ind_output = sys.MM.sym.output.order(i,sys.MM.sym.output.order(i,:)~=0);
    if length(ind_output)==1
        char_derived_name = [char_derived_name; sys.system.output.name(ind_output)];
    elseif length(ind_output)>1
        var_output = ['C_',char(sys.system.output.name(ind_output(1)))];
        for j=2:length(ind_output)
            var_output = [var_output,['_',char(sys.system.output.name(ind_output(j)))]];
        end
        var_output = sym(var_output);
        char_derived_name = [char_derived_name; var_output];
    end
end
char_derived_name = char(transpose(char_derived_name));
char_derived_name = char_derived_name(10:end-3);
char_derived_name = textscan(char_derived_name,'%s','delimiter',',');
char_derived_name = char_derived_name{1};
for i = 1:length(char_derived_name)
    char_observables_1{i} = [char(char_derived_name(i)),'_au'];
end

char_derived_function = char(transpose(sys.MM.sym.output.function));
char_derived_function = char_derived_function(10:end-3);
char_derived_function = ['"' char_derived_function '"'];
char_derived_function = strrep(char_derived_function,',','",');
char_derived_function = strrep(char_derived_function,', ',', "');
char_derived_function = textscan(char_derived_function,'%s','delimiter',',');
char_derived_function = char(char_derived_function{1});
str_derived = [char(char_derived_name) repmat('        C   "nmol/l"      "conc."      ',length(sys.MM.sym.output.moments),1) ...
    char_derived_function repmat(' \n',length(sys.MM.sym.output.moments),1)];

str_observables_1 = [char(char_observables_1) repmat('        C   "au"      "conc."     1     0     "offset_',length(char_observables_1),1) ...
    char(char_observables_1) repmat(' + scale_',length(char_observables_1),1) char(char_observables_1) ...
    repmat(' * ',length(char_observables_1),1) char(char_observables_1) repmat('" \n',length(char_observables_1),1)];


% if ~isempty(ind_input)
if (sys.observeInput==1)
    str_observables_2 = [char(sys.system.input.name) '_au        C   "au"      "conc."     1     0     "' ...
        char(sys.system.input.name) '" \n'];
end

str_errors = [char(char_observables_1) repmat('      "sd_',length(char_observables_1),1) ...
    char(char_observables_1) repmat('" \n', length(char_observables_1),1)];

char_ic = char(transpose(sys.MM.sym.state.M0));
char_ic = char_ic(10:end-3);
char_ic = textscan(char_ic,'%s','delimiter',',');
char_ic = char(char_ic{1});
str_conditions_1 = [repmat('init_',size(char_states,1),1) char_states ...
    repmat('       "',size(char_states,1),1) char_ic repmat('" \n',size(char_states,1),1)];
str_conditions_2 = [char(sys.conditions) repmat(' \n',length(sys.conditions),1)];

fprintf(fid,str_description);
fprintf(fid,str_compartments');
fprintf(fid,'\nSTATES \n');
fprintf(fid,str_states');
fprintf(fid,str_input);
fprintf(fid,'\nODES \n');
fprintf(fid,str_odes');
fprintf(fid,'\nDERIVED \n');
fprintf(fid,str_derived');
if (sys.observables == 1)
    fprintf(fid,'\nOBSERVABLES \n');
    fprintf(fid,str_observables_1');
    if (sys.observeInput==1)
        fprintf(fid,str_observables_2);
    end
    fprintf(fid,'\nERRORS \n');
    fprintf(fid,str_errors');
end
fprintf(fid,'\nCONDITIONS \n');
fprintf(fid,str_conditions_1');
fprintf(fid,str_conditions_2');
fclose(fid);

% rehash to ensure that the function is known/ used.
rehash