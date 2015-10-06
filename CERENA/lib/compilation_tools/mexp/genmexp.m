% genmexp generates symbolically egenerates model definition files for the System size expansion and moment equations,
% which can be subsequently by the d2d-toolbox or to generate simulation .mex files using the cvodewrapper
%
% USAGE:
% ======
% genmexp(modelname,modelfile,expansion)
%
% INPUTS:
% =======
% modelname ... specifies the name of the model. this specification will be used for the names of the model definition
%               files
% modelfile ... specifies the location of the model files which define the reaction network
% expansion ... specifies the type of moment expansion scheme to be used. options are
%   'RRE' ... Reaction Rate Equations. First order System Size Expansion of the mean
%   'EMRE' ... Effective Mesoscopic Rate Equations. Second order System Size Expansion of the mean
%   'LNA' ... Linear Noise Approximation. First order System Size Expansion of the mean and variance
%   'IOS' ... Inverse Omega Square Method. Second order System Size Expansion of the mean and variance
%   'MEC_X_CC_Y_SC' ... Centered Moment Equations for moments of order 1
%   to X with closure CC and output order Y and scale SC
%                  The options for SC are
%                  c ... concentration
%                  a ... absolute
%                  The options for the closure CC are
%                  MF ... Mean field
%                  LD ... Low Dispersion
%                  ZC ... Zero Cumulants
%                  DM ... Derivative Matching
%                  LN ... Log Normal
%   'MEUC_X_CC_Y' ... Unentered Moment Equations for moments of order 1 to X with closure CC and output order Y
%                  The options for the closure CC are
%                  MF ... Mean field
%                  LD ... Low Dispersion
%                  ZC ... Zero Cumulants
%                  DM ... Derivative Matching
%   'CMEC_X_CC_Y_SC' ... Centered Conditional Moment Equations for moments of order 1
%   to X with closure CC and output order Y and scale SC
%                  The options for SC are
%                  c ... concentration
%                  a ... absolute
%                  The options for the closure CC are
%                  MF ... Mean field
%                  LD ... Low Dispersion
%                  ZC ... Zero Cumulants
%                  DM ... Derivative Matching

function [System] = genmexp(modelname,modelfile,expansion)

% check directory of file
modeldir = fileparts(modelfile);

if(strcmp(modeldir,''))
    if strcmp(modelfile(end-1:end),'.m')
        modeldir =  strrep(which(modelfile),modelfile,'');
    else
        modeldir =  strrep(which(modelfile),[modelfile '.m'],'');
    end
elseif(strcmp(modeldir(1),'.'))
    % path was provided as relative path
    modeldir = [ pwd modeldir(2:end) ];
end

tmpstr = which('cvodewrap.m');
cvodewrap_flag = not(strcmp(tmpstr,''));
if(cvodewrap_flag)
    disp('cvodewrap found in path! Will generate files')
    if(exist([modeldir expansion '_' modelname '_syms.m'],'file'))
        reply = input('Modelfile with the same name already exists for cvodewrap, do you want to it overwrite? (y/n) [n]', 's');
        if(not(any([strcmp(reply,'y'),strcmp(reply,'Y'),strcmp(reply,'yes')])))
            cvodewrap_flag = 0;
        end
    end
else
    disp('cvodewrap not found in path! Will not generate files')
end

tmpstr = which('idawrap.m');
idawrap_flag = not(strcmp(tmpstr,''));
if(idawrap_flag)
    disp('idawrap found in path! Will generate files')
    if(exist([modeldir expansion '_' modelname '_syms.m'],'file'))
        reply = input('Modelfile with the same name already exists for idawrap, do you want to it overwrite? (y/n) [n]', 's');
        if(not(any([strcmp(reply,'y'),strcmp(reply,'Y'),strcmp(reply,'yes')])))
            idawrap_flag = 0;
        end
    end
else
    disp('idawrap not found in path! Will not generate files')
end

% check for d2d
tmpstr = which('arInit.m');
d2d_flag = not(strcmp(tmpstr,['']));
if(d2d_flag)
    reply = input('Do you want d2d-files to be generated?(y/n) [n]', 's');
    if(not(any([strcmp(reply,'y'),strcmp(reply,'Y'),strcmp(reply,'yes')])))
        d2d_flag = 0;
    else
        disp('d2d-toolbox found in path! Will generate files')
        if(exist([modeldir expansion '_' modelname '.def' ],'file'))
            reply = input('Modelfile with the same name already exists for d2d-toolbox, do you want to it overwrite? (y/n) [n]', 's');
            if(not(any([strcmp(reply,'y'),strcmp(reply,'Y'),strcmp(reply,'yes')])))
                d2d_flag = 0;
            end
        end
    end
else
    disp('d2d-toolbox not found in path! Will not generate files')
end


reply = input('Do you want the MATLAB simulation files to be generated? (y/n) [n]', 's');
if(any([strcmp(reply,'y'),strcmp(reply,'Y'),strcmp(reply,'yes')]))
    matlab_flag = 1;
else
    matlab_flag = 0;
end

% if neither of them is found, abort. Maybe we should then just generate respective files in the same directory?
if(not(any([idawrap_flag,cvodewrap_flag,d2d_flag,matlab_flag])))
    disp('No files to generate! Aborting ...');
    return
end

% check modelfile
tmpstr = which(modelfile);
if(strcmp(tmpstr,['']))
    disp('Specified modelfile was not found in path! Aborting ...');
end

%% SSE
if(any([strcmp(expansion,'RRE'),strcmp(expansion,'EMRE'),strcmp(expansion,'LNA'),strcmp(expansion,'IOS')]))
    
    idawrap_flag = 0;
    switch(expansion)
        case 'RRE'
            RREFLAG = 1;
            LNAFLAG = 0;
            EMREFLAG = 0;
            IOSFLAG = 0;
        case 'LNA'
            RREFLAG = 1;
            LNAFLAG = 1;
            EMREFLAG = 0;
            IOSFLAG = 0;
        case 'EMRE'
            RREFLAG = 1;
            LNAFLAG = 1;
            EMREFLAG = 1;
            IOSFLAG = 0;
        case 'IOS'
            RREFLAG = 1;
            LNAFLAG = 1;
            EMREFLAG = 1;
            IOSFLAG = 1;
    end
    
    %% INITIALIZATION
    
    eval(modelfile);
    switch System.scaleIndicator
        case 'microscopic'
%             System.scaleConversion = 'Micro_to_Macro';
            System.scaleConversion = 'None';
        case 'macroscopic'
%             System.scaleConversion = 'None';
            System.scaleConversion = 'Macro_to_Micro';
    end
    System = completeSystem(System);
    System.name = expansion;
    %% Generate SSE equations
    System = getSSE(System);
    par = System.parameter.variable;
    kappa = System.kappa;
    nmx = System.SSE.nmx;
    states = System.SSE.states;
    xdot = System.SSE.xdot;
    x0 = System.SSE.x0;
    output = System.SSE.output;
end

%% MOMENT EQUATIONS
if(strcmp(expansion(1:2),'ME'))
    idawrap_flag = 0;
    if(strcmp(expansion(3:4),'C_'))
        % CENTERED MOMENTS
        
        % parse moment order
        options.moment_order = str2double(expansion(5));
        scale = expansion(12);
        switch(scale)
            case 'c'
                options.scale = 'concentration';
            case 'a'
                options.scale = 'absolute';
        end
        
        % parse closure scheme
        closure = expansion(7:8);
        switch(closure)
            case 'MF'
                options.moment_closure = 'mean-field';
            case 'LD'
                options.moment_closure = 'low-dispersion';
            case 'ZC'
                options.moment_closure = 'zero-cumulants';
            case 'DM'
                options.moment_closure = 'derivative-matching';
            case 'LN'
                options.moment_closure = 'log-normal';
            otherwise
                disp('Invalid Expansion String')
                return
        end
        if length(expansion) > 12
            options.moment_order_max = str2double(expansion(14));
        end
        
        options.moment_order_output = str2double(expansion(10));
        
        eval(modelfile);
        switch System.scaleIndicator
            case 'microscopic'
                System.scaleConversion = 'None';
            case 'macroscopic'
                System.scaleConversion = 'Macro_to_Micro';
        end
        System = completeSystem(System);
        System.name = 'MM';
        
        MM = getMM_centered(System,options);
        System.MM = MM;
    elseif(strcmp(expansion(3:5),'UC_'))
        % UNCENTERED MOMENTS
        
        % parse moment order
        options.moment_order = str2double(expansion(6));
        
        % parse closure scheme
        closure = expansion(8:9);
        switch(closure)
            case 'MF'
                options.moment_closure = 'mean-field';
            case 'LD'
                options.moment_closure = 'low-dispersion';
            case 'ZC'
                options.moment_closure = 'zero-cumulants';
            case 'DM'
                options.moment_closure = 'derivative-matching';
            otherwise
                disp('Invalid Expansion String')
                return
        end
        
        options.moment_order_output = str2double(expansion(end));
        
        eval(modelfile);
        System = completeSystem(System);
        
        MM = getMM_uncentered(System,options);
        System.MM = MM;
    end
    
    par = System.parameter.variable;
    kappa = System.kappa;
    kappavar = System.kappa.variable;
    states = MM.sym.state.moments;
    System.states = states;
    System.MM.output_order = options.moment_order_output;
    xdot = MM.sym.state.derivative;
    output = MM.sym.output.function;
    x0 = MM.sym.state.M0;
    
    % replace parameters
    c = sym('c',[length(par),1]);
    c_kappa = sym('c_kappa',[length(kappavar),1]);
    xdot = subs(xdot,c,par);
    x0 = subs(x0,c,par);
    output = subs(output,c,par);
    xdot = subs(xdot,c_kappa,kappavar);
    x0 = subs(x0,c_kappa,kappavar);
    output = subs(output,c_kappa,kappavar);
    syms t
    xdot = subs(xdot,System.time,t);
    x0 = subs(x0,System.time,t);
    if strcmp(options.moment_closure,'derivative-matching')
        ind_0 = find(x0 == 0);
        x0(ind_0) = 1e-10;
    end
    output = subs(output,System.time,t);
    
end

%% CONDITIONAL MOMENT EQUATIONS
if(strcmp(expansion(1:3),'CME'))
    cvodewrap_flag = 0;
    if(strcmp(expansion(4:5),'C_'))
        % CENTERED MOMENTS
        
        % parse moment order
        options.moment_order = str2double(expansion(6));
        scale = expansion(13);
        switch(scale)
            case 'c'
                options.scale = 'concentration';
            case 'a'
                options.scale = 'absolute';
        end
        
        % parse closure scheme
        closure = expansion(8:9);
        switch(closure)
            case 'MF'
                options.moment_closure = 'mean-field';
            case 'LD'
                options.moment_closure = 'low-dispersion';
            case 'ZC'
                options.moment_closure = 'zero-cumulants';
            case 'DM'
                options.moment_closure = 'derivative-matching';
            otherwise
                disp('Invalid Expansion String')
                return
        end
        
        options.moment_order_output = str2double(expansion(11));
        if length(expansion) > 13
            options.moment_order_max = str2double(expansion(15));
        end
        
        eval(modelfile);
        switch System.scaleIndicator
            case 'microscopic'
                System.scaleConversion = 'None';
            case 'macroscopic'
                System.scaleConversion = 'Macro_to_Micro';
        end
        System = completeSystem(System);
        System.name = 'MCM';
        
        CMM = getcMM_centered(System,System.state.xmin,System.state.xmax,[],options);
        System.CMM = CMM;
    elseif(strcmp(expansion(3:5),'UC_'))
        % UNCENTERED MOMENTS
        
        % parse moment order
        options.moment_order = str2double(expansion(6));
        
        
        % parse closure scheme
        closure = expansion(8:9);
        switch(closure)
            case 'MF'
                options.moment_closure = 'mean-field';
            case 'LD'
                options.moment_closure = 'low-dispersion';
            case 'ZC'
                options.moment_closure = 'zero-cumulants';
            case 'DM'
                options.moment_closure = 'derivative-matching';
            otherwise
                disp('Invalid Expansion String')
                return
        end
        
        options.moment_order_output = str2double(expansion(end));
        
        eval(modelfile);
        System = completeSystem(System);
        
        CMM = getcMM_uncentered(System,options);
        System.CMM = CMM;
    end
    
    par = System.parameter.variable;
    kappa = System.kappa;
    kappavar = System.kappa.variable;
    states = CMM.state.sym.M;
    System.states = states;
    System.CMM.output_order = options.moment_order_output;
    nmx = size(System.CMM.state_moments.order,1);
    f = CMM.derivatives.sym.VF;
    MMat = CMM.derivatives.sym.mass_matrix;
    output = CMM.output.function;
    x0 = CMM.state.sym.M0;
    if strcmp(options.moment_closure,'derivative-matching')
        if ~isempty(x0(x0==0))
            x0(x0==0) = 1e-10;
        end
    end
    try
        dx0 = subs(f,states,x0)./ subs(MMat,states,x0);
    catch
        if ~isempty(x0(x0==0))
            x0(x0==0) = 1e-10;
        end
        dx0 = subs(f,states,x0)./ subs(MMat,states,x0);
    end
    
    % replace parameters
    c = sym('c',[length(par),1]);
    c_kappa = sym('c_kappa',[length(kappavar),1]);
    f = subs(f,c,par);
    MMat = subs(MMat,c,par);
    x0 = subs(x0,c,par);
    dx0 = subs(dx0,c,par);
    output = subs(output,c,par);
    f = subs(f,c_kappa,kappavar);
    MMat = subs(MMat,c_kappa,kappavar);
    x0 = subs(x0,c_kappa,kappavar);
    dx0 = subs(dx0,c_kappa,kappavar);
    output = subs(output,c_kappa,kappavar);
    syms t
    f = subs(f,System.time,t);
    MMat = subs(MMat,System.time,t);
    x0 = subs(x0,System.time,t);
    dx0 = subs(dx0,System.time,t);
    output = subs(output,System.time,t);
end

if(strcmp(expansion,'FSP'))
    idawrap_flag = 0;
    eval(modelfile);
    switch System.scaleIndicator
        case 'microscopic'
            System.scaleConversion = 'None';
        case 'macroscopic'
            System.scaleConversion = 'Macro_to_Micro';
    end
    System = completeSystem(System);
    System.name = 'FSP';
    [System] = getFSP(System,System.state.xmin,System.state.xmax);
    states = sym('p',[size(System.index,1),1]);
    System.states = states;
    % Initial condition
    model.p0 = zeros(size(System.index,1),1);
    [~,ix] =ismember(System.state.mu0',System.index,'Rows');
    model.p0(ix) = 1;
    par = System.parameter.variable;
    kappa = System.kappa;
    kappavar = System.kappa.variable;
    xdot = sym(zeros(length(states),1));
    ind100 = 1:100:length(states);
    if exist([modeldir '_FSP_xdot.mat'])>0
        load('FSP_xdot')
    else
        for i = ind100(1:end-1)
            row_i = zeros(100,length(states));
            tic
            for j = 1:length(System.A_i)
                row_i = row_i + System.sigma_A_i{j}(par)*System.A_i{j}(i:i+99,:);
            end
            toc
            tic
            xdot(i:i+99) = row_i * states;
            toc
        end
        row_i = zeros(length(states(ind100(end):end)),length(states));
        for j = 1:length(System.A_i)
            row_i = row_i + System.sigma_A_i{j}(par)*System.A_i{j}(ind100(end):end,:);
        end
    end
    xdot(ind100(end):end) = row_i * states;
    [output,output_ind] = getMomentsFSP_centered_sym(states,System.index,2);
    System.output.order = output_ind;
    x0 = model.p0;
    c = sym('c',[length(par),1]);
    c_kappa = sym('c_kappa',[length(kappavar),1]);
    xdot = subs(xdot,c,par);
    x0 = subs(x0,c,par);
    output = subs(output,c,par);
    xdot = subs(xdot,c_kappa,kappavar);
    x0 = subs(x0,c_kappa,kappavar);
    output = subs(output,c_kappa,kappavar);
    syms t
    xdot = subs(xdot,System.time,t);
    x0 = subs(x0,System.time,t);
    output = subs(output,System.time,t);
end

%% Saving the System struct
if ~exist([modeldir,'results'],'dir')
    mkdir([modeldir,'results'])
    addpath([modeldir,'results'])
end
save([modeldir,'results/',modelname,'_',expansion,'_System.mat'],'System')
rehash
%% MODEL FILES CVODEWRAP

if(exist('xdot','var')||exist('f','var'))
    %% ode15s function
    if (matlab_flag)
        if(cvodewrap_flag)
            
            xdot_mat = xdot;
            output_mat = output;
            x0_mat = transpose(x0);
            % Symbolic array for the prespecified initial conditions of mu0 in the modelDef
            fmu0_sym = sym('fmu0',[length(System.state.mu0),1]);
            % Symbolic array for the prespecified initial conditions of C0 in the modelDef
            fC0_sym = sym('fC0',[length(System.state.C0),1]);
            for j = 1:length(states)
                x{j,1} = sprintf('x(%i)',j);
            end
            for j = 1:length(par)
                Theta{j,1} = sprintf('Theta(%i)',j);
            end
            for j = 1:length(kappa.variable)
                kap{j,1} = sprintf('kappa(%i)',j);
            end
            x = sym(x);
            Theta = sym(Theta);
            theta = subs(Theta,'Theta','theta');
            kap = sym(kap);
            xdot_mat = mysubs(xdot_mat,states,x);
            xdot_mat = mysubs(xdot_mat,par,theta);
            xdot_mat = mysubs(xdot_mat,kappa.variable,kap);
            if(any([strcmp(expansion,'RRE'),strcmp(expansion,'EMRE'),strcmp(expansion,'LNA'),strcmp(expansion,'IOS')]))
                char_xdot_mat = char(xdot_mat);
            else
                char_xdot_mat = char(transpose(xdot_mat));
            end
            str_xdot_mat = char_xdot_mat(9:end-2);
            output_mat = mysubs(output_mat,states,x);
            output_mat = mysubs(output_mat,par,theta);
            output_mat = mysubs(output_mat,kappa.variable,kap);
            if(any([strcmp(expansion,'RRE'),strcmp(expansion,'EMRE'),strcmp(expansion,'LNA'),strcmp(expansion,'IOS')]))
                char_output_mat = char(output_mat);
            else
                char_output_mat = char(transpose(output_mat));
            end
            str_output_mat = char_output_mat(9:end-2);
            str_output_mat = strrep(str_output_mat,',',';');
            for j = 1:length(states)
                str_output_mat = strrep(str_output_mat,char(x(j)),['x(:,',num2str(j,'%d'),')']);
            end
            str_output_mat = strrep(str_output_mat,'*','.*');
            str_output_mat = strrep(str_output_mat,'/','./');
            str_output_mat = strrep(str_output_mat,'^','.^');
            if isempty(str_output_mat)
                str_output_mat = '[]';
            end
            x0_mat = mysubs(x0_mat,fmu0_sym,System.state.fmu0);
            x0_mat = mysubs(x0_mat,fC0_sym,System.state.fC0);
            x0_mat = mysubs(x0_mat,states,x);
            x0_mat = mysubs(x0_mat,par,theta);
            x0_mat = mysubs(x0_mat,kappa.variable,kap);
            char_x0_mat = char(x0_mat);
            str_x0_mat = char_x0_mat(9:end-2);
            
            clear([modeldir expansion '_' modelname '_matlab.m']);
            fid = fopen([modeldir expansion '_' modelname '_matlab.m'],'w');
            if (any([strcmp(expansion,'RRE'),strcmp(expansion,'EMRE'),strcmp(expansion,'LNA'),strcmp(expansion,'IOS')]))
                str_FUN_Sim = ['%% function [tout,X,MX,Y] = ' expansion '_' modelname '_matlab(t,theta,kappa) \n'...
                    ];
            else
                str_FUN_Sim = ['%% function [tout,X,Y] = ' expansion '_' modelname '_matlab(t,theta,kappa) \n'...
                    ];
            end
            str_FUN_Sim = [str_FUN_Sim ...
                'function varargout = ' expansion '_' modelname '_matlab(varargin) \n\n'...
                't = varargin{1};\n'...
                'theta = varargin{2};\n'...
                'if(nargin>=3)\n'...
                '    kappa=varargin{3};\n'...
                '   if(length(kappa)==',num2str(kappa.nk1),')\n'...
                '    kappa(',num2str(kappa.nk1+1),':',num2str(length(kappa.variable)),')=0;\n'...
                '   end\n'...
                'else\n'...
                '    kappa = zeros(1,',num2str(length(kappa.variable)),');\n'...
                'end\n'...
                '%% Initial conditions\n'...
                'x0 = x0fun(theta,kappa);\n'...
                '\n'...
                '%% Simulation\n'...
                '[tout,X] = ode15s(@(t,x) rhs(t,x,theta,kappa),t,x0);\n'...
                'Y = rhsO(t,X,theta,kappa);\n'...
                '\n'...
                '%% Assign output\n'...
                'varargout{1} = tout;\n'...
                'if nargout >= 2\n'...
                'varargout{2} = X;\n'...
                'end\n'...
                ];
            if (any([strcmp(expansion,'RRE'),strcmp(expansion,'EMRE'),strcmp(expansion,'LNA'),strcmp(expansion,'IOS')]))
                str_FUN_Sim = [str_FUN_Sim ...
                    'if nargout >= 3\n'...
                    '    %% Moments of species\n'...
                    '    varargout{3} = Y(:,1:',num2str(nmx),');\n'...
                    'end\n'...
                    'if nargout >= 4\n'...
                    '    %% Moments of output variables\n'...
                    '    varargout{4} = Y(:,',num2str(nmx+1),':end);\n'...
                    'end\n'...
                    'if nargout >= 5\n'...
                    '    error(''Too many output arguments.'');\n'...
                    'end\n'...
                    '\n\n'...
                    ];
            else
                str_FUN_Sim = [str_FUN_Sim ...
                    'if nargout >= 3\n'...
                    '    varargout{3} = Y;\n'...
                    'end\n'...
                    'if nargout >= 4\n'...
                    '    error(''Too many output arguments.'');\n'...
                    'end\n'...
                    '\n\n'...
                    ];
            end
            str_FUN_Sim = [str_FUN_Sim ...
                '%%%% RIGHT-HAND SIDE\n'...
                'function [dxdt] = rhs(t,x,theta,kappa) \n\n'...
                'dxdt = ' strrep(str_xdot_mat,',',';...\n        ') ';\n\n'...
                '%%%% OUTPUT MAP\n'...
                'function y = rhsO(t,x,theta,kappa) \n\n'...
                'y = ' strrep(str_output_mat,';',',...\n        ') ';\n'...
                '\n\n'...
                '%%%% INITIAL CONDITIONS FOR STATE\n'...
                'function x0 = x0fun(theta,kappa) \n\n'...
                'x0 = ' strrep(str_x0_mat,',',';...\n      ') ';\n'...
                '\n\n'...
                ];
            fprintf(fid,str_FUN_Sim);
            fclose(fid);
            rehash
        end
        if(idawrap_flag)
            f_mat = f;
            MMat_mat = MMat;
            output_mat = output;
            x0_mat = transpose(x0);
            % Symbolic array for the prespecified initial conditions of mu0 in the modelDef
            fmu0_sym = sym('fmu0',[length(System.state.mu0),1]);
            % Symbolic array for the prespecified initial conditions of C0 in the modelDef
            fC0_sym = sym('fC0',[length(System.state.C0),1]);
            for j = 1:length(states)
                x{j,1} = sprintf('x(%i)',j);
            end
            for j = 1:length(par)
                Theta{j,1} = sprintf('Theta(%i)',j);
            end
            for j = 1:length(kappa.variable)
                kap{j,1} = sprintf('kappa(%i)',j);
            end
            x = sym(x);
            Theta = sym(Theta);
            theta = subs(Theta,'Theta','theta');
            kap = sym(kap);
            f_mat = mysubs(f_mat,states,x);
            f_mat = mysubs(f_mat,par,theta);
            f_mat = mysubs(f_mat,kappa.variable,kap);
            char_f_mat = char(transpose(f_mat));
            str_f_mat = char_f_mat(9:end-2);
            MMat_mat = mysubs(MMat_mat,states,x);
            MMat_mat = mysubs(MMat_mat,par,theta);
            MMat_mat = mysubs(MMat_mat,kappa.variable,kap);
            char_MMat_mat = char(transpose(MMat_mat));
            str_MMat_mat = char_MMat_mat(9:end-2);
            output_mat = mysubs(output_mat,states,x);
            output_mat = mysubs(output_mat,par,theta);
            output_mat = mysubs(output_mat,kappa.variable,kap);
            char_output_mat = char(output_mat);
            str_output_mat = char_output_mat(9:end-2);
            str_output_mat = strrep(str_output_mat,',',';');
            for j = 1:length(states)
                str_output_mat = strrep(str_output_mat,char(x(j)),['x(:,',num2str(j,'%d'),')']);
            end
            str_output_mat = strrep(str_output_mat,'*','.*');
            str_output_mat = strrep(str_output_mat,'/','./');
            str_output_mat = strrep(str_output_mat,'^','.^');
            if isempty(str_output_mat)
                str_output_mat = '[]';
            end
            x0_mat = mysubs(x0_mat,fmu0_sym,System.state.fmu0);
            x0_mat = mysubs(x0_mat,fC0_sym,System.state.fC0);
            x0_mat = mysubs(x0_mat,states,x);
            x0_mat = mysubs(x0_mat,par,theta);
            x0_mat = mysubs(x0_mat,kappa.variable,kap);
            char_x0_mat = char(x0_mat);
            str_x0_mat = char_x0_mat(9:end-2);
            
            clear([modeldir expansion '_' modelname '_matlab.m']);
            fid = fopen([modeldir expansion '_' modelname '_matlab.m'],'w');
            str_FUN_Sim = ['%% function [tout,X,MX,Y] = ' expansion '_' modelname '_matlab(t,theta,kappa) \n'...
                'function varargout = ' expansion '_' modelname '_matlab(varargin) \n\n'...
                't = varargin{1};\n'...
                'theta = varargin{2};\n'...
                'if(nargin>=3)\n'...
                '    kappa=varargin{3};\n'...
                '   if(length(kappa)==',num2str(kappa.nk1),')\n'...
                '    kappa(',num2str(kappa.nk1+1),':',num2str(length(kappa.variable)),')=0;\n'...
                '   end\n'...
                'else\n'...
                '    kappa = zeros(1,',num2str(length(kappa.variable)),');\n'...
                'end\n'...
                'x0 = [];\n'...
                '%% Initial conditions\n'...
                '    x0 = x0fun(theta,kappa);\n'...
                '\n'...
                '%% Simulation\n'...
                'options = odeset(''Mass'',@(t,x) MMat(t,x,theta,kappa),...\n'...
                '''MStateDependence'',''strong'',''RelTol'',1e-8,''AbsTol'',1e-8);\n'...
                '[tout,X] = ode15s(@(t,x) rhs(t,x,theta,kappa),t,x0,options);\n'...
                'Y = rhsO(t,X,theta,kappa);\n'...
                '\n'...
                '%% Assign output\n'...
                'varargout{1} = tout;\n'...
                'if nargout >= 2\n'...
                'varargout{2} = X;\n'...
                'end\n'...
                'if nargout >= 3\n'...
                '    %% Overall moments of species\n'...
                '    varargout{3} = Y(:,1:',num2str(nmx),');\n'...
                'end\n'...
                'if nargout >= 4\n'...
                '    %% Moments of output variables\n'...
                '    varargout{4} = Y(:,',num2str(nmx+1),':end);\n'...
                'end\n'...
                'if nargout >= 5\n'...
                '    error(''Too many output arguments.'');\n'...
                'end\n'...
                '\n\n'...
                '%%%% RIGHT-HAND SIDE\n'...
                'function [dxdt] = rhs(t,x,theta,kappa) \n\n'...
                'dxdt = ' strrep(str_f_mat,',',';...\n        ') ';\n\n'...
                '%%%% MASS MATRIX\n'...
                'function [M] = MMat(t,x,theta,kappa) \n\n'...
                'M = spdiags(' strrep(str_MMat_mat,',',';...\n        ') ',0,' num2str(length(f),'%d') ',' num2str(length(f),'%d') ');\n\n'...
                '%%%% OUTPUT MAP\n'...
                'function y = rhsO(t,x,theta,kappa) \n\n'...
                'y = ' strrep(str_output_mat,';',',...\n        ') ';\n'...
                '\n\n'...
                '%%%% INITIAL CONDITIONS FOR STATE\n'...
                'function x0 = x0fun(theta,kappa) \n\n'...
                'x0 = ' strrep(str_x0_mat,',',';...\n      ') ';\n'...
                '\n\n'...
                ];
            fprintf(fid,str_FUN_Sim);
            fclose(fid);
            rehash
        end
    end
    %%
    if(cvodewrap_flag || idawrap_flag)
        model = expansion;
        clear([modeldir expansion '_' modelname '_syms.m']);
        fid = fopen([modeldir expansion '_' modelname '_syms.m'],'w');
        %         fprintf(fid,['function [model] = ' model '_' modelname '_syms()\n']);
        fprintf(fid,['%% function [model] = ' model '_' modelname '_syms(f0_user)\n']);
        fprintf(fid,['function [model] = ' model '_' modelname '_syms(varargin)\n']);
        fprintf(fid,['\n']);
        fprintf(fid,['%% CVODES OPTIONS\n']);
        fprintf(fid,['\n']);
        fprintf(fid,['% absolute tolerance\n']);
        fprintf(fid,['model.atol = 1e-8;\n']);
        fprintf(fid,['% relative tolerance\n']);
        fprintf(fid,['model.rtol = 1e-8;\n']);
        fprintf(fid,['% maximal number of steps\n']);
        fprintf(fid,['model.maxsteps = 1e4;\n']);
        fprintf(fid,['\n']);
        fprintf(fid,['%% STATES\n']);
        fprintf(fid,['\n']);
        
        char_states = char(transpose(states));
        
        fprintf(fid,['syms ' strrep(char_states(10:end-3),',','') '\n']);
        fprintf(fid,['\n']);
        fprintf(fid,['x = [\n']);
        fprintf(fid,[char_states(10:end-3) ' ...\n']);
        fprintf(fid,['];\n']);
        fprintf(fid,['\n']);
        fprintf(fid,['%% PARAMETERS\n']);
        fprintf(fid,['\n']);
        
        tmpstr = '';
        for j=1:length(par);
            tmpstr = [tmpstr char(par(j)) ' '];
        end
        fprintf(fid,['syms ' tmpstr '\n']);
        fprintf(fid,['\n']);
        fprintf(fid,['%% KAPPA (constant parameters)\n']);
        fprintf(fid,['\n']);
        
        tmpstr = '';
        for j=1:length(kappa.variable);
            tmpstr = [tmpstr char(kappa.variable(j)) ' '];
        end
        fprintf(fid,['syms ' tmpstr '\n']);
        fprintf(fid,['\n']);
        fprintf(fid,['syms t\n']);
        fprintf(fid,['\n']);
        
        tmpstr = '';
        for j=1:length(par)-1;
            tmpstr = [tmpstr char(par(j)) ','];
        end
        tmpstr = [tmpstr char(par(end)) ''];
        fprintf(fid,['p = [' tmpstr '];\n']);
        fprintf(fid,['\n']);
        
        tmpstr = '';
        for j=1:length(kappa.variable)-1;
            tmpstr = [tmpstr char(kappa.variable(j)) ','];
        end
        tmpstr = [tmpstr char(kappa.variable(end)) ''];
        fprintf(fid,['k = [' tmpstr '];\n']);
        fprintf(fid,['\n']);
        
        fprintf(fid,['if nargin > 0\n']);
        fprintf(fid,['   f0_user = varargin{1};\n']);
        fprintf(fid,['   if ~isnumeric(f0_user)\n']);
        fprintf(fid,['      p_user = setdiff(symvar(f0_user),p);\n']);
        fprintf(fid,['      %% ADDITIONAL PARAMETERS IN INITIAL CONDITIONS\n']);
        fprintf(fid,['      p = [p,p_user];\n']);
        fprintf(fid,['   end\n']);
        for j = 1:length(System.state.fmu0)
            fprintf(fid,['	fmu0',num2str(j),' = f0_user(',num2str(j),'); \n']);
        end
        for j = 1:length(System.state.fC0)
            fprintf(fid,['	fC0',num2str(j),' = f0_user(',num2str(length(System.state.fmu0)+j),'); \n']);
        end
        fprintf(fid,['else\n']);
        for j = 1:length(System.state.fmu0)
            fprintf(fid,['	fmu0',num2str(j),' = ',char(System.state.fmu0(j)),'; \n']);
        end
        for j = 1:length(System.state.fC0)
            fprintf(fid,['	fC0',num2str(j),' = ',char(System.state.fC0(j)),'; \n']);
        end
        fprintf(fid,['end\n']);
        
        fprintf(fid,['%% INPUT \n']);
        fprintf(fid,['\n']);
        fprintf(fid,['u = sym.empty(0,0);\n']);
        fprintf(fid,['\n']);
        fprintf(fid,['%% SYSTEM EQUATIONS\n']);
        fprintf(fid,['\n']);
        
        if (cvodewrap_flag)
            fprintf(fid,['xdot = sym(zeros(size(x)));\n']);
            fprintf(fid,['\n']);
            
            Nx = length(xdot);
            
            for j = 1:Nx
                %                 fprintf(fid,['xdot(' num2str(j) ') = ' char(simplify(xdot(j))) ';\n']);
                fprintf(fid,['xdot(' num2str(j) ') = ' char(xdot(j)) ';\n']);
            end
        end
        if (idawrap_flag)
            fprintf(fid,['f = sym(zeros(size(x)));\n']);
            fprintf(fid,['\n']);
            
            Nf = length(f);
            
            for j = 1:Nf
                %                 fprintf(fid,['f(' num2str(j) ') = ' char(simplify(f(j))) ';\n']);
                fprintf(fid,['f(' num2str(j) ') = ' char(f(j)) ';\n']);
            end
            tmpstr = char(simplify(MMat));
            fprintf(fid, ['M = diag(sym(',tmpstr(8:end-1), '));\n']);
        end
        
        
        fprintf(fid,['%% INITIAL CONDITIONS\n']);
        fprintf(fid,['\n']);
        fprintf(fid,['x0 = sym(zeros(size(x)));\n']);
        fprintf(fid,['\n']);
        
        if(any([strcmp(expansion,'RRE'),strcmp(expansion,'EMRE'),strcmp(expansion,'LNA'),strcmp(expansion,'IOS')]))
            for j = 1:length(System.state.mu0)
                if( ~isequaln(System.state.mu0(j),sym(0)))
                    fprintf(fid,['x0(' num2str(j) ') = ' char(simplify(x0(j))) ';\n']);
                end
            end
        else
            for j = 1:length(x0)
                %                 fprintf(fid,['x0(' num2str(j) ') = ' char(simplify(x0(j))) ';\n']);
                tmpchar = char(x0(j));
                fprintf(fid,['x0(' num2str(j) ') = ' tmpchar ';\n']);
            end
        end
        if (idawrap_flag)
            fprintf(fid,['dx0 = sym(zeros(size(x)));\n']);
            for j = 1:length(x0)
                %                 fprintf(fid,['dx0(' num2str(j) ') = ' char(simplify(dx0(j))) ';\n']);
                fprintf(fid,['dx0(' num2str(j) ') = ' char(dx0(j)) ';\n']);
            end
            
        end
        
        fprintf(fid,['\n']);
        fprintf(fid,['%% OBSERVABLES\n']);
        fprintf(fid,['\n']);
        fprintf(fid,['y = sym(zeros(' num2str(length(output)) ',1));\n']);
        fprintf(fid,['\n']);
        
        for j = 1:length(output)
            fprintf(fid,['y(' num2str(j) ') = ' char(output(j)) ';\n']);
        end
        
        fprintf(fid,['\n']);
        fprintf(fid,['%% SYSTEM STRUCT\n']);
        fprintf(fid,['\n']);
        if any([strcmp(expansion,'RRE'),strcmp(expansion,'EMRE'),strcmp(expansion,'LNA'),strcmp(expansion,'IOS')])
            fprintf(fid,['model.sym.nmx = ',num2str(nmx),';\n']);
        elseif(strcmp(expansion(1:3),'CME'))
            fprintf(fid,['model.sym.nmx = ',num2str(nmx),';\n']);
        else
            fprintf(fid,['model.sym.nmx = 0;\n']);
        end
        fprintf(fid,['model.sym.x = x;\n']);
        fprintf(fid,['model.sym.u = u;\n']);
        if (cvodewrap_flag)
            fprintf(fid,['model.sym.xdot = xdot;\n']);
        end
        if (idawrap_flag)
            fprintf(fid,['model.sym.f = f;\n']);
            fprintf(fid,['model.sym.M = M;\n']);
        end
        fprintf(fid,['model.sym.p = p;\n']);
        fprintf(fid,['model.sym.k = k;\n']);
        fprintf(fid,['model.sym.x0 = x0;\n']);
        if (idawrap_flag)
            fprintf(fid,['model.sym.dx0 = dx0;\n']);
        end
        fprintf(fid,['model.sym.y = y;\n']);
        fprintf(fid,['%% Additional fields for the prespecified length of kappa\n']);
        fprintf(fid,['model.sym.nk1 = ',num2str(kappa.nk1),';\n']);
        fprintf(fid,['end']);
        fclose(fid);
    end
    
    %% MODEL FILES D2D
    
    if(d2d_flag)
        model = expansion;
        clear([modeldir expansion '_' modelname '.def']);
        fid = fopen([modeldir expansion '_' modelname '.def'],'w');
        str_description = ['DESCRIPTION \n' ...
            '"' modelname '" \n\n'...
            'PREDICTOR \n'...
            't               T   min         time   0   80\n\n'...
            'COMPARTMENTS \n'];
        str_compartments = ['\n'];
        
        char_states = char(transpose(states));
        char_states = char_states(10:end-3);
        char_states = textscan(char_states,'%s','delimiter',',');
        char_states = char(char_states{1});
        str_states = [char_states repmat('        C   "nmol/l"      "conc."',length(states),1) ...
            repmat(' \n',length(states),1)];
        char_odes = char(transpose(xdot));
        char_odes = char_odes(10:end-3);
        char_odes = ['"' char_odes '"'];
        char_odes = strrep(char_odes,',','",');
        char_odes = strrep(char_odes,', ',', "');
        char_odes = textscan(char_odes,'%s','delimiter',',');
        char_odes = char(char_odes{1});
        str_odes = [char_odes  repmat(' \n',length(states),1)];
        
        % if ~isempty(ind_input)
        char_ic = char(transpose([System.state.mu0./System.state.volume;sym(zeros(size(states,1)-size(System.state.mu0,1),1))]));
        char_ic = char_ic(10:end-3);
        char_ic = textscan(char_ic,'%s','delimiter',',');
        char_ic = char(char_ic{1});
        str_conditions_1 = [repmat('init_',size(char_states,1),1) char_states ...
            repmat('       "',size(char_states,1),1) char_ic repmat('" \n',size(char_states,1),1)];
        
        fprintf(fid,str_description);
        fprintf(fid,str_compartments');
        fprintf(fid,'\nSTATES \n');
        fprintf(fid,str_states');
        %                 fprintf(fid,str_input);
        fprintf(fid,'\nODES \n');
        fprintf(fid,str_odes');
        fprintf(fid,'\nCONDITIONS \n');
        fprintf(fid,str_conditions_1');
        fclose(fid);
    end
    return
end
disp('Invalid Expansion String')
return

end

% better subs
function out = mysubs(in, old, new)
if(~isnumeric(in) && ~isempty(old) && ~isempty(findsym(in)))
    matVer = ver('MATLAB');
    if(str2double(matVer.Version)>=8.1)
        out = subs(in, old(:), new(:));
    else
        out = subs(in, old(:), new(:), 0);
    end
else
    out = in;
end
end