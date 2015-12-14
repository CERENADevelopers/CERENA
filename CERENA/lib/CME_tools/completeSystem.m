%function System = completeSystem(System)
function System = completeSystem(varargin)
if nargin >= 1
    System = varargin{1};
else
    error('At least one input argument is required!')
end
options.sym_kappa = true;
if nargin >= 2
    options = setdefault(varargin{2},options);
end
tmpsize = size(System.compartments);
if tmpsize(1)==1
    System.compartments = transpose(System.compartments);
end
tmpsize = size(System.volumes);
if tmpsize(1)==1
    System.volumes = transpose(System.volumes);
end
%% STOICHIOMETRY
if isfield(System,'state')
    System.state.number   = length(System.state.variable);
    % converting to column vector
    tmpsize = size(System.state.variable);
    if tmpsize(1)==1
        System.state.variable = transpose(System.state.variable);
    end
    tmpsize = size(System.state.compartment);
    if tmpsize(1)==1
        System.state.compartment = transpose(System.state.compartment);
    end
    if isfield(System.state,'xmin')
        tmpsize = size(System.state.xmin);
        if tmpsize(1)==1
            System.state.xmin = transpose(System.state.xmin);
        end
    end
    if isfield(System.state,'xmax')
        tmpsize = size(System.state.xmax);
        if tmpsize(1)==1
            System.state.xmax = transpose(System.state.xmax);
        end
    end
    tmpsize = size(System.state.mu0);
    if tmpsize(1)==1
        System.state.mu0 = transpose(System.state.mu0);
    end
    %% Initial Conditions
    System.state.mu0 = sym(System.state.mu0);
    % Initial conditions for (co)variances
    if ~isfield(System.state,'C0')
        System.state.C0 = sym(zeros(System.state.number*(System.state.number+1)/2,1));
    else
        System.state.C0 = sym(System.state.C0);
        tmpsize = size(System.state.C0);
        if tmpsize(1)==1
            System.state.C0 = transpose(System.state.C0);
        end
    end
    if isfield(System,'kappa')
        kappa = System.kappa.variable;
    else
        kappa = [];
    end
    % Length of the kappa provided by the user
    nk1 = length(kappa);
    if options.sym_kappa
        % Adding indicators for initial conditions of mu0 to kappa
        kmu0_ind = sym('indmu',[length(System.state.mu0),1]);
        % Adding initial conditions of mu0 to kappa
        kmu0_sym = sym('kmu0',[length(System.state.mu0),1]);
        % Adding indicators for initial conditions of C0 to kappa
        kC0_ind = sym('indC',[length(System.state.C0),1]);
        % Adding initial conditions of C0 to kappa
        kC0_sym = sym('kC0',[length(System.state.C0),1]);
        % Adding all to kappa
        kappa = [kappa;kmu0_ind;kC0_ind;kmu0_sym;kC0_sym];
    end
    System.kappa.variable = kappa;
    System.kappa.nk1 = nk1;
    
    if options.sym_kappa
        % Symbolic array for the prespecified initial conditions of mu0 in the
        % modelDef
        fmu0_sym = sym('fmu0',[length(System.state.mu0),1]);
        % Symbolic array for the prespecified initial conditions of C0 in the
        % modelDef
        fC0_sym = sym('fC0',[length(System.state.C0),1]);
        
        % Saving the initial conditions prespecified in the modelDef
        System.state.fmu0 = System.state.mu0;
        System.state.fC0 = System.state.C0;
        System.state.fmu0_sym = fmu0_sym;
        System.state.fC0_sym = fC0_sym;
        
        System.state.mu0 = kmu0_ind .* kmu0_sym + (1 - kmu0_ind).* fmu0_sym;
        System.state.C0 = kC0_ind .* kC0_sym + (1 - kC0_ind).* fC0_sym;
    else
        System.state.fmu0 = [];
        System.state.fC0 = [];
        System.state.fmu0_sym = [];
        System.state.fC0_sym = [];    
    end
    
    %%%%%%% Other fields
    System.state.volume = sym(zeros(length(System.state.variable),1));
    for i = 1:System.state.number
        indComp = find(strcmp(System.compartments,System.state.compartment(i)));
        System.state.volume(i) = System.volumes(indComp);
    end
    % Default naming of species
    if ~isfield(System.state,'name')
        for i =1:System.state.number
            System.state.name{i,1} = char(System.state.variable(i));
        end
    end
    % Default type of species
    if ~isfield(System.state,'type')
        System.state.type = repmat({'moment'},[System.state.number,1]);
    end
end

if isfield(System,'parameter')
    tmpsize = size(System.parameter.variable);
    if tmpsize(1)==1
        System.parameter.variable = transpose(System.parameter.variable);
    end
    if ~isfield(System.parameter,'name')
        for i =1:length(System.parameter.variable)
            System.parameter.name{i,1} = char(System.parameter.variable(i));
        end
    else
        tmpsize = size(System.parameter.name);
        if tmpsize(1)==1
            System.parameter.name = transpose(System.parameter.name);
        end
    end
end

% Constant parameters
if isfield(System,'kappa')
    tmpsize = size(System.kappa.variable);
    if tmpsize(1)==1
        System.kappa.variable = transpose(System.kappa.variable);
    end
    if ~isfield(System.kappa,'name')
        for i =1:length(System.kappa.variable)
            System.kappa.name{i,1} = char(System.kappa.variable(i));
        end
    else
        tmpsize = size(System.kappa.name);
        if tmpsize(1)==1
            System.kappa.name = transpose(System.kappa.name);
        end
    end
end
% if isfield(system.parameter,'MacroscopicValue')
%     system.parameter.MicroscopicValue = system.parameter.MacroscopicValue;
% end
% Loop: reactions
if isfield(System,'reaction')
    S_e = zeros(length(System.state.variable),length(System.reaction));
    S_p = zeros(length(System.state.variable),length(System.reaction));
    for k = 1:length(System.reaction)
        % Educt
        for j = 1:length(System.reaction(k).educt)
            ind = find(System.state.variable == System.reaction(k).educt(j));
            S_e(ind,k) = S_e(ind,k) + 1;
        end
        
        % Product
        for j = 1:length(System.reaction(k).product)
            ind = find(System.state.variable == System.reaction(k).product(j));
            S_p(ind,k) = S_p(ind,k) + 1;
        end
        
        % Stoichiometry
        System.reaction(k).stoichiometry = S_p(:,k)-S_e(:,k);
        symvarReac = symvar(System.reaction(k).propensity);
        ind = ismember(symvarReac,System.parameter.variable);
        System.reaction(k).parameter = symvarReac(ind);
    end
    %% Derivation of rate for micro/macro propensities
    if isfield(System,'scaleIndicator')
        switch System.scaleIndicator
            case 'microscopic'
                for k = 1:length(System.reaction)
                    if any(S_e(:,k)>1)
                        tmpRate = coeffs(System.reaction(k).propensity,symvar(System.reaction(k).educt));
                        tmpInd = S_e(S_e(:,k)>1,k);
                        System.reaction(k).rate = tmpRate(end)*prod(factorial(tmpInd));
                    else
                        System.reaction(k).rate = coeffs(System.reaction(k).propensity,symvar(System.reaction(k).educt));
                    end
                end
            case 'macroscopic'
                for k = 1:length(System.reaction)
                    System.reaction(k).rate = coeffs(System.reaction(k).propensity,symvar(System.reaction(k).educt));
                end
        end
    end
    %% Microscopic and Macroscopic interconversion
%     if isfield(System,'scaleConversion')
%         switch System.scaleConversion
%             case 'Micro_to_Macro'
%                 for k = 1:length(System.reaction)
%                     if any(S_e(:,k)>1)
%                         System.reaction(k).propensity = System.reaction(k).rate * prod(System.reaction(k).educt);
%                     end
%                     indEduc = find(S_e(:,k)>0);
%                     tmpProp = System.reaction(k).propensity;
%                     tmpProp = subs(tmpProp,System.state.variable,System.state.variable.*System.state.volume);
%                     System.reaction(k).propensity = tmpProp/System.state.volume(indEduc(1));
%                     System.reaction(k).rate = coeffs(System.reaction(k).propensity,symvar(System.reaction(k).educt));
%                 end
%             case 'Macro_to_Micro'
%                 for k = 1:length(System.reaction)
%                     indEduc = find(S_e(:,k)>0);
%                     if any(S_e(:,k)>1)
%                         tmpProp = System.reaction(k).propensity;
%                         tmpProp = subs(tmpProp,System.state.variable,System.state.variable./System.state.volume);
%                         System.reaction(k).propensity = tmpProp * System.state.volume(indEduc(1));
%                         System.reaction(k).propensity = simplify(System.reaction(k).propensity);
%                         tmpProp = coeffs(System.reaction(k).propensity,symvar(System.reaction(k).educt));
%                         System.reaction(k).rate = tmpProp;
%                         %                         tmpProp = system.reaction(k).rate;
%                         for iEduc = indEduc
%                             tmpProp = tmpProp * expand(nchoosek(System.state.variable(iEduc),S_e(iEduc,k)));
%                         end
%                         System.reaction(k).propensity = tmpProp;
%                     else
%                         tmpProp = System.reaction(k).propensity;
%                         tmpProp = subs(tmpProp,System.state.variable,System.state.variable./System.state.volume);
%                         if ~isempty(indEduc)
%                             System.reaction(k).propensity = tmpProp * System.state.volume(indEduc(1));
%                         else
%                             indProd = find(S_p(:,k)>0);
%                             System.reaction(k).propensity = tmpProp * System.state.volume(indProd(1));
%                         end
%                         System.reaction(k).rate = coeffs(System.reaction(k).propensity,symvar(System.reaction(k).educt));
%                     end
%                 end
%         end
%     end
    %% More general - not assuming mass action on the reactants
    if isfield(System,'scaleConversion')
        switch System.scaleConversion
            case 'Micro_to_Macro'
                for k = 1:length(System.reaction)
                    if any(S_e(:,k)>1)
                        System.reaction(k).propensity = System.reaction(k).rate * prod(System.reaction(k).educt);
                    end
                    indEduc = find(S_e(:,k)>0);
                    tmpProp = System.reaction(k).propensity;
                    tmpProp = subs(tmpProp,System.state.variable,System.state.variable.*System.state.volume);
                    System.reaction(k).propensity = tmpProp/System.state.volume(indEduc(1));
                    System.reaction(k).rate = coeffs(System.reaction(k).propensity,symvar(System.reaction(k).educt));
                end
            case 'Macro_to_Micro'
                for k = 1:length(System.reaction)
                    tmpProp = System.reaction(k).propensity;
                    svtmpProp = symvar(tmpProp);
                    indsv = ismember(svtmpProp,System.state.variable);
                    svstate = svtmpProp(indsv);
                    if ~isempty(svstate)
                        tmpProp = subs(tmpProp,svstate.^2,svstate.*(svstate-1)/2);
                    end
                    indEduc = find(S_e(:,k)>0);
                    tmpProp = subs(tmpProp,System.state.variable,System.state.variable./System.state.volume);
                        if ~isempty(indEduc)
                            System.reaction(k).propensity = tmpProp * System.state.volume(indEduc(1));
                        else
                            indProd = find(S_p(:,k)>0);
                            System.reaction(k).propensity = tmpProp * System.state.volume(indProd(1));
                        end
                end
        end
    end
    
    System.eductStoichiometry = S_e;
    System.productStoichiometry = S_p;
    System.stoichiometry = S_p - S_e;
end

if isfield(System,'output')
    System.output.number   = length(System.output.variable);
    tmpsize = size(System.output.variable);
    if tmpsize(1)==1
        System.output.variable = transpose(System.output.variable);
    end
    tmpsize = size(System.output.function);
    if tmpsize(1)==1
        System.output.function = transpose(System.output.function);
    end
    if ~isfield(System.output,'name')
        for i =1:System.output.number
            System.output.name{i,1} = char(System.output.variable(i));
        end
    else
        tmpsize = size(System.output.name);
        if tmpsize(1)==1
            System.output.name = transpose(System.output.name);
        end
    end
end
if isfield(System,'input')
    System.input.number   = length(System.input.variable);
    tmpsize = size(System.input.variable);
    if tmpsize(1)==1
        System.input.variable = transpose(System.input.variable);
    end
    tmpsize = size(System.input.function);
    if tmpsize(1)<tmpsize(2)
        System.input.function = transpose(System.input.function);
    end
    if ~isfield(System.input,'name')
        for i =1:System.input.number
            System.input.name{i,1} = char(System.input.variable(i));
        end
    else
        tmpsize = size(System.input.name);
        if tmpsize(1)==1
            System.input.name = transpose(System.input.name);
        end
    end
end