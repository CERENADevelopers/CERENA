function system = completeSystemFSP(system)

tmpsize = size(system.compartments);
if tmpsize(1)==1
    system.compartments = transpose(system.compartments);
end
tmpsize = size(system.volumes);
if tmpsize(1)==1
    system.volumes = transpose(system.volumes);
end
%% STOICHIOMETRY
if isfield(system,'state')
    system.state.number   = length(system.state.variable);
    % converting to column vector
    tmpsize = size(system.state.variable);
    if tmpsize(1)==1
        system.state.variable = transpose(system.state.variable);
    end
    tmpsize = size(system.state.compartment);
    if tmpsize(1)==1
        system.state.compartment = transpose(system.state.compartment);
    end
    if isfield(system.state,'xmin')
        tmpsize = size(system.state.xmin);
        if tmpsize(1)==1
            system.state.xmin = transpose(system.state.xmin);
        end
    end
    if isfield(system.state,'xmax')
        tmpsize = size(system.state.xmax);
        if tmpsize(1)==1
            system.state.xmax = transpose(system.state.xmax);
        end
    end
    tmpsize = size(system.state.mu0);
    if tmpsize(1)==1
        system.state.mu0 = transpose(system.state.mu0);
    end
    %% Initial Conditions
    % Initial conditions for (co)variances
    if ~isfield(system.state,'C0')
        system.state.C0 = sym(zeros(system.state.number*(system.state.number+1)/2,1));
    else
        tmpsize = size(system.state.C0);
        if tmpsize(1)==1
            system.state.C0 = transpose(system.state.C0);
        end
    end
    
    %%%%%%% Other fields
    system.state.volume = sym(zeros(length(system.state.variable),1));
    for i = 1:system.state.number
        indComp = find(strcmp(system.compartments,system.state.compartment(i)));
        system.state.volume(i) = system.volumes(indComp);
    end
    % Default naming of species
    if ~isfield(system.state,'name')
        for i =1:system.state.number
            system.state.name{i,1} = char(system.state.variable(i));
        end
    end
    % Default type of species
    if ~isfield(system.state,'type')
        system.state.type = repmat({'moment'},[system.state.number,1]);
    end
end

if isfield(system,'parameter')
    tmpsize = size(system.parameter.variable);
    if tmpsize(1)==1
        system.parameter.variable = transpose(system.parameter.variable);
    end
    if ~isfield(system.parameter,'name')
        for i =1:length(system.parameter.variable)
            system.parameter.name{i,1} = char(system.parameter.variable(i));
        end
    else
        tmpsize = size(system.parameter.name);
        if tmpsize(1)==1
            system.parameter.name = transpose(system.parameter.name);
        end
    end
end

% Constant parameters
if isfield(system,'kappa')
    tmpsize = size(system.kappa.variable);
    if tmpsize(1)==1
        system.kappa.variable = transpose(system.kappa.variable);
    end
    if ~isfield(system.kappa,'name')
        for i =1:length(system.kappa.variable)
            system.kappa.name{i,1} = char(system.kappa.variable(i));
        end
    else
        tmpsize = size(system.kappa.name);
        if tmpsize(1)==1
            system.kappa.name = transpose(system.kappa.name);
        end
    end
    % Adding kappa to parameters as for matlab-based kappa and parameters
    % are the same
    system.parameter.variable = [system.parameter.variable;system.kappa.variable];
    system.parameter.name = [system.parameter.name;system.kappa.name];
end

% Loop: reactions
if isfield(system,'reaction')
    S_e = zeros(length(system.state.variable),length(system.reaction));
    S_p = zeros(length(system.state.variable),length(system.reaction));
    for k = 1:length(system.reaction)
        % Educt
        
        for j = 1:length(system.reaction(k).educt)
            ind = find(system.state.variable == system.reaction(k).educt(j));
            S_e(ind,k) = S_e(ind,k) + 1;
        end
        
        % Product
        
        for j = 1:length(system.reaction(k).product)
            ind = find(system.state.variable == system.reaction(k).product(j));
            S_p(ind,k) = S_p(ind,k) + 1;
        end
        
        % Stoichiometry
        system.reaction(k).stoichiometry = S_p(:,k)-S_e(:,k);
        symvarReac = symvar(system.reaction(k).propensity);
        ind = ismember(symvarReac,system.parameter.variable);
        system.reaction(k).parameter = symvarReac(ind);
    end
     %% Derivation of rate frommicro/macro propensities
    if isfield(system,'scaleIndicator')
        switch system.scaleIndicator
            case 'microscopic'
                for k = 1:length(system.reaction)
                    if any(S_e(:,k)>1)
                        tmpRate = coeffs(system.reaction(k).propensity,symvar(system.reaction(k).educt));
                        tmpInd = S_e(S_e(:,k)>1,k);
                        system.reaction(k).rate = tmpRate(end)*prod(factorial(tmpInd));
                    else
                        system.reaction(k).rate = coeffs(system.reaction(k).propensity,symvar(system.reaction(k).educt));
                    end
                end
            case 'macroscopic'
                for k = 1:length(system.reaction)
                    system.reaction(k).rate = coeffs(system.reaction(k).propensity,symvar(system.reaction(k).educt));
                end
        end
    end
    %% Microscopic and Macroscopic interconversion
    if isfield(system,'scaleConversion')
        switch system.scaleConversion
            case 'Micro_to_Macro'
                for k = 1:length(system.reaction)
                    if any(S_e(:,k)>1)
                        system.reaction(k).propensity = system.reaction(k).rate * prod(system.reaction(k).educt);
                    end
                    indEduc = find(S_e(:,k)>0);
                    tmpProp = system.reaction(k).propensity;
                    tmpProp = subs(tmpProp,system.state.variable,system.state.variable.*system.state.volume);
                    system.reaction(k).propensity = tmpProp/system.state.volume(indEduc(1));
                end
            case 'Macro_to_Micro'
                for k = 1:length(system.reaction)
                    indEduc = find(S_e(:,k)>0);
                    if any(S_e(:,k)>1)
                        tmpProp = system.reaction(k).rate;
                        for iEduc = indEduc
                            tmpProp = tmpProp * expand(nchoosek(system.state.variable(iEduc),S_e(iEduc,k)));
                        end
                        system.reaction(k).propensity = tmpProp;
                    end
                    tmpProp = system.reaction(k).propensity;
                    tmpProp = subs(tmpProp,system.state.variable,system.state.variable./system.state.volume);
                    system.reaction(k).propensity = tmpProp * system.state.volume(indEduc(1));
                end
        end
    end
    system.eductStoichiometry = S_e;
    system.productStoichiometry = S_p;
    system.stoichiometry = S_p - S_e;
end

if isfield(system,'output')
    system.output.number   = length(system.output.variable);
    tmpsize = size(system.output.variable);
    if tmpsize(1)==1
        system.output.variable = transpose(system.output.variable);
    end
    tmpsize = size(system.output.function);
    if tmpsize(1)==1
        system.output.function = transpose(system.output.function);
    end
    if ~isfield(system.output,'name')
        for i =1:system.output.number
            system.output.name{i,1} = char(system.output.variable(i));
        end
    else
        tmpsize = size(system.output.name);
        if tmpsize(1)==1
            system.output.name = transpose(system.output.name);
        end
    end
end
if isfield(system,'input')
    system.input.number   = length(system.input.variable);
    tmpsize = size(system.input.variable);
    if tmpsize(1)==1
        system.input.variable = transpose(system.input.variable);
    end
    tmpsize = size(system.input.function);
    if tmpsize(1)==1
        system.input.function = transpose(system.input.function);
    end
    if ~isfield(system.input,'name')
        for i =1:system.input.number
            system.input.name{i,1} = char(system.input.variable(i));
        end
    else
        tmpsize = size(system.input.name);
        if tmpsize(1)==1
            system.input.name = transpose(system.input.name);
        end
    end
end
%% FSP
[system] = getFSP(system,system.state.xmin,system.state.xmax);
% Initial condition
system.x0 = zeros(size(system.index,1),1);
[~,ix] =ismember(system.state.mu0',system.index,'Rows');
system.x0(ix) = 1;
