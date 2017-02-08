%function System = completeSystem(System)
function System = completeSystemGraph(varargin)
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
    
    System.eductStoichiometry = S_e;
    System.productStoichiometry = S_p;
    System.stoichiometry = S_p - S_e;
end

