% function getSSE(System,options)
function System = getSSE(varargin)
if nargin >= 1
    System = varargin{1};
else
    error('At least one input argument is required!')
end
options = [];
if nargin >= 2
    options = setdefault(varargin{2},options);
end

expansion = System.name;
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

% propensity vector
    %       ( k_j/Omega                zeroth order
    % f_j = ( k_j*\phi_v               first order
    %       ( k_j*Omega*\phi_v*\phi_w  second order
    R = size(System.reaction,2);
    N = System.state.number;
    % assemble states and parameters
    par = System.parameter.variable;
    kappa = System.kappa;
%     nk1 = System.kappa.nk1;
    phi = System.state.variable;
    s = System.eductStoichiometry;
    S = System.stoichiometry;
    
    disp(['assembling diffusion matrices ...'])
    
    a = sym(zeros(1,R));
    f0 = sym(zeros(1,R));
    f1 = sym(zeros(1,R));
    
    %     for j = 1:R
    %         a(j) = System.reaction(j).propensity;
    %         f0(j) = System.reaction(j).MacroscopicPropensity;
    %     end
    
    for j = 1:R
        a(j) = System.reaction(j).propensity;
        f0(j) = System.reaction(j).rate*prod(phi.^s(:,j));
        f1(j) = -System.reaction(j).rate/2*prod(phi.^(s(:,j)-1).*s(:,j).*(s(:,j)-1));
    end
    
    Di = sym(zeros(N,1));
    D1i = sym(zeros(N,1));
    for i = 1:N
        Di(i) = sum(S(i,:).*f0);
        D1i(i) = sum(S(i,:).*f1);
    end
    Di = simplify(Di);
    D1i = simplify(D1i);
    
    if(any([LNAFLAG,EMREFLAG,IOSFLAG]))
        Dij = sym(zeros(N,N));
        D1ij = sym(zeros(N,N));
        for i = 1:N
            for j = 1:N
                Dij(i,j) = sum(S(i,:).*S(j,:).*f0);
                D1ij(i,j) = sum(S(i,:).*S(j,:).*f1);
            end
        end
        Dij = simplify(Dij);
        D1ij = simplify(D1ij);
    end
    
    if(any([IOSFLAG]))
        Dijk = sym(zeros(N,N,N));
        for i = 1:N
            for j = 1:N
                for k = 1:N
                    Dijk(i,j,k) = sum(S(i,:).*S(j,:).*S(k,:).*f0);
                end
            end
        end
    end
    
    % assemble jacobians
    %
    % Jst...z_ij...r = d/dphi_s d/dphi_t ... d/dphi_z D_ij...r
    %
    disp(['assembling jacobian matrices ...'])
    
    if(any([LNAFLAG,EMREFLAG,IOSFLAG]))
        Ja_i = sym(zeros(N,N));
        J1a_i = sym(zeros(N,N));
        Ja_i = transpose(simplify(jacobian(Di,phi)));
        J1a_i = transpose(simplify(jacobian(D1i,phi)));
    end
    
    if(any([EMREFLAG,IOSFLAG]))
        Jab_i = sym(zeros(N,N,N));
        for i = 1:N
            Jab_i(:,:,i) = simplify(hessian(Di(i),phi));
            for j = 1:N
                Jabg_i(:,:,j,i) = diff(Jab_i(:,:,i),phi(j));
            end
        end
        
    end
    
    if(any([IOSFLAG]))
        Ja_ij = sym(zeros(N,N,N));
        for i = 1:N
            for j = 1:N
                Ja_ij(:,i,j) = simplify(jacobian(Dij(i,j),phi));
            end
        end
    end
    
    
    
    if(any([IOSFLAG]))
        Jab_ij = sym(zeros(N,N,N,N));
        for i = 1:N
            for j = 1:N
                Jab_ij(:,:,i,j) = simplify(hessian(Dij(i,j),phi));
            end
        end
    end
    
    %% EXPANSION
    % generate fluctuation states
    
    % LNA [e_r e_k]_0
    
    if(LNAFLAG)
        disp(['generating LNA states ...'])
%         LNA_states = sym(zeros(0,1));
        LNA_states = sym([]);
        I_LNA = [];
        
        for i = 1:N
            for j = 1:N
                eval(['syms COV' num2str(i) num2str(j) '_LNA']);
                if(i>=j)
                    eval(['LNA_states = [LNA_states COV' num2str(i) num2str(j) '_LNA];']);
                    I_LNA = [I_LNA;j,i];
                end
            end
        end
    end
    
    if(EMREFLAG)
        % EMRE [e_r]_1
        disp(['generating EMRE states ...'])
%         EMRE_states = sym(zeros(0,1));
        EMRE_states = sym([]);
        for i = 1:N
            eval(['syms ' char(phi(i)) '_CORR_EMRE']);
            eval(['EMRE_states = [EMRE_states ' char(phi(i)) '_CORR_EMRE];']);
        end
    end
    
    
    if(IOSFLAG)
        % IOS [e_r e_k]_2
        disp(['generating IOS states ...'])
%         IOS_states_1 = sym(zeros(0,1));
        IOS_states_1 = sym([]);
        for i = 1:N
            for j = 1:N
                eval(['syms COV' num2str(i) num2str(j) '_CORR_IOS']);
                if(i>=j)
                    eval(['IOS_states_1 = [IOS_states_1 COV' num2str(i) num2str(j) '_CORR_IOS];']);
                end
            end
        end
        
        % IOS [e_r e_k e_l]_1
%         IOS_states_2 = sym(zeros(0,1));
        IOS_states_2 = sym([]);
        I_IOS_2 = [];
        for i = 1:N
            for j = 1:N
                for k = 1:N
                    eval(['syms SKEW' num2str(i) num2str(j) num2str(k) '_IOS']);
                    if(all([j<=i,k<=j]))
                        eval(['IOS_states_2 = [IOS_states_2 SKEW' num2str(i) num2str(j) num2str(k) '_IOS];']);
                        I_IOS_2 = [I_IOS_2;k,j,i];
                    end
                end
            end
        end
    end
    
    
    % state numbers
    
    N_RRE = N;
    if(LNAFLAG)
        N_LNA = length(LNA_states);
    end
    if(EMREFLAG)
        N_EMRE = length(EMRE_states);
    end
    if(IOSFLAG)
        N_IOS_1 = length(IOS_states_1);
        N_IOS_2 = length(IOS_states_2);
    end
    % rhs of ode
    
    % RRE
    % dt phi
    disp(['assembling RRE equations ...'])
    xdot(1:N_RRE) = S*transpose(f0);
    
    % LNA
    % dt [e_r e_k]_0
    if(LNAFLAG)
        disp(['assembling LNA equations ...'])
        i_xdot = N_RRE + 1;
        for i = 1:N
            for j = 1:i
                xdot(i_xdot) = 0;
                for alpha = 1:N
                    eval(['xdot(i_xdot) = xdot(i_xdot) ' ...
                        '+ Ja_i(alpha,i)*COV' num2str(alpha) num2str(j) '_LNA ' ...
                        '+ Ja_i(alpha,j)*COV' num2str(alpha) num2str(i) '_LNA; ']);
                end
                eval(['xdot(i_xdot) = xdot(i_xdot) ' ...
                    '+ Dij(i,j);']);
                i_xdot = i_xdot+1;
            end
        end
    end
    
    % EMRE
    % dt [e_r]_1
    if(EMREFLAG)
        disp(['assembling EMRE equations ...'])
        for i = 1:N
            xdot(i_xdot) = 0;
            for alpha = 1:N
                eval(['xdot(i_xdot) = xdot(i_xdot) ' ...
                    '+ Ja_i(alpha,i)*' char(phi(alpha)) '_CORR_EMRE;']);
                for beta = 1:N
                    eval(['xdot(i_xdot) = xdot(i_xdot) + 1/2*Jab_i(alpha,beta,i)*COV' num2str(alpha) num2str(beta) '_LNA;']);
                end
            end
            eval(['xdot(i_xdot) = xdot(i_xdot) ' ...
                '+ D1i(i);']);
            i_xdot = i_xdot+1;
        end
    end
    
    % IOS 1
    % dt [e_r e_k]_2
    if(IOSFLAG)
        disp(['assembling IOS equations ...'])
        for i = 1:N
            for j = 1:i
                xdot(i_xdot) = 0;
                for alpha = 1:N
                    eval(['xdot(i_xdot) = xdot(i_xdot) ' ...
                        '+ Ja_i(alpha,i)*COV' num2str(alpha) num2str(j) '_CORR_IOS ' ...
                        '+ D1i(i)*' char(phi(j)) '_CORR_EMRE '...
                        '+ J1a_i(alpha,i)*COV' num2str(alpha) num2str(j) '_LNA' ...
                        '+ 1/2*Ja_ij(alpha,i,j)*' char(phi(alpha)) '_CORR_EMRE' ...
                        '+ Ja_i(alpha,j)*COV' num2str(alpha) num2str(i) '_CORR_IOS ' ...
                        '+ D1i(j)*' char(phi(i)) '_CORR_EMRE '...
                        '+ J1a_i(alpha,j)*COV' num2str(alpha) num2str(i) '_LNA' ...
                        '+ 1/2*Ja_ij(alpha,j,i)*' char(phi(alpha)) '_CORR_EMRE;' ...
                        ]);
                    for beta = 1:N
                        eval(['xdot(i_xdot) = xdot(i_xdot) ' ...
                            '+ 1/4*Jab_ij(alpha,beta,i,j)*COV' num2str(alpha) num2str(beta) '_LNA' ...
                            '+ 1/2*Jab_i(alpha,beta,i)*SKEW' num2str(alpha) num2str(beta) num2str(j) '_IOS' ...
                            '+ 1/4*Jab_ij(alpha,beta,j,i)*COV' num2str(alpha) num2str(beta) '_LNA' ...
                            '+ 1/2*Jab_i(alpha,beta,j)*SKEW' num2str(alpha) num2str(beta) num2str(i) '_IOS;' ...
                            ]);
                        for gamma = 1:N
                            eval(['xdot(i_xdot) = xdot(i_xdot) ' ...
                                '+ 1/6*Jabg_i(alpha,beta,gamma,i)*(COV' num2str(alpha) num2str(beta) '_LNA*COV' num2str(gamma) num2str(j) '_LNA ' ...
                                '+ COV' num2str(alpha) num2str(gamma) '_LNA*COV' num2str(beta) num2str(j) '_LNA ' ...
                                '+ COV' num2str(alpha) num2str(j) '_LNA*COV' num2str(beta) num2str(gamma) '_LNA)' ...
                                '+ 1/6*Jabg_i(alpha,beta,gamma,j)*(COV' num2str(alpha) num2str(beta) '_LNA*COV' num2str(gamma) num2str(i) '_LNA ' ...
                                '+ COV' num2str(alpha) num2str(gamma) '_LNA*COV' num2str(beta) num2str(i) '_LNA ' ...
                                '+ COV' num2str(alpha) num2str(i) '_LNA*COV' num2str(beta) num2str(gamma) '_LNA);' ...
                                ]);
                        end
                    end
                    
                end
                i_xdot = i_xdot+1;
            end
        end
        % IOS 2
        % dt [e_r e_k e_l]_1
        for i = 1:N
            for j = 1:i
                for k = 1:j
                    xdot(i_xdot) = 0;
                    for alpha = 1:N
                        eval(['xdot(i_xdot) = xdot(i_xdot) ' ...
                            '+ Ja_i(alpha,i)*SKEW' num2str(alpha) num2str(j) num2str(k) '_IOS ' ...
                            '+ Ja_ij(alpha,j,k)*COV' num2str(alpha) num2str(i) '_LNA '...
                            '+ Ja_i(alpha,j)*SKEW' num2str(alpha) num2str(i) num2str(k) '_IOS ' ...
                            '+ Ja_ij(alpha,i,k)*COV' num2str(alpha) num2str(j) '_LNA '...
                            '+ Ja_i(alpha,k)*SKEW' num2str(alpha) num2str(i) num2str(j) '_IOS ' ...
                            '+ Ja_ij(alpha,i,j)*COV' num2str(alpha) num2str(k) '_LNA; '...
                            ]);
                        for beta = 1:N
                            eval(['xdot(i_xdot) = xdot(i_xdot) ' ...
                                '+ 1/2*Jab_i(alpha,beta,i)*(COV' num2str(alpha) num2str(beta) '_LNA*COV' num2str(j) num2str(k) '_LNA ' ...
                                '+ COV' num2str(alpha) num2str(j) '_LNA*COV' num2str(beta) num2str(k) '_LNA ' ...
                                '+ COV' num2str(alpha) num2str(k) '_LNA*COV' num2str(beta) num2str(j) '_LNA) ' ...
                                '+ 1/2*Jab_i(alpha,beta,j)*(COV' num2str(alpha) num2str(beta) '_LNA*COV' num2str(i) num2str(k) '_LNA ' ...
                                '+ COV' num2str(alpha) num2str(i) '_LNA*COV' num2str(beta) num2str(k) '_LNA ' ...
                                '+ COV' num2str(alpha) num2str(k) '_LNA*COV' num2str(beta) num2str(i) '_LNA) ' ...
                                '+ 1/2*Jab_i(alpha,beta,k)*(COV' num2str(alpha) num2str(beta) '_LNA*COV' num2str(i) num2str(j) '_LNA ' ...
                                '+ COV' num2str(alpha) num2str(i) '_LNA*COV' num2str(beta) num2str(j) '_LNA ' ...
                                '+ COV' num2str(alpha) num2str(j) '_LNA*COV' num2str(beta) num2str(i) '_LNA); ' ...
                                ]);
                        end
                    end
                    eval(['xdot(i_xdot) = xdot(i_xdot)' ...
                        '+ 1/3*Dijk(i,j,k)' ...
                        '+ Dij(j,k)*' char(phi(i)) '_CORR_EMRE' ...
                        '+ D1i(i)*COV' num2str(j) num2str(k) '_LNA' ...
                        '+ 1/3*Dijk(j,i,k)' ...
                        '+ Dij(i,k)*' char(phi(j)) '_CORR_EMRE' ...
                        '+ D1i(j)*COV' num2str(i) num2str(k) '_LNA' ...
                        '+ 1/3*Dijk(k,i,j)' ...
                        '+ Dij(i,j)*' char(phi(k)) '_CORR_EMRE' ...
                        '+ D1i(k)*COV' num2str(i) num2str(j) '_LNA;' ...
                        ]);
                    i_xdot = i_xdot+1;
                    
                end
            end
        end
    end
    
    % exploit symmetry
    if(any([IOSFLAG,LNAFLAG]))
        for i = 1:N
            for j = i+1:N
                if(LNAFLAG)
                    eval(['xdot = subs(xdot,COV' num2str(i) num2str(j) '_LNA,COV' num2str(j) num2str(i) '_LNA);']);
                end
                if(IOSFLAG)
                    eval(['xdot = subs(xdot,COV' num2str(i) num2str(j) '_CORR_IOS,COV' num2str(j) num2str(i) '_CORR_IOS);']);
                end
            end
        end
    end
    
    if(any(IOSFLAG))
        for i = 1:N
            for j = 1:N
                for k = 1:N
                    ind = sort([i,j,k],2,'descend');
                    eval(['xdot = subs(xdot,SKEW' num2str(i) num2str(j) num2str(k) '_IOS,SKEW' num2str(ind(1)) num2str(ind(2)) num2str(ind(3)) '_IOS);']);
                end
            end
        end
    end
    
    % substitute input function
    if isfield(System,'input')
        xdot = subs(xdot,System.input.variable,System.input.function);
    end
    
    % substitute time
    syms t
    xdot = subs(xdot,System.time,t);
    
    % find System volume
    
    %         for j=1:System.state.number
    %             for k = 1:length(System.compartments)
    %                 if(strcmp(System.state.compartment(j),System.compartments(k)))
    %                     System.state.volume(j) = System.volumes(k);
    %                 end
    %             end
    %         end
    
    
    idx_covar = 1;
    for j = 1:System.state.number
        for k = 1:j
            covar_vol(idx_covar) = sqrt(System.state.volume(j))*sqrt(System.state.volume(k));
            idx_covar = idx_covar+1;
        end
    end
    
    idx_skew = 1;
    for j = 1:System.state.number
        for k = 1:j
            for l = 1:k
                skew_vol(idx_skew) = System.state.volume(j)*System.state.volume(k)*System.state.volume(l);
                idx_skew = idx_skew+1;
            end
        end
    end
    
    % transformation of variables
    for i = 1:N_RRE
        xdot(i) = xdot(i)/System.state.volume(i);
    end
    xdot = subs(xdot,phi,phi.*System.state.volume);
    
    
    % transformation of variables, we have to do transformation here
    % two times, once for change of variables from mol. numbers to
    % concentrations and once to account for the prefactor in the
    % System size expansion
    if(LNAFLAG)
        for i = N_RRE+1:N_RRE+N_LNA
            xdot(i) = xdot(i)/(covar_vol(i-N_RRE).^2);
        end
        xdot = subs(xdot,LNA_states,LNA_states.*(covar_vol.^2));
    end
    
    if(EMREFLAG)
        for i = N_RRE+N_LNA+1:N_RRE+N_LNA+N_EMRE
            xdot(i) = xdot(i)/System.state.volume(i-(N_RRE+N_LNA));
        end
        xdot = subs(xdot,EMRE_states,EMRE_states.*transpose(System.state.volume));
    end
    
    if(IOSFLAG)
        for i = N_RRE+N_LNA+N_EMRE+1:N_RRE+N_LNA+N_EMRE+N_IOS_1
            xdot(i) = xdot(i)/(covar_vol(i-(N_RRE+N_LNA+N_EMRE))^2);
        end
        xdot = subs(xdot,IOS_states_1,IOS_states_1.*(covar_vol.^2));
        
        
        
        for i = N_RRE+N_LNA+N_EMRE+N_IOS_1+1:N_RRE+N_LNA+N_EMRE+N_IOS_1+N_IOS_2
            xdot(i) = xdot(i)/skew_vol(i-(N_RRE+N_LNA+N_EMRE+N_IOS_1));
        end
        xdot = subs(xdot,IOS_states_2,IOS_states_2.*skew_vol);
    end
    
    xdot = simplify(xdot);
    
    switch(expansion)
        case 'RRE'
            states = phi;
        case 'LNA'
            states = [phi;transpose(LNA_states)];
        case 'EMRE'
            states = [phi;transpose(LNA_states);transpose(EMRE_states)];
        case 'IOS'
            states = [phi;transpose(LNA_states);transpose(EMRE_states);transpose(IOS_states_1);transpose(IOS_states_2)];
    end
    System.states = states;
    
    x0 = sym(zeros(size(states)));
    for j = 1:length(System.state.mu0)
        x0(j)= simplify(System.state.mu0(j)/System.state.volume(j));
    end
    %     % inital conditions for MATLAB-based simulation files
    %     x0_mat = sym(zeros(size(states)));
    %     for j = 1:length(System.state.fmu0)
    %         x0_mat(j)= simplify(System.state.fmu0(j)/System.state.volume(j));
    %     end
    %     % Finding the indices of paramteric initial conditions
    %     ind_mu0_sym = [];
    %     for j = 1:length(System.state.mu0)
    %         if ~isempty(symvar(System.state.mu0(j)))
    %             ind_mu0_sym = [ind_mu0_sym;j];
    %         end
    %     end
    %     % Adding initial conditions to kappa
    %     x0sym = sym('x0',size(states));
    %     % Length of the kappa provided by the user
    %     nk1 = length(kappa);
    %     kappa = [kappa;x0sym];
    
    if(any([EMREFLAG,IOSFLAG]))
        for j = 1:N_RRE
            output(j) =  phi(j) + EMRE_states(j);
        end
    else
        for j = 1:N_RRE
            output(j) =  phi(j);
        end
    end
    if(IOSFLAG)
        %         for j = 1:N_RRE
        %             output(N_RRE+j) =  LNA_states(j*(j+1)/2) + IOS_states_1(j*(j+1)/2)-EMRE_states(j)^2;
        %         end
        for j = 1:N_LNA
            output(N_RRE+j) =  LNA_states(j) + IOS_states_1(j) - EMRE_states(I_LNA(j,1)) * EMRE_states(I_LNA(j,2));
        end
        output = [output,IOS_states_2];
    elseif(LNAFLAG)
        %         for j = 1:N_RRE
        %             output(N_RRE+j) =  LNA_states(j*(j+1)/2) ;
        %         end
        for j = 1:N_LNA
            output(N_RRE+j) =  LNA_states(j) ;
        end
    end
    %     % Corrected covariances
    %     if (IOSFLAG)
    %         for j = 1:N_LNA
    %             cov_IOS(j) = LNA_states(j) + IOS_states_1(j) - EMRE_states(I_LNA(j,1)) * EMRE_states(I_LNA(j,2));
    %         end
    %         mom_IOS = [output(1:N),cov_IOS,IOS_states_2];
    %     end
    switch expansion
        case 'RRE'
            System.state.order = 1:System.state.number;
        case 'LNA'
            System.state.order = [[zeros(N_RRE,1),[1:System.state.number]'];I_LNA];
        case 'EMRE'
            System.state.order = [[zeros(N_RRE,1),[1:System.state.number]'];I_LNA];
        case 'IOS'
            I_LNA_ext = [zeros(N_LNA,1),I_LNA];
            System.state.order = [[zeros(N_RRE,2),[1:System.state.number]'];I_LNA_ext;I_IOS_2];
    end
    %% Output map for the output variables
    if isfield(System,'output')
        switch expansion
            case 'RRE'
                options.moment_order = 1;
                options.moment_order_output = 1;
                [H,My,Iy,n_Iy,M,I] = getOutputMoments(System,[],[],[],[],[],[],[],options);
                I_reorder = I;
                H = mysubs(H,M,states);
                for j = 1:length(H)
                    output(N_RRE+j) = H(j);
                end
            case 'LNA'
                options.moment_order = 2;
                options.moment_order_output = 2;
                [H,My,Iy,n_Iy,M,I] = getOutputMoments(System,[],[],[],[],[],[],[],options);
                [~,loc] = ismember(I_LNA,I,'rows');
                M_reorder = [M(1:N);M(loc)];
                I_reorder = [I(1:N,:);I(loc,:)];
                M_LNA = transpose(output);
                H = mysubs(H,M_reorder,M_LNA);
                for j = 1:length(H)
                    output(N_RRE+N_LNA+j) = H(j);
                end
            case 'EMRE'
                options.moment_order = 2;
                options.moment_order_output = 2;
                [H,My,Iy,n_Iy,M,I] = getOutputMoments(System,[],[],[],[],[],[],[],options);
                [~,loc] = ismember(I_LNA,I,'rows');
                M_reorder = [M(1:N);M(loc)];
                I_reorder = [I(1:N,:);I(loc,:)];
                M_EMRE = transpose(output);
                H = mysubs(H,M_reorder,M_EMRE);
                for j = 1:length(H)
                    output(N_RRE+N_LNA+j) = H(j);
                end
            case 'IOS'
                options.moment_order = 3;
                options.moment_order_output = 3;
                [H,My,Iy,n_Iy,M,I] = getOutputMoments(System,[],[],[],[],[],[],[],options);
                I_LNA_ext = [zeros(N_LNA,1),I_LNA];
                [~,loc_1] = ismember(I_LNA_ext,I,'rows');
                [~,loc_2] = ismember(I_IOS_2,I,'rows');
                M_reorder = [M(1:N);M(loc_1);M(loc_2)];
                I_reorder = [I(1:N,:);I(loc_1,:);I(loc_2,:)];
                M_IOS = transpose(output);
                H = mysubs(H,M_reorder,M_IOS);
                for j = 1:length(H)
                    output(N_RRE+N_LNA+N_IOS_2+j) = H(j);
                end
        end
        System.output.order = Iy;
%         System.state.order = I_reorder;
        c = sym('c',[length(par),1]);
        output = subs(output,c,par);
        output = subs(output,System.time,t);
    end
    
    %% Assemble output argument
    if strcmp(expansion,'RRE')
        System.SSE.nmx = N_RRE;
    elseif any([strcmp(expansion,'EMRE'),strcmp(expansion,'LNA')])
        System.SSE.nmx = N_RRE+N_LNA;
    elseif strcmp(expansion,'IOS')
        System.SSE.nmx = N_RRE+N_LNA+N_IOS_2;
    end
    System.SSE.states = states;
    System.SSE.xdot = xdot;
    System.SSE.x0 = x0;
    System.SSE.output = output;
    
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