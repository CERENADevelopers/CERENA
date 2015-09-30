% getSSA.m
% constructs the system X(theta) of stochastic population given species and
% reactions.
%     
%
% NOTE: The propensity functions are assumed to be linear in theta!
%
% USAGE:
% ======
% [sys] = getFSP(reactions,species,reaction,xmin,xmax)
% [sys] = getFSP(x,theta,reaction,xmin,xmax,g)
%
% INPUT:
% ======
% species ... columnvector containing the symbolic expressions of the species.
% theta ... rowvector containing the symbolic expressions of the parameter.
% reaction ... struct containing the reaction details:
%   (i).propensity ... rate of the i.th propensity function (e.g. = k1*x1*x2).
%   (i).stoichiometry ... update of after i.th reaction occures.
%           This has to be a columnvector of the some dimension as the x
%           (e.g. = [0,1,0]').
%   (i).parameter ... parameter of i.th propensity function (e.g. = k1).
%
% OUTPUT:
% =======
% sys.propensities ... function handle providing propensity-vector v(theta)
% sys.stoichiometry ... stoichiometry matrix.
%
%       
% 28/10/2011 ... Sebastian Waider

% function [sys] = getSSA(species, theta, reaction)
function [sys] = getFSP(varargin)

%% CHECK / ASSIGN INPUTS
% Check required inputs
    if ((nargin~=3) || isempty(varargin{1}) || isempty(varargin{2}) || isempty(varargin{3}))
    % Error message:
    error('Three nonempty inputs are required');
    end
%% INITIALIZE ROUTINE
species = varargin{1};
theta= varargin{2};
reaction= varargin{3};


%% GENERATE SYSTEM

% Loop: propensity functions
for i = 1:length(reaction)
    % Find right probensity function:
    
    % Construct numerical function for rates:
    rate_wo_par = sym2fun(simplify(...
        reaction(i).propensity/reaction(i).parameter),species);
    % Initialize matrix:
    system.A_i{i} = sparse(zeros(nindex));
    % Construction of A^{(i)} column-wise:
    for j = 1:nindex
       system.A_i{i}(j,j) = - rate_wo_par(index(j,:)');
       
       j2 = (index(:,1) == (index(j,1)+system.reaction(i).stoichiometry(1)));
       for k = 2:nx
           j2 = j2 .* (index(:,k) == (index(j,k)+system.reaction(i).stoichiometry(k)));
       end
       j2 = find(j2);
       % Assign reaction
       if ~isempty(j2)
           system.A_i{i}(j2,j) = rate_wo_par(index(j,:)');
       end
    end
    % Assign coefficient:
    system.sigma_A_i{i} = eval(['@(theta) theta(' num2str(find(system.theta == system.reaction(i).parameter)) ');']);
end

%% CONSTRUCTION OF FUNCTION HANDLE 'A'
% Construct functional expression:
expr = '@(theta) [system.sigma_A_i{1}(theta)*system.A_i{1}';
for i = 2:length(system.theta)
    expr = [expr '+ system.sigma_A_i{' num2str(i) ' }(theta)*system.A_i{' num2str(i) '}'];
end
expr = [expr '];'];
% Construct function handle
system.A = eval(expr);

%% CONSTRUCTION OF B, C AND D MATRIX
% B Matrix:
system.B = @(theta) [];
% C Matrix:
system.C = @(theta) [];
% D Matrix:
system.D = @(theta) [];