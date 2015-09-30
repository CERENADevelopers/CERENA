% testGradient.m calculates finite difference approximations to the
%   gradient to check an analytical version.
%
%   backward differences: g_fd_f = (f(theta+eps*e_i) - f(theta))/eps
%   forward differences:  g_fd_b = (f(theta) - f(theta-eps*e_i))/eps
%   central differences:  g_fd_c = (f(theta+eps*e_i) - f(theta-eps*e_i))/(2*eps)
%
%   in order to work with tensors of order n the gradient must be returned as tensor of
%   order n+1 where the n+1th tensor dimension indexes the parameters with respect to which
%   the differentiation was carried out
%
% USAGE:
% ======
% [...] = testGradient(theta,fun,eps,il,ig)
% [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(...)
%
% INPUTS:
% =======
% theta ... parameter vector at which gradient is evaluated
% fun ... function of theta for which gradients are checked
% eps ... epsilon used for finite difference approximation of gradient
% il ... argout index at which function values are returned
% ig ... argout index at which gradient values are returned
%
% OUTPUTS:
% ========
% g ... gradient computed by f
% g_fd_f ... backward differences
% g_fd_b ... forward differences
% g_fd_c ... central differences
%
% 2014/06/11 Jan Hasenauer
% 2015/01/16 Fabian Froehlich

%function [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(theta,fun,eps,il,ig)
function [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(varargin)
    
    theta = varargin{1};
    fun = varargin{2};
    eps = varargin{3};
    if nargin >= 4
        il = varargin{4};
    else
        il = 1;
    end
    
    if nargin >= 5
        ig = varargin{5};
    else
        ig = 2;
    end
    theta = theta(:);
    
    
    str_1 = '[';
    % Evaluation of function and gradient
    i = 1;
    while true
        if(i==il)
            str_1 = [str_1 'l'];
        elseif(i==ig)
            str_1 = [str_1 'g'];
        else
            str_1 = [str_1 '~'];
        end
        if i == max(il,ig);
            break;
        else
            str_1 = [str_1 ','];
        end
        i=i+1;
        
    end
    
    eval([str_1 '] = fun(theta);']);
    
    % Computation of finite difference gradient
    
    
    
    g_fd_f = nan(size(g));
    g_fd_b = nan(size(g));
    g_fd_c = nan(size(g));
    
    
    str_2 = '[';
    % Evaluation of function and gradient
    i = 1;
    while true
        if(i==il)
            str_2 = [str_2 'l'];
        else
            str_2 = [str_2 '~'];
        end
        if i == max(il);
            break;
        else
            str_2 = [str_2 ','];
        end
        i=i+1;
        
    end
    
    for i = 1:length(theta)
        % function evaluation
        [zeros(i-1,1);eps;zeros(length(theta)-i,1)];
        size(theta);
        size([zeros(i-1,1);eps;zeros(length(theta)-i,1)]);
        
        eval([str_2 '_i_f] = fun(theta+[zeros(i-1,1);eps;zeros(length(theta)-i,1)]);']);
        eval([str_2 '_i_b] = fun(theta-[zeros(i-1,1);eps;zeros(length(theta)-i,1)]);']);
        
        % forward differences
        sg = size(g);
        if(length(theta)==1)
            eval(['g_fd_f(' repmat(':,',1,numel(size(g))) 'i) = (l_i_f-l)/eps;'])
            
            % backward differences
            eval(['g_fd_b(' repmat(':,',1,numel(size(g))) 'i) = -(l_i_b-l)/eps;'])
            
            % central differences
            eval(['g_fd_c(' repmat(':,',1,numel(size(g))) 'i) = (l_i_f-l_i_b)/(2*eps);'])
        elseif(sg(end)==1)
            eval(['g_fd_f(' repmat(':,',1,numel(size(g))-2) 'i) = (l_i_f-l)/eps;'])
            
            % backward differences
            eval(['g_fd_b(' repmat(':,',1,numel(size(g))-2) 'i) = -(l_i_b-l)/eps;'])
            
            % central differences
            eval(['g_fd_c(' repmat(':,',1,numel(size(g))-2) 'i) = (l_i_f-l_i_b)/(2*eps);'])            
        else
            eval(['g_fd_f(' repmat(':,',1,numel(size(g))-1) 'i) = (l_i_f-l)/eps;'])
            
            % backward differences
            eval(['g_fd_b(' repmat(':,',1,numel(size(g))-1) 'i) = -(l_i_b-l)/eps;'])
            
            % central differences
            eval(['g_fd_c(' repmat(':,',1,numel(size(g))-1) 'i) = (l_i_f-l_i_b)/(2*eps);'])
        end
    end
    
    % Output
    [g,g_fd_f,g_fd_b,g_fd_c];
    
