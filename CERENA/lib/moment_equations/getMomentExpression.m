function expr = getMomentExpression(state_names,order,options)
    
    switch options.notation
        % Problem specific formulation
        case 'specific'
            % 0th order reaction
            if sum(order) == 0
                expr = '';
            % no moment closure required
            elseif sum(order) <= options.moment_order
                expr = '\\mathbb{E}[';
                lind = find(order>=1);
                for l = lind(:)'
                    if l > lind(1)
                        expr = [expr '\\'];
                    end
                    expr = [expr state_names{l}];
                    if order(l) >= 2
                        expr = [expr '^' num2str(order(l),'%d')];
                    end
                end
                expr = [expr ']'];
            % moment closure required
            elseif sum(order) > options.moment_order
                switch options.moment_closure
                    case 'none'
                        expr = '\\mathbb{E}[';
                        lind = find(order>=1);
                        for l = lind(:)'
                            if l > lind(1)
                                expr = [expr '\\'];
                            end
                            expr = [expr state_names{l}];
                            if order(l) >= 2
                                expr = [expr '^' num2str(order(l),'%d')];
                            end
                        end
                        expr = [expr ']'];
                    otherwise
                        error('This is not supported so far.')
                end
            end
            
        % Problem independent notation
        case 'general'
            % 0th order reaction
            if sum(order) == 0
                expr = '';
            % no moment closure required
            elseif sum(order) <= options.moment_order
                expr = 'M_{';
                lind = find(order>=1);
                for l = lind(:)'
                    for l2 = 1:order(l)
                        expr = [expr num2str(l,'%d')];
                    end
                end
                expr = [expr '}'];
            % moment closure required
            elseif sum(order) > options.moment_order
                switch options.moment_closure
                    case 'none'
                        expr = 'M_{';
                        lind = find(order>=1);
                        for l = lind(:)'
                            for l2 = 1:order(l)
                                expr = [expr num2str(l,'%d')];
                            end
                        end
                        expr = [expr '}'];
                    otherwise
                        error('This is not supported so far.')
                end
            end
            
    end    
end