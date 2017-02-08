% function System = getReactionGraph(System,options)
function System = getReactionGraph(varargin)
if nargin>=1
    System = varargin{1};
else
    error('at least one input is required!')
end
options.reduction_order = 1;
if nargin >=2
    options = setdefault(varargin{2},options);
end
if options.moment_order == 3
    I_o3 = [];
end
rGraph = eye(System.state.number);
rGraph_o2 = zeros(System.state.number);
rGraph_o3 = zeros(System.state.number);
% rGraphParam = sym(zeros(System.state.number));
for ir= 1:length(System.reaction)
    ie = System.eductStoichiometry(:,ir);
    ie = find(ie~=0);
    ip = System.productStoichiometry(:,ir);
    ip = find(ip~=0);
    svProp = symvar(System.reaction(ir).propensity);
    if ~isempty(svProp)
        eductAll = intersect(svProp,System.state.variable);
        ieAll = find(ismember(System.state.variable,eductAll));
        ieProp = setdiff(ieAll,ie);
        ie = unique([ie;ieProp]);
        ip = unique([ip;ieProp]);
    end
    % edge between educts and products
    for j = 1:length(ie)
        for k = 1:length(ip)
            rGraph(ie(j),ip(k)) = 1;
            rGraph(ip(k),ie(j)) = 1;
            
%             rGraphParam(ie(j),ip(k)) = System.reaction(ir).parameter(1);
%             rGraphParam(ip(k),ie(j)) = System.reaction(ir).parameter(1);
        end
    end
    % edge between educts
    if length(ie) >= 2
        rEdges = nchoosek(ie,2);
        for j = 1:size(rEdges,1)
            rGraph(rEdges(j,1),rEdges(j,2)) = 1;
            rGraph(rEdges(j,2),rEdges(j,1)) = 1;
            
%             rGraphParam(rEdges(j,1),rEdges(j,2)) = System.reaction(ir).parameter(1);
%             rGraphParam(rEdges(j,1),rEdges(j,2)) = System.reaction(ir).parameter(1);
        end
    end
    
    % 3rd-order moments
    if options.moment_order == 3
        %                     v = ie;
        %                     M = ip;
        M1 = [];
        for i_ip = 1:length(ip)
            
            v = [ie;ip(i_ip)];
            M = [ie;ip(i_ip)];
            for j = 1:2
                M = [kron(M,ones(length(v),1)),kron(ones(size(M,1),1),v)];
            end
            M1 = [M1;M];
        end
        M1 = unique(sort(M1,2),'rows');
        
        v = ie;
        M = ie;
        for j = 1:2
            M = [kron(M,ones(length(v),1)),kron(ones(size(M,1),1),v)];
        end
        M2 = unique(sort(M,2),'rows');
        if ~isempty(M1)
            if ~isempty(M2)
                I_o3 = [I_o3;union(M1,M2,'rows')];
            else
                I_o3 = [I_o3;M1];
            end
        else
            if ~isempty(M2)
                I_o3 = [I_o3;M2];
            end
        end
        %         if length(ie) > 1
        %             twoEducts = nchoosek(ie,2);
        %             ie_ip = [ie;ip];
        %             for j=1:size(twoEducts,1)
        %                 tmp = [repmat(twoEducts(j,:),length(ie_ip),1),ie_ip];
        %                 tmp = unique(sort(tmp,2),'rows');
        %                 I_o3 = [I_o3;tmp];
        %             end
        %         end
    end
end

if options.reduction_order >= 2
    
    for i = 1:size(rGraph,1)
        tmpind = find(rGraph(i,:));
        tmpind = tmpind(tmpind~=i);
        tmpind = tmpind';
        tmpind = [kron(tmpind,ones(length(tmpind),1)),kron(ones(length(tmpind),1),tmpind)];
        tmpind = unique(sort(tmpind,2),'rows');
        if ~isempty(tmpind)
            rGraph_o2((tmpind(:,2)-1)*size(rGraph_o2,1)+tmpind(:,1)) = 1 - rGraph((tmpind(:,2)-1)*size(rGraph,1)+tmpind(:,1));
        end
        %         rGraph((tmpind(:,2)-1)*size(rGraph,1)+tmpind(:,1)) = 1;
    end
end

if options.reduction_order >= 3
    
    for i = 1:size(rGraph,1)
        tmpind1 = find(rGraph(i,:));
        tmpind2 = find(rGraph_o2(i,:));
        tmpind1 = tmpind1(tmpind1~=i);
        tmpind2 = tmpind2(tmpind2~=i);
        tmpind1 = tmpind1';
        tmpind2 = tmpind2';
        tmpind = [kron(tmpind1,ones(length(tmpind2),1)),kron(ones(length(tmpind1),1),tmpind2)];
        tmpind = unique(sort(tmpind,2),'rows');
        if ~isempty(tmpind)
            rGraph_o3((tmpind(:,2)-1)*size(rGraph,1)+tmpind(:,1)) = 1 - rGraph_o2((tmpind(:,2)-1)*size(rGraph,1)+tmpind(:,1)) - rGraph((tmpind(:,2)-1)*size(rGraph,1)+tmpind(:,1));
        end
    end
end

rGraph = rGraph + rGraph_o2 + rGraph_o3;
System.reactionGraph = rGraph;
% System.reactionGraphParameter = rGraphParam;
if options.moment_order == 3
    System.I_o3 = unique(I_o3,'rows');
end