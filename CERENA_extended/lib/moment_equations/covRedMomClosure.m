function hoC = covRedMomClosure(hoI, system)
if ~isnumeric(hoI)
    tmpI = textscan(char(hoI),'%s %d %d','Delimiter','_');
    hoI = [tmpI{2},tmpI{3}];
end
hoC = sym(zeros(size(hoI,1),1));
for i = 1:size(hoI,1)
    hoI_ = hoI(i,:);
    hoI_ = hoI_(hoI_~=0);
%     ind_hoI_ = find(ismember(system.reactionDistanceInd,hoI_,'rows'));
%     dist_hoI_ = system.reactionDistance(ind_hoI_);
    i1 = hoI_(1);
    i2 = hoI_(2);
    if system.reactionGraph(hoI_(1),hoI_(2))==0 %if dist_hoI_>1
        [cost,shortestPath] = dijkstra(system.reactionGraph, i1, i2);
        middleNode = shortestPath(ceil(end/2));
        if (middleNode == i1 || middleNode == i2)
            [cost,shortestPath] = dijkstra(system.reactionGraph, i2, i1);
            middleNode = shortestPath(ceil(end/2));
            if (middleNode == i1 || middleNode == i2) %%%CHANGE LATER!!!
                disp('middle node set to one!')
                middleNode = 1;
            end
        end
        if middleNode < i1
            C1 = getC([middleNode,i1]);
        else
            C1 = getC([i1,middleNode]);
        end
        
        if middleNode < i2
            C2 = getC([middleNode,i2]);
        else
            C2 = getC([i2,middleNode]);
        end
%         C3 = getC([middleNode,middleNode]);
        mu_1 = getMu(i1);
        mu_2 = getMu(i2);
        mu_3 = getMu(middleNode);
        par_C1 = system.reactionGraphParameter(middleNode,i2);
        par_C2 = system.reactionGraphParameter(middleNode,i1);
%         hoC(i) = C1 * C2 / (C3+1e-10);
        hoC(i) = -(par_C1*mu_2*C1 + par_C2*mu_1*C2) / ((par_C1 + par_C2)*(mu_3+1e-10));
    else
        hoC(i) = 0;
    end
end
