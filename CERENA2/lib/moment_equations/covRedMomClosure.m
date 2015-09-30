function hoC = covRedMomClosure(hoI, system)
hoC = sym(zeros(size(hoI,1),1));
for i = 1:size(hoI,1)
    hoI_ = hoI(i,:);
    hoI_ = hoI_(hoI_~=0);
    ind_hoI_ = find(ismember(system.reactionDistanceInd,hoI_,'rows'));
    dist_hoI_ = system.reactionDistance(ind_hoI_);
    i1 = hoI_(1);
    i2 = hoI_(2);
    if dist_hoI_>1
        [cost shortestPath] = dijkstra(system.reactionGraph, i1, i2);
        middleNode = shortestPath(ceil(end/2));
        if (middleNode == i1 || middleNode == i2)
            [cost shortestPath] = dijkstra(system.reactionGraph, i2, i1);
            middleNode = shortestPath(ceil(end/2));
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
        C3 = getC([middleNode,middleNode]);
        hoC(i) = C1 * C2 / (C3+1e-10);
    else
        hoC(i) = 0;
    end
end
