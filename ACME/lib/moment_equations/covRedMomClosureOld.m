function hoC = covRedMomClosure(hoI, system)
hoC = sym(zeros(size(hoI,1),1));
for i = 1:size(hoI,1)
    hoI_ = hoI(i,:);
    hoI_ = hoI_(hoI_~=0);
    ind_hoI_ = find(ismember(system.reactionDistanceInd,hoI_,'rows'));
    dist_hoI_ = system.reactionDistance(ind_hoI_);
    i1 = hoI_(1);
    i2 = hoI_(2);
    ind_middlePair = [find(system.reactionDistanceInd(:,1) == i1);find(system.reactionDistanceInd(:,2) == i1)];
    ind_middlePair = setdiff(ind_middlePair,ind_hoI_);
    dist1 = 100;
    dist2 = 100;
    while((dist1>dist_hoI_ || dist2 > dist_hoI_) && ~isempty(ind_middlePair))
        [~,imin] = min(system.reactionDistance(ind_middlePair));
        middleNode = setdiff(system.reactionDistanceInd(ind_middlePair(imin),:),i1);
        if middleNode < i2
            tmp = find(ismember(system.reactionDistanceInd,[middleNode,i2],'rows'));
            tmp2 = find(ismember(hoI,[0,0,middleNode,i2],'rows'));
        else
            tmp = find(ismember(system.reactionDistanceInd,[i2,middleNode],'rows'));
            tmp2 = find(ismember(hoI,[0,0,i2,middleNode],'rows'));
        end
        if isempty(tmp2)
            dist1 = system.reactionDistance(ind_middlePair(imin));
            dist2 = system.reactionDistance(tmp);
        else
            dist1 = 100;
            dist2 = 100;
        end
        ind_middlePair = setdiff(ind_middlePair,ind_middlePair(imin));
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
end
