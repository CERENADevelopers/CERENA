function i = findM(K,k)

i = find(all(bsxfun(@minus,K,k) == 0,2));