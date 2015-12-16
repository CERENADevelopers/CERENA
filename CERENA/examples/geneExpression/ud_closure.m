function [hoM_closure,flag,moment_closure] = ud_closure(I, order, System)
mu = sym('mu_',[System.state.number,1]);
hoM_closure = sym(zeros(size(I,1),1));
flag = zeros(size(I,1),1);
for i = 1:size(I,1)
    Ii = I(i,(find(I(i,:)~=0)));
    if Ii(2)==Ii(3)
        hoMi = sym(strrep(['C_',num2str(Ii)],'  ','_'));
        % generate all possible partitions in YI
        B = partitions(Ii);
        % generating cumulant corresponding to this moment
        K = sym(0);
        for b = 1:length(B)
            parts = [];
            for ib = 1:length(B{b})
                alpha_ib = convertI2alpha(B{b}{ib},System.state.number);
                uncentMom_ib = convertUncent2Cent(B{b}{ib},alpha_ib,mu,System.state.number);
                parts = [parts,uncentMom_ib];
            end
            npart = length(parts);
            K = K + factorial(npart-1) * (-1)^(npart-1) * prod(parts);
        end
        hoM_closure(i) = solve(K,hoMi);
        flag(i) = 1;
    else
        flag(i) = 0;
    end
end
moment_closure = 'low-dispersion';

