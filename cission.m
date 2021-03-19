function [SF] = cission(Sig,sb,sn);

SF = 0;
for i = 1 : 3
    for j = 1 : 3
        SF = SF + Sig(i,j)*0.5*(sb(i)*sn(j)+sb(j)*sn(i));
    end
end