function [c,ceq,gradc,gradceq] = discrete_nonlcon(x,p)

[ceq1,gradceq1] = dynConst_discrete(x,p);
[ceq2,gradceq2] = hipCon(x,p);

[c1,ceq3,gradc1,gradceq3]=grf_cons_d(x,p);

[ceq4,gradceq4]=initYPosCons(x,p);
[c2,gradc2]=yposCon(x,p);
c=[c1;c2];
gradc=[gradc1,gradc2];
ceq = [ceq1;ceq2;ceq3;ceq4];
gradceq = [gradceq1,gradceq2,gradceq3,gradceq4];

end