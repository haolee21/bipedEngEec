function [c,ceq,gradc,gradceq] = five_link_nonlcon(x,p)

% x = gpuArray(x);
[c1,gradc1]=yposCon(x);
[c2,gradc2]=headPosCon(x,p);
c = [c1;c2];
gradc=[gradc1,gradc2];

[ceq1,gradceq1]=dynConst(x,p);

[ceq2,gradceq2]=gaitLenCon(x,p);

[ceq3,gradceq3]=initYPosCons(x,p);
% [c,ceq1,gradc,gradceq1] = gather(c,ceq,gradc,gradceq);

ceq = [10*ceq1;ceq2;ceq3];
gradceq = [10*gradceq1,gradceq2,gradceq3];

end