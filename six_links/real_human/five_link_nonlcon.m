function [c,ceq,gradc,gradceq] = five_link_nonlcon(x,p)

% x = gpuArray(x);
[c1,gradc1]=yposCon(x,p);
[c2,gradc2]=headPosCon(x,p);
[c3,gradc3]=ankPosCon(x,p);

c = [c1;c2;c3];
gradc=[gradc1,gradc2,gradc3];

[ceq1,gradceq1]=dynConst(x,p);

%[ceq2,gradceq2]=gaitLenCon(x,p);

[ceq3,gradceq3]=initYPosCons(x,p);

[ceq4,gradceq4]=hipCon(x,p);
% [c,ceq1,gradc,gradceq1] = gather(c,ceq,gradc,gradceq);

ceq = [ceq1;ceq3;ceq4];
gradceq = [gradceq1,gradceq3,gradceq4];

end