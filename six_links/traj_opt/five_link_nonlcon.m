function [c,ceq,gradc,gradceq] = five_link_nonlcon(x,p)

% x = gpuArray(x);
% [c1,gradc1]=yposCon(x,p);
% [c2,gradc2]=headPosCon(x,p);
% [c3,gradc3]=ankPosCon(x,p);
[c3,gradc3]=toe_heel(x,p);
c = c3;
gradc=gradc3;

[ceq1,gradceq1]=dynConst(x,p);

%[ceq2,gradceq2]=gaitLenCon(x,p);

%[ceq3,gradceq3]=initYPosCons(x,p);

[ceq4,gradceq4]=hipCon(x,p);
% [c,ceq1,gradc,gradceq1] = gather(c,ceq,gradc,gradceq);
[ceq5,gradceq5]=ext_const(x,p);
ceq = [ceq1;ceq4;ceq5];
gradceq = [gradceq1,gradceq4,gradceq5];

end