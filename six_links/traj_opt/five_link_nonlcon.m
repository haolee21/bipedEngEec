function [c,ceq,gradc,gradceq] = five_link_nonlcon(x,p)

% x = gpuArray(x);
% [c1,gradc1]=yposCon(x,p);
% [c2,gradc2]=heelPosCon(x,p);
% [c2,gradc2]=headPosCon(x,p);
% [c3,gradc3]=ankPosCon(x,p);
% [c3,gradc3]=toe_heel(x,p);
[c1,gradc1]=hipVelCon(x,p);
%[c2,gradc2]=gaitLenCon(x,p);
c = c1;
gradc=gradc1;

%[ceq1,gradceq1]=dynConst(x,p);
[ceq1,gradceq1] = dynConst2(x,p);



%[ceq3,gradceq3]=initYPosCons(x,p);

[ceq4,gradceq4]=hipCon(x,p);
% [c,ceq1,gradc,gradceq1] = gather(c,ceq,gradc,gradceq);
% [ceq5,gradceq5]=ext_const(x,p);
ceq = [ceq1;ceq4];
gradceq = [gradceq1,gradceq4];

end