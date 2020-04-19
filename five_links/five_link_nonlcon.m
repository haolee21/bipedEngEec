function [c,ceq,gradc,gradceq] = five_link_nonlcon(x,p)

% x = gpuArray(x);
[c,gradc]=yposCon(x);
[ceq1,gradceq1]=dynConst(x,p);

[ceq2,gradceq2]=gaitLenCon(x,p);

[ceq3,gradceq3]=initYPosCons(x,p);
% [c,ceq1,gradc,gradceq1] = gather(c,ceq,gradc,gradceq);

ceq = [ceq1;ceq2*10;50*ceq3];
gradceq = [gradceq1,gradceq2*10,50*gradceq3];

end