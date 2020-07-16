function [c,ceq,gradc,gradceq] = five_link_nonlcon(x,p)

% x = gpuArray(x);
[c4,gradc4]=yposCon(x,p);
% [c2,gradc2]=heelPosCon(x,p);
% [c2,gradc2]=headPosCon(x,p);
% [c3,gradc3]=ankPosCon(x,p);
% [c3,gradc3]=toe_heel(x,p);
[c1,gradc1]=hipVelCon(x,p);
%[c2,gradc2]=gaitLenCon(x,p);
[ceq1,gradceq1] = dynConst2(x,p);

[c2,~,gradc2,~] = grf_cons(x,p);

[ceq4,gradceq4]=hipCon(x,p);

[c3,gradc3]=slack_cons(x,p);

c = [c1;c2;c3;c4];
gradc=[gradc1,gradc2,gradc3,gradc4];

%[ceq1,gradceq1]=dynConst(x,p);

% [ceq3,gradceq3]=fri_cons(x,p);


%[ceq3,gradceq3]=initYPosCons(x,p);

% [c,ceq1,gradc,gradceq1] = gather(c,ceq,gradc,gradceq);
% [ceq5,gradceq5]=ext_const(x,p);
ceq = [ceq1;ceq4];
gradceq = [gradceq1,gradceq4];

end