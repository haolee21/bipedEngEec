function [dObj,dObjGrad]=objFun(x,p)
q = x(1:p.numJ,:);
dq= x(p.numJ+1:2*p.numJ,:);
u = x(p.numJ*2+1:3*p.numJ);

dObj = 0.5*sum(u.^2);

dObjGrad = [zeros(p.numJ*2,size(x,2));ones(p.numJ,size(x,2))];

end