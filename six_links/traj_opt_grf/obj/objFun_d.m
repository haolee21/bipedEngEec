function [dObj,dObjGrad]=objFun_d(x,p)
% q = x(1:p.numJ,:);
% dq= x(p.numJ+1:2*p.numJ,:);
u = x(p.numJ+1:2*p.numJ,:);
%s_var = x(p.numJ*3+5:p.numJ*3+6,:);
dObj = 0.5*sum(p.jointW.*sum(u.^2,2).');
%dObj = dObj+0.5*sum(s_var.^2,'all');
dObjGrad = [zeros(p.numJ,size(x,2));diag(p.jointW)*u;zeros(6,size(x,2))];


end