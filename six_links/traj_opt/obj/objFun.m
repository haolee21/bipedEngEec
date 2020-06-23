function [dObj,dObjGrad]=objFun(x,p)
% q = x(1:p.numJ,:);
% dq= x(p.numJ+1:2*p.numJ,:);
u = x(p.numJ*2+1:3*p.numJ,:);

dObj = 0.5*sum(p.jointW.*sum(u.^2,2).');



dObjGrad = [zeros(p.numJ*2,size(x,2));diag(p.jointW)*u];
% dObjGrad(2,:) = q(2,:);
% dObjGrad(5,:) = q(5,:);
% dObjGrad = dObjGrad+fri_grad*p.fri_coeff;
end