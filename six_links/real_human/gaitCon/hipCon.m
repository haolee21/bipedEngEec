function [ceq,ceq_grad] = hipCon(x,p)
x1 = x(:,1);
xm = x(:,end);
ceq = hip_x_pos(x1.')-hip_x_pos(xm.')-p.hipLen;

ceq_grad = zeros(size(x,1),size(x,2));

ceq_grad(1:3,1)=hip_x_grad(x1.').';
ceq_grad(1:3,end)=-hip_x_grad(xm.').';
ceq_grad = reshape(ceq_grad,[size(x,1)*size(x,2),1]);
end