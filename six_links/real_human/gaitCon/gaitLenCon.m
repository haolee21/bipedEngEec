function [ceq,ceq_grad] = gaitLenCon(x,p)
x1 = x(:,1);
xm = x(:,end);
ceq = end_x_pos(xm(1),xm(2),xm(3),xm(4),xm(5),xm(6))-end_x_pos(x1(1),x1(2),x1(3),x1(4),x1(5),x1(6))+p.gaitLen;

ceq_grad = zeros(size(x,1),size(x,2));

ceq_grad(1:p.numJ,1)=-1*end_x_grad(x1(1),x1(2),x1(3),x1(4),x1(5),x1(6)).';
ceq_grad(1:p.numJ,end)=end_x_grad(xm(1),xm(2),xm(3),xm(4),xm(5),xm(6)).';


ceq_grad = reshape(ceq_grad,[size(x,1)*size(x,2),1]);


end