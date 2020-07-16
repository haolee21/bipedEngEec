function [ceq,ceq_grad] = hipCon(x,p)
x1 = x(:,1);
xm = x(:,end);
ceq = hipPos_x(x1.')-hipPos_x(xm.')-p.hipLen;
if(~isnumeric(x))
    ceq_grad = sym(zeros(size(x,1),size(x,2)));
else
    
    ceq_grad = zeros(size(x,1),size(x,2));
end

ceq_grad(1:p.numJ,1)=dHipPos_x(x1.').';
ceq_grad(1:p.numJ,end)=-dHipPos_x(xm.').';
ceq_grad = reshape(ceq_grad,[size(x,1)*size(x,2),1]);
end