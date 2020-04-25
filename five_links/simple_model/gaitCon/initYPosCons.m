function [ceq,ceq_grad]=initYPosCons(x,p)
x1=x(:,1).';
xm=x(:,end).';

initYPos=end_y_pos(x1);
endYPos=end_y_pos(xm);

ceq=[p.init_y-initYPos;p.init_y-endYPos];


ceq_grad1 = zeros(size(x,1),size(x,2));
ceq_grad1(1:p.numJ,1)=end_y_grad(x1).';
ceq_grad2 = zeros(size(x,1),size(x,2));
ceq_grad2(1:p.numJ,end)=end_y_grad(xm).';
ceq_grad1 = reshape(ceq_grad1,[size(x,1)*size(x,2),1]);
ceq_grad2 = reshape(ceq_grad2,[size(x,1)*size(x,2),1]);
ceq_grad=[-ceq_grad1,-ceq_grad2];
end