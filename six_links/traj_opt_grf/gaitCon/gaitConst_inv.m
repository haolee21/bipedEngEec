function [c,ceq,gradc,gradceq]=gaitConst_inv(x,p)
x1 = x(:,1);
xm = x(:,end);
ceq1 = hip_x_pos(x1.')-hip_x_pos(xm.')-p.hipLen; 
ceq_grad1 = zeros(size(x,1),size(x,2));
ceq_grad1(1:3,1)=hip_x_grad(x1.').';
ceq_grad1(1:3,end)=-hip_x_grad(xm.').';
ceq_grad1 = reshape(ceq_grad1(1:2*p.numJ,:),[size(x,1)*size(x,2),1]);



% init/end pos constraint

initYPos=end_y_pos(x1.');
endYPos=end_y_pos(xm.');

ceq2=[p.init_y-initYPos;p.init_y-endYPos];


ceq_grad2_1 = zeros(size(x,1),size(x,2));
ceq_grad2_1(1:p.numJ,1)=end_y_grad(x1.').';
ceq_grad2_2 = zeros(size(x,1),size(x,2));
ceq_grad2_2(1:p.numJ,end)=end_y_grad(xm.').';
ceq_grad2_1 = reshape(ceq_grad2_1,[size(x,1)*size(x,2),1]);
ceq_grad2_2 = reshape(ceq_grad2_2,[size(x,1)*size(x,2),1]);
ceq_grad2=[-ceq_grad2_1,-ceq_grad2_2];


c =[];
gradc=[];
ceq = [ceq1;ceq2];

gradceq = [ceq_grad1,ceq_grad2];



end