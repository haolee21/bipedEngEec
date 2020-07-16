function [ceq,ceq_grad]=initYPosCons(x,p)
x1=x(:,1).';
xm=x(:,end).';

initYPos=toePos_y(x1);
endYPos=toePos_y(xm);

% ceq=[p.init_y-initYPos;p.init_y-endYPos];

ceq = endYPos-initYPos;


% ceq_grad1 = zeros(size(x,1),size(x,2));
% ceq_grad1(1:p.numJ,1)=-end_y_grad(x1).';
% ceq_grad2 = zeros(size(x,1),size(x,2));
% ceq_grad2(1:p.numJ,end)=-end_y_grad(xm).';
ceq_grad = zeros(size(x,1),size(x,2));
ceq_grad(1:p.numJ,1)=-dToePos_y(x1).';
ceq_grad(1:p.numJ,end)=dToePos_y(xm).';


% ceq_grad1 = reshape(ceq_grad1,[size(x,1)*size(x,2),1]);
% ceq_grad2 = reshape(ceq_grad2,[size(x,1)*size(x,2),1]);
% ceq_grad=[ceq_grad1,ceq_grad2];

ceq_grad = reshape(ceq_grad,[size(x,1)*size(x,2),1]);
end