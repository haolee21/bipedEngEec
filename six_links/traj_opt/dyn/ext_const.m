function [ceq,ceq_grad] = ext_const(x,p)
%% external force at beginning need to be the same as the end

% we only need to check Fn since Fs is just related to it
x1 = x(:,1).';
x2 = x(:,end).';
[~,~,Fn_toe1,~,dFn_toe1] = toe_grf(x1,p);
[~,~,Fn_toe2,~,dFn_toe2] = toe_grf(x2,p);

[~,~,Fn_heel1,~,dFn_heel1] = heel_grf(x1,p);
[~,~,Fn_heel2,~,dFn_heel2] = heel_grf(x2,p);


ceq = [Fn_toe1-Fn_toe2,Fn_heel1-Fn_heel2].';


ceq_grad = zeros(size(x,1),size(x,2),2);

ceq_grad(:,1,1) = [dFn_toe1;zeros(p.numJ,1)];
ceq_grad(:,end,1) = [-dFn_toe2;zeros(p.numJ,1)];
ceq_grad(:,1,2) = [dFn_heel1;zeros(p.numJ,1)];
ceq_grad(:,end,2) = [-dFn_heel2;zeros(p.numJ,1)];

ceq_grad = reshape(ceq_grad,[size(x,1)*size(x,2),2]);








end