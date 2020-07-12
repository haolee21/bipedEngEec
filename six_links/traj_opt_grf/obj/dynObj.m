function [dObj,dObjGrad] = dynObj(x,p,u,f_toe,f_heel)

dObjGrad = zeros(size(x,1)*size(x,2),1);
dObj =0;


[f,dfx] = f_dyn(x,p,f_toe,f_heel,u);
dObj = dObj+f;
dObjGrad = dObjGrad+dfx;
% hipLen constraint

dObj = dObj+1000*(hip_x_pos(x(:,1).')-hip_x_pos(x(:,end).')-p.hipLen)^2;

if nargout>1
    dObjGrad(1:3,1)=dObjGrad(1:3,1)+1000*hip_x_grad(x(:,1).').'*(hip_x_pos(x(:,1).')-hip_x_pos(x(:,end).')-p.hipLen)*2;
    dObjGrad(end-2*p.numJ+1:end-2*p.numJ+3,1)=dObjGrad(end-2*p.numJ+1:end-2*p.numJ+3,1)-1000*hip_x_grad(x(:,end).').'*(hip_x_pos(x(:,1).')-hip_x_pos(x(:,end).')-p.hipLen)*2;
end


% set initial toe, heel heigh
dObj = dObj ...
    + 0.5*((end_y_pos(x(:,1).')-p.init_y)^2 ...
    +(heel_y_pos(x(:,1).')-p.init_y)^2 ...
    +(end_y_pos(x(:,end).')-p.init_y)^2 ...
    +(heel_y_pos(x(:,end).')-p.init_y)^2)*p.initPos_w;
if nargout>1
    dObjGrad(1:p.numJ,1) = dObjGrad(1:p.numJ,1)...
        +end_y_grad(x(:,1).').'*end_y_pos(x(:,1).'-p.init_y)*p.initPos_w...
        +heel_y_grad(x(:,1).').'*heel_y_pos(x(:,1).'-p.init_y)*p.initPos_w;
    dObjGrad(end-p.numJ*2+1:end-p.numJ,1) = dObjGrad(end-p.numJ*2+1:end-p.numJ,1)...
        +end_y_grad(x(:,end).').'*end_y_pos(x(:,end).'-p.init_y)*p.initPos_w...
        +heel_y_grad(x(:,end).').'*heel_y_pos(x(:,end).'-p.init_y)*p.initPos_w;
    
    
end


% put init/end pos constraint here

x = reshape(x,[size(x,1)*size(x,2),1]);
err = p.Aeq*x-p.beq;
dObj = dObj + err.'*err*0.5*p.end_pos_w;
dObjGrad = dObjGrad+p.Aeq.'*err*p.end_pos_w;

% initial toe, heel y velocity =0
% J_toe_1 = six_J(x(1,1),x(2,1),x(3,1),x(4,1),x(5,1),x(6,1));
% J_heel_1 = six_J2(x(1,1),x(2,1),x(3,1),x(4,1),x(5,1),x(6,1));
% toe_vel1 = J_toe_1(2,:)*x(p.numJ+1:p.numJ*2,1);
% heel_vel1 = J_heel_1(2,:)*x(p.numJ+1:p.numJ*2,1);
%
%
% J_toe_2 = six_J(x(1,end),x(2,end),x(3,end),x(4,end),x(5,end),x(6,end));
% J_heel_2 = six_J2(x(1,end),x(2,end),x(3,end),x(4,end),x(5,end),x(6,end));
% toe_vel2 = J_toe_2(2,:)*x(p.numJ+1:p.numJ*2,end);
% heel_vel2 = J_heel_2(2,:)*x(p.numJ+1:p.numJ*2,end);
%
% dObj = dObj+(toe_vel1^2+heel_vel1^2+toe_vel2^2+heel_vel2^2)*p.init_vel_w*0.5;
%
% if nargout>1
%
%     %start
%     grad1 = dJ_q1(x(:,1).');
%     grad2 = dJ_q2(x(:,1).');
%     grad3 = dJ_q3(x(:,1).');
%     grad4 = dJ_q4(x(:,1).');
%     grad5 = dJ_q5(x(:,1).');
%     grad6 = dJ_q6(x(:,1).');
%     dJ_toe_q1 = [grad1(2,:);grad2(2,:);grad3(2,:);grad4(2,:);grad5(2,:);grad6(2,:)]*x(p.numJ+1:2*p.numJ,1);
%
%     grad1 = dJ2_q1(x(:,1).');
%     grad2 = dJ2_q2(x(:,1).');
%     grad3 = dJ2_q3(x(:,1).');
%     grad4 = dJ2_q4(x(:,1).');
%     grad5 = dJ2_q5(x(:,1).');
%     grad6 = dJ2_q6(x(:,1).');
%     dJ_heel_q1 = [grad1(2,:);grad2(2,:);grad3(2,:);grad4(2,:);grad5(2,:);grad6(2,:)]*x(p.numJ+1:2*p.numJ,1);
%
%     %end
%     grad1 = dJ_q1(x(:,end).');
%     grad2 = dJ_q2(x(:,end).');
%     grad3 = dJ_q3(x(:,end).');
%     grad4 = dJ_q4(x(:,end).');
%     grad5 = dJ_q5(x(:,end).');
%     grad6 = dJ_q6(x(:,end).');
%     dJ_toe_q2 = [grad1(2,:);grad2(2,:);grad3(2,:);grad4(2,:);grad5(2,:);grad6(2,:)]*x(p.numJ+1:2*p.numJ,end);
%
%     grad1 = dJ2_q1(x(:,end).');
%     grad2 = dJ2_q2(x(:,end).');
%     grad3 = dJ2_q3(x(:,end).');
%     grad4 = dJ2_q4(x(:,end).');
%     grad5 = dJ2_q5(x(:,end).');
%     grad6 = dJ2_q6(x(:,end).');
%     dJ_heel_q2 = [grad1(2,:);grad2(2,:);grad3(2,:);grad4(2,:);grad5(2,:);grad6(2,:)]*x(p.numJ+1:2*p.numJ,end);
%
%     %the gradient are on both position and velocity
%
%     dObjGrad(1:p.numJ*2) = [dJ_toe_q1;J_toe_1(2,:).']*p.init_vel_w*toe_vel1...
%                           +[dJ_heel_q1;J_heel_1(2,:).']*p.init_vel_w*heel_vel1...
%                           +dObjGrad(1:p.numJ*2);
%     dObjGrad(end-p.numJ*2+1:end) = [dJ_toe_q2;J_toe_2(2,:).']*p.init_vel_w*toe_vel2...
%                                   +[dJ_heel_q2;J_heel_2(2,:).']*p.init_vel_w*heel_vel2...
%                                   +dObjGrad(end-p.numJ*2+1:end);
% end


% make velocity smooth, this term should get smaller as time goes on

% dObj = dObj + p.smooth_vel*0.5*sum(x(p.numJ+1:2*p.numJ,:).^2,'all');
% if nargout>1
%     for i=1:size(x,2)
%         dObjGrad((i-1)*p.numJ*2+p.numJ+1:i*p.numJ*2,:) = p.smooth_vel* x(p.numJ+1:2*p.numJ,i)+ dObjGrad((i-1)*p.numJ*2+p.numJ+1:i*p.numJ*2,:);
%     end
% end

end