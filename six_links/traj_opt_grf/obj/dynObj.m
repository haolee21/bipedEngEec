function [dObj,dObjGrad] = dynObj(x,p,u,f_toe,f_heel)

% gnd_i1 = floor(size(x,2)*0.45);
% gnd_i2 = floor(size(x,2)*0.55);

% if nargin<3
%     q = x(1:p.numJ,:);
%     dq = x(p.numJ+1:2*p.numJ,:);
%     ddq = gradient(dq)/p.sampT;
%     
%     c.c_toe = p.k0/p.k1*ones(1,p.dataLen);
%     c.c_heel = p.k0/p.k1*ones(1,p.dataLen);
%     J.toe = cell(1,p.dataLen);
%     J.heel = cell(1,p.dataLen);
%     tau= zeros(p.numJ,p.dataLen);
%     for i =1:p.dataLen
%         cur_q = q(:,i);
%         cur_dq = dq(:,i);
%         
%         % check toe pos
%         toe_h = end_y_pos(cur_q.');
%         heel_h = heel_y_pos(cur_q.');
%         
%         J_toe = six_J(cur_q(1),cur_q(2),cur_q(3),cur_q(4),cur_q(5),cur_q(6));
%         J_heel = six_J2(cur_q(1),cur_q(2),cur_q(3),cur_q(4),cur_q(5),cur_q(6));
%         
%         J.toe{1,i} = J_toe(1:2,:);
%         J.heel{1,i} = J_heel(1:2,:);
%         
%         G = six_G(cur_q(1),cur_q(2),cur_q(3),cur_q(4),cur_q(5),cur_q(6));
%         V = six_V(cur_q(2),cur_q(3),cur_q(4),cur_q(5),cur_q(6),cur_dq(1),cur_dq(2),cur_dq(3),cur_dq(4),cur_dq(5),cur_dq(6));
%         M = six_M(cur_q(2),cur_q(3),cur_q(4),cur_q(5),cur_q(6));
%         tau(:,i) = u_no_ext([cur_q;cur_dq],ddq(:,i),p);
%         
%         toe_vel = J_toe*dq(:,i);
%         heel_vel = J_heel*dq(:,i);
%         if(toe_h<p.toe_th)
%             cur_c=p.k*(p.toe_th-toe_h)^2-p.cmax*(p.toe_th-toe_h)/p.dmax*toe_vel(2);
%             if(cur_c>0)
%                 c.c_toe(1,i) = p.k0/(cur_c^2+p.k1);
%             end
%             
%         end
%         if(heel_h<p.toe_th)
%             cur_c = p.k*(p.toe_th-heel_h)^2-p.cmax*(p.toe_th-heel_h)/p.dmax*heel_vel(2);
%             if(cur_c>0)
%                 c.c_heel(1,i)=p.k0/(cur_c^2+p.k1);
%                 
%             end
%         end
%         
%     end
%     %% solve u,f based on q,dq,ddq
%     
%     tau = reshape(tau,[size(tau,1)*size(tau,2),1]);
%     JT_toe_cell = cellfun(@transpose,J.toe,'UniformOutput',false);
%     JT_heel_cell = cellfun(@transpose,J.heel,'UniformOutput',false);
%     JT_toe_big = blkdiag(JT_toe_cell{:});
%     JT_heel_big = blkdiag(JT_heel_cell{:});
%     R = diag(repmat(p.jointW,1,p.dataLen));
%     
%     W_toe_cell = arrayfun(p.c_matFun,c.c_toe,'UniformOutput',false);
%     W_toe = blkdiag(W_toe_cell{:});
%     W_heel_cell = arrayfun(p.c_matFun,c.c_heel,'UniformOutput',false);
%     W_heel = blkdiag(W_heel_cell{:});
%     cvx_solver sedumi
%     cvx_begin quiet
%     variables f_toe(2*p.dataLen,1) f_heel(2*p.dataLen,1) u(p.numJ*p.dataLen,1)
%     minimize( 1000000*norm(JT_toe_big*f_toe+JT_heel_big*f_heel+u-tau)+f_toe.'*W_toe*f_toe+f_heel.'*W_heel*f_heel+u.'*R*u)
%     subject to
%     p.fri_mat*f_toe<=0;
%     p.fri_mat*f_heel<=0;
%     cvx_end
%     % reshape to fit the format later
%     f_toe = reshape(f_toe,[2,p.dataLen]);
%     f_heel= reshape(f_heel,[2,p.dataLen]);
%     u = reshape(u,[p.numJ,p.dataLen]);
%     
%     
%     
%     
% end


dObjGrad = zeros(size(x,1)*size(x,2),1);
dObj =0;

for i=1:size(x,2)-1
    gnd_cost =0;
%     if nargout>1
        % dynamic loss
        [f1,df1x] = f_dyn(x(:,i),p,f_toe(:,i),f_heel(:,i),u(:,i));
        [f2,df2x] = f_dyn(x(:,i+1),p,f_toe(:,i+1),f_heel(:,i+1),u(:,i+1));
        x_half = 0.5*(x(1:p.numJ*2,i)+x(1:p.numJ*2,i+1))+p.sampT/8*(f1-f2);
        u_half = 0.5*(u(:,i)+u(:,i+1));
        f_toe_half = 0.5*(f_toe(:,i)+f_toe(:,i+1));
        f_heel_half = 0.5*(f_heel(:,i)+f_heel(:,i+1));
        [f_half,df_half_x] = f_dyn(x_half,p,f_toe_half,f_heel_half,u_half);
        err = x(:,i+1)-x(:,i)-p.sampT/6*(f1+4*f_half+f2);
        
        dV1dx = -eye(p.numJ*2)-p.sampT/6*(df1x+(2*eye(p.numJ*2)+p.sampT/2*df1x)*df_half_x);
        dV2dx = eye(p.numJ*2)-p.sampT/6*(df2x+(2*eye(p.numJ*2)-p.sampT/2*df2x)*df_half_x);
        
        dObjGrad((i-1)*p.numJ*2+1:i*p.numJ*2)=p.dynCoeff*2*dV1dx*err+ dObjGrad((i-1)*p.numJ*2+1:i*p.numJ*2);
        dObjGrad(i*p.numJ*2+1:(i+1)*p.numJ*2)=p.dynCoeff*2*dV2dx*err+dObjGrad(i*p.numJ*2+1:(i+1)*p.numJ*2);
        
        
%         if(i>=gnd_i1 && i<=gnd_i2)
%             gnd_cost =  (p.gndclear - end_y_pos(x(:,i).'))*p.gndCoeff;
%             dObjGrad((i-1)*p.numJ*2+1:(i-1)*p.numJ*2+p.numJ)=dObjGrad((i-1)*p.numJ*2+1:(i-1)*p.numJ*2+p.numJ)-end_y_grad(x(:,i).').'*p.gndCoeff;
%             
%             
%         end
        
        
        
        
%     else
%         f1 = f_dyn(x(:,i),p,f_toe(:,i),f_heel(:,i),u(:,i));
%         f2 = f_dyn(x(:,i+1),p,f_toe(:,i+1),f_heel(:,i+1),u(:,i+1));
%         x_half = 0.5*(x(1:p.numJ*2,i)+x(1:p.numJ*2,i+1))+p.sampT/8*(f1-f2);
%         u_half = 0.5*(u(:,i)+u(:,i+1));
%         f_toe_half = 0.5*(f_toe(:,i)+f_toe(:,i+1));
%         f_heel_half = 0.5*(f_heel(:,i)+f_heel(:,i+1));
%         f_half = f_dyn(x_half,p,f_toe_half,f_heel_half,u_half);
%         
%         err = x(:,i+1)-x(:,i)-p.sampT/6*(f1+4*f_half+f2);
        
        
%     end
    dObj = dObj+p.dynCoeff*(err.'*err)+gnd_cost;
    
end
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