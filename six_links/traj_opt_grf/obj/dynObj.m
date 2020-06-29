function [dObj,dObjGrad] = dynObj(x,p,u,f_toe,f_heel)


dObjGrad = zeros(size(x,1)*size(x,2),1);
dObj =0;
for i=1:size(x,2)-1
    
    % dynamic loss
    [f1,df1x] = f_dyn(x(:,i),p,f_toe(:,i),f_heel(:,i),u(:,i));
    [f2,df2x] = f_dyn(x(:,i+1),p,f_toe(:,i+1),f_heel(:,i+1),u(:,i));
    x_half = 0.5*(x(1:p.numJ*2,i)+x(1:p.numJ*2,i+1))+p.sampT/8*(f1-f2);
    u_half = 0.5*(u(:,i)+u(:,i+1));
    f_toe_half = 0.5*(f_toe(:,i)+f_toe(:,i+1));
    f_heel_half = 0.5*(f_heel(:,i)+f_toe(:,i+1));
    [f_half,df_half_x] = f_dyn(x_half,p,f_toe_half,f_heel_half,u_half);
    
    err = x(:,i+1)-x(:,i)-p.sampT/6*(f1+4*f_half+f2);
    dObj = dObj+err.'*err;
    
    
    dV1dx = -eye(p.numJ*2)-p.sampT/6*(df1x+(2*eye(p.numJ*2)+p.sampT/2*df1x)*df_half_x);
    dV2dx = eye(p.numJ*2)-p.sampT/6*(df2x+(2*eye(p.numJ*2)-p.sampT/2*df2x)*df_half_x);
    
    dObjGrad((i-1)*p.numJ*2+1:i*p.numJ*2)=2*dV1dx*err+ dObjGrad((i-1)*p.numJ*2+1:i*p.numJ*2);
    dObjGrad(i*p.numJ*2+1:(i+1)*p.numJ*2)=2*dV2dx*err+dObjGrad(i*p.numJ*2+1:(i+1)*p.numJ*2);
    
      
end
% hipLen constraint
% 
dObj = dObj+100000*(hip_x_pos(x(:,1).')-hip_x_pos(x(:,end).')-p.hipLen)^2;

dObjGrad(1:3,1)=dObjGrad(1:3,1)+100000*hip_x_grad(x(:,1).').'*(hip_x_pos(x(:,1).')-hip_x_pos(x(:,end).')-p.hipLen)*2;
dObjGrad(end-2*p.numJ+1:end-2*p.numJ+3,1)=dObjGrad(end-2*p.numJ+1:end-2*p.numJ+3,1)-100000*hip_x_grad(x(:,end).').'*(hip_x_pos(x(:,1).')-hip_x_pos(x(:,end).')-p.hipLen)*2;



end