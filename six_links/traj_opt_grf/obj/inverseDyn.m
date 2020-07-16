function [dObj,gradObj] = inverseDyn(u,q,dq,p)

dObj =0;
gradObj = zeros(size(u,1)*size(u,2),1);
for i=1:size(u,2)-1
    [f1,df1du] = f_forward(u(:,i),q(:,i),dq(:,i),p);
    [f2,df2du] = f_forward(u(:,i+1),q(:,i+1),dq(:,i+1),p);
    
    s_half = 0.5*([q(:,i);dq(:,i)]+[q(:,i+1);dq(:,i+1)])+p.sampT/8*(f1-f2);
    
    u_half = 0.5*(u(:,i)+u(:,i+1));
    
    [f_half,df_half_du]=f_forward(u_half,s_half(1:p.numJ),s_half(p.numJ+1:end),p);
    err = f2-f1-p.sampT/6*(f1+4*f_half+f2);
    dObj = dObj+err.'*err*0.5;
    gradObj((i-1)*p.numJ+1:i*p.numJ,1)=-1/6*p.sampT*(df1du+2*df_half_du)*err+gradObj((i-1)*p.numJ+1:i*p.numJ,1);
    gradObj(i*p.numJ+1:(i+1)*p.numJ,1)=-1/6*p.sampT*(df2du+2*df_half_du)*err+gradObj(i*p.numJ+1:(i+1)*p.numJ,1);
    
end


end