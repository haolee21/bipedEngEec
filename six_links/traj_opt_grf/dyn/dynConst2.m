function [ceq,gradceq]=dynConst2(x,p)

% both ceq and gradceq need to be reshaped later
ceq = zeros(2*p.numJ*(size(x,2)-1)-1,1);
gradceq = zeros(size(x,1),size(x,2),size(ceq,1));



for i=1:size(x,2)-1
    % since we use x as optimized variable, we use 's' for joint states 
    s1 = x(1:2*p.numJ,i);
    s2 = x(1:2*p.numJ,i+1);
    [f1,df1s,df1u,df1f] = f_x2(x(:,i),p);
    [f2,df2s,df2u,df2f] = f_x2(x(:,i+1),p);
    s_half = 0.5*(x(1:p.numJ*2,i)+x(1:p.numJ*2,i+1))+p.sampT/8*(f1-f2);
    u_half = 0.5*(x(2*p.numJ+1:3*p.numJ,i)+x(2*p.numJ+1:3*p.numJ,i+1));
    F_half = 0.5*(x(3*p.numJ+1:3*p.numJ+4,i)+x(3*p.numJ+1:3*p.numJ+4,i+1));

    x_half = [s_half;u_half;F_half];
    
    [f_half,df_half_s,df_half_u,df_half_f] = f_x2(x_half,p);
    
    ceq(2*p.numJ*(i-1)+1:2*p.numJ*i,1) = s2-s1-p.sampT*(f1+4*f_half+f2)/6;
%     ceq(2*p.numJ*(i-1)+1:2*p.numJ*i,1) = s2-s1-p.sampT*(f1+f2)/2;
    
    
    gradceq(1:2*p.numJ,i,  2*p.numJ*(i-1)+1:2*p.numJ*i) = -eye(p.numJ*2)-p.sampT/6*(df1s+(2*eye(p.numJ*2)+p.sampT/2*df1s)*df_half_s);
    gradceq(1:2*p.numJ,i+1,2*p.numJ*(i-1)+1:2*p.numJ*i)= eye(p.numJ*2)-p.sampT/6*(df2s+(2*eye(p.numJ*2)-p.sampT/2*df2s)*df_half_s);

%     gradceq(1:2*p.numJ,i,  2*p.numJ*(i-1)+1:2*p.numJ*i) = -eye(2*p.numJ)-p.sampT/2*df1s;
%     gradceq(1:2*p.numJ,i+1,2*p.numJ*(i-1)+1:2*p.numJ*i)= eye(2*p.numJ)-p.sampT/2*df2s;
%     
    gradceq(2*p.numJ+1:3*p.numJ,i,  2*p.numJ*(i-1)+1:2*p.numJ*i) = -1/6*p.sampT*(df1u+2*df_half_u);
    gradceq(2*p.numJ+1:3*p.numJ,i+1,2*p.numJ*(i-1)+1:2*p.numJ*i) = -1/6*p.sampT*(df2u+2*df_half_u);
%     gradceq(2*p.numJ+1:3*p.numJ,i,  2*p.numJ*(i-1)+1:2*p.numJ*i) = -p.sampT/2*df1u;
%     gradceq(2*p.numJ+1:3*p.numJ,i+1,2*p.numJ*(i-1)+1:2*p.numJ*i) = -p.sampT/2*df2u;
%     
    gradceq(3*p.numJ+1:3*p.numJ+4,i,  2*p.numJ*(i-1)+1:2*p.numJ*i) = -1/6*p.sampT*(df1f+2*df_half_f);
    gradceq(3*p.numJ+1:3*p.numJ+4,i+1,2*p.numJ*(i-1)+1:2*p.numJ*i) = -1/6*p.sampT*(df2f+2*df_half_f);
%     gradceq(3*p.numJ+1:3*p.numJ+4,i,  2*p.numJ*(i-1)+1:2*p.numJ*i) = -1/2*p.sampT*df1f;
%     gradceq(3*p.numJ+1:3*p.numJ+4,i+1,2*p.numJ*(i-1)+1:2*p.numJ*i) = -1/2*p.sampT*df2f;
%     
    

    
    
    
end


gradceq = reshape(gradceq,[size(x,1)*size(x,2),size(ceq,1)]);



end