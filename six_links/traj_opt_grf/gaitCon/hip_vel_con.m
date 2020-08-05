function [c,gradc]=hip_vel_con(x,p)
x1 = x(:,1);
x2 = x(:,2);
x1_e = x(:,end-1);
x2_e = x(:,end);
dq1 = (x2-x1)/p.sampT;
dqe = (x2_e-x1_e)/p.sampT;
q1 = (x1+x2)/2;
qe = (x1_e+x2_e)/2;

c = hip_vel_x(q1.',dq1.',p.sampT)-hip_vel_x(qe.',dqe.',p.sampT);

if nargout>1
    gradc = zeros(size(x,1),size(x,2));
    gradc(1:p.numJ,1)=dHip_vel_q1(q1.',dq1.',p.sampT);
    gradc(1:p.numJ,2)=dHip_vel_q2(q1.',dq1.',p.sampT);
    gradc(1:p.numJ,end-1)=-dHip_vel_q1(qe.',dqe.',p.sampT);
    gradc(1:p.numJ,end)=-dHip_vel_q2(qe.',dqe.',p.sampT);
    gradc = reshape(gradc,[size(x,1)*size(x,2),1]);
end

% rewrite the hip_vel_x condition as minimum hip_vel speed

% c = zeros(size(x,2)-1,1);
% gradc=zeros(size(x,1),size(x,2),size(c,1));
% for i=1:size(x,2)-1
%     x1 = x(:,i);
%     x2 = x(:,i+1);
%     q = (x1+x2)/2;
%     dq = (x2-x1)/p.sampT;
%     c(i)=hip_vel_x(q.',dq.',p.sampT)-p.hip_vel;
%     if nargout>1
%         gradc(1:p.numJ,i,i)=dHip_vel_q1(q.',dq.',p.sampT);
%         gradc(1:p.numJ,i+1,i)=dHip_vel_q2(q.',dq.',p.sampT);
%     end
%     
%     
% end
% gradc = reshape(gradc,[size(x,1)*size(x,2),size(c,1)]);
end