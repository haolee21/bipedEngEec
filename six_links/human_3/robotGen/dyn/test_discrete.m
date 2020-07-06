% this script is for verifying the accuracy of discrete lagrangian
clear;
addpath ('/home/lowlimb/cdrive/UCLA/lab/Exoskeleton/Dynamic Model/walkSim_ankle/six_links/traj_opt_grf/plotRobot')
addpath( '/home/lowlimb/cdrive/UCLA/lab/Exoskeleton/Dynamic Model/walkSim_ankle/six_links/human_3/robotGen')
% we first need to simulate a randon impulse
sampT = 0.01;
totT= 1;
ut = 0:sampT:totT;
param = load('param_contact_test.mat').param;
u = [sin(2*pi*ut+rand(1));
     sin(2*pi*ut+rand(1));
     sin(2*pi*ut+rand(1));
     sin(2*pi*ut+rand(1));
     sin(2*pi*ut+rand(1));
     sin(2*pi*ut+rand(1))]*5;
u = u+rand(size(u,1),size(u,2))*0;
fext = rand(2,size(u,2))*1;

param.knee_stiff =76.325; % I use max moment (MVC/angle), since the stiffness of the paper is too high
param.ank_stiff=408.65;

q0 = [80;0;0;0;0;-90]/180*pi;
dq0=[0;0;0;0;0;0];
opts = odeset('RelTol',1e-8,'AbsTol',1e-6);
[t,y] = ode45(@(t,x)sim_dyn(t,x,ut,u,fext,param),[0,1],[q0;dq0],opts);

y = interp1(t,y,ut,'spline').';

test = [y;u;fext;zeros(2,size(ut,2))];
drawRobot_self(test,param,figure(1));

% test discrete lagrangian
q = y(1:6,:);
dq = y(7:12,:);
err = zeros(6,length(ut)-2);
for i=1:size(ut,2)-2
    q1 = q(:,i);
    q2 = q(:,i+1);  % we centered at q2
    q3 = q(:,i+2);
    u1 = u(:,i);
    u2 = u(:,i+1);
    u3 = u(:,i+2);
    f1 = fext(:,i);
    f2 = fext(:,i+1);
    f3 = fext(:,i+2);
    tend_ank1_1 = [pi/2-q1(1,1),0,0,0,0,0]*param.ank_stiff;
    tend_ank1_2 = [pi/2-q2(1,1),0,0,0,0,0]*param.ank_stiff;
    tend_ank1_3 = [pi/2-q3(1,1),0,0,0,0,0]*param.ank_stiff;
    
    tend_ank2_1 = [0,0,0,0,0,-pi/2-q1(6,1)]*param.ank_stiff;
    tend_ank2_2 = [0,0,0,0,0,-pi/2-q2(6,1)]*param.ank_stiff;
    tend_ank2_3 = [0,0,0,0,0,-pi/2-q3(6,1)]*param.ank_stiff;
    
    
    dL1 = [dL1_1(q2.',q3.',sampT);dL1_2(q2.',q3.',sampT);dL1_3(q2.',q3.',sampT);dL1_4(q2.',q3.',sampT);dL1_5(q2.',q3.',sampT);dL1_6(q2.',q3.',sampT)]*sampT;
    dL2 = [dL2_1(q1.',q2.',sampT);dL2_2(q1.',q2.',sampT);dL2_3(q1.',q2.',sampT);dL2_4(q1.',q2.',sampT);dL2_5(q1.',q2.',sampT);dL2_6(q1.',q2.',sampT)]*sampT;
    
    u_sum = (u1+2*u2+u3)*sampT/4;
    
    J1 = six_J(q1(1),q1(2),q1(3),q1(4),q1(5),q1(6));
    J2 = six_J(q2(1),q2(2),q2(3),q2(4),q2(5),q2(6));
    J3 = six_J(q3(1),q3(2),q3(3),q3(4),q3(5),q3(6));
    
    f_sum = (J1(1:2,:).'*f1+2*J2(1:2,:).'*f2+J3(1:2,:).'*f3)*sampT/4;
    
    tend_ank1 = (tend_ank1_1+2*tend_ank1_2+tend_ank1_3).'/4*sampT;
    tend_ank2 = (tend_ank2_1+2*tend_ank2_2+tend_ank2_3).'/4*sampT;
    err(:,i) = dL1+dL2+u_sum+f_sum+tend_ank1+tend_ank2;
    
    
end


% compare it with original method
err2 = zeros(12,length(ut)-1);
for i=1:size(ut,2)-1
    q1 = q(:,i);
    dq1 = dq(:,i);
    q2 = q(:,i+1);
    dq2 = dq(:,i+1);
    u1 = u(:,i);
    u2 = u(:,i+1);
    M1 = six_M(q1(2),q1(3),q1(4),q1(5),q1(6));
    V1 = six_V(q1(2),q1(3),q1(4),q1(5),q1(6),dq1(1),dq1(2),dq1(3),dq1(4),dq1(5),dq1(6));
    G1 = six_G(q1(1),q1(2),q1(3),q1(4),q1(5),q1(6));
    
    M2 = six_M(q2(2),q2(3),q2(4),q2(5),q2(6));
    V2 = six_V(q2(2),q2(3),q2(4),q2(5),q2(6),dq2(1),dq2(2),dq2(3),dq2(4),dq2(5),dq2(6));
    G2 = six_G(q2(1),q2(2),q2(3),q2(4),q2(5),q2(6));
    ddq1 = (u1.'-V1-G1)/M1;
    ddq2 = (u2.'-V2-G2)/M2;
    
    err2(:,i) = [q2;dq2]-[q1;dq1]-([dq2;ddq2.']+[dq1;ddq1.'])*sampT/2;
    
    
    
end