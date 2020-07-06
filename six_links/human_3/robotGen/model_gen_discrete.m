% Generate robot with discrete Lagrangian


% percentage is from de Leva (1996), table 4
clear;

addpath ../

numJ = 6;
totH=1.8; %total height % still not add up to 100%, but these are just numbers
l_foot = 0.1476*totH; %use de Leva number to get the percentage 
l_calf = 0.2514*totH; %from Plagenhoef
l_thigh =0.24334*totH;
l_torso = 0.4223*totH; %head + torso, will have to adjust lc_torso later
l_heel = 0.0425*totH;  %from Rudolfs Drillis (1964) Table.1, G
DH = [0, 0, 0, 0, 0, 0;...
      0, l_calf, 0, 0, 0, 0;...
      0, l_thigh, 0, 0, 0, 0;...
      0,  0, 0, 0, 0, 0;...
      0,l_thigh, 0, 0, 0, 0;...
      0,l_calf,0,0,0,0];
robot=DH2Robot(DH,1);
robot.gravity=[0,9.81,0];


%kg
totM = 75;
m_foot = 0.0137*totM;
m_calf =0.0433*totM;
m_thigh = 0.1416*totM;
m_head = 0.0694*totM;
m_trunk = 0.4346*totM;
m_torso = m_trunk+m_head; %this total mass is not 100% though

model.totM = totM;
model.l_heel = l_heel;
model.l_foot = l_foot;

%CoM pos
lc_thigh2 = 0.433*l_thigh;
lc_calf2 = 0.433*l_calf;
lc_calf1 = l_calf-lc_calf2;
lc_thigh1 = l_thigh-lc_thigh2;

lc_foot = 0.4415*l_foot; % de Leva 
% include head when calculating, suppose neck is fixed
lc_torso = (0.4486*l_torso*m_torso+(0.1166*totH*(1-0.5976)+l_torso)*m_head)/(m_head+m_trunk); 




I_calf = [0,0,0;     %ROG data from "Comparison of direct collocation optimal control"
          0,0,0;
          0,0,m_calf*l_calf^2*0.302^2];
% I_calf1(3,3) = m_calf*l_calf^2/12;



I_thight = [0,0,0;
            0,0,0;
            0,0,m_thigh*l_thigh^2*0.323^2];
% I_thight(3,3) = m_thight*l_thight^2/12;
I_torso = [0,0,0;
           0,0,0;
           0,0,m_torso*l_torso^2*0.496^2];
% I_torso(3,3) = m_torso*l_torso^2/12;
I_foot = zeros(3);


% we need extra struct to storage link length, since for the feet, we need
% heel and foot long
end_eff = [l_foot,l_heel,l_calf];


robot.links(1).m = m_calf;
robot.links(2).m = m_thigh;
robot.links(3).m = m_torso;
robot.links(4).m = m_thigh;
robot.links(5).m = m_calf;
robot.links(6).m = m_foot;


robot.links(1).r = [lc_calf1,0,0]';
robot.links(2).r = [lc_thigh1,0,0]';
robot.links(3).r = [lc_torso,0,0]';
robot.links(4).r = [lc_thigh2,0,0]';
robot.links(5).r = [lc_calf2,0,0]';
robot.links(6).r = [lc_foot,0,0]';


robot.links(1).I = I_calf;
robot.links(2).I = I_thight;
robot.links(3).I = I_torso;
robot.links(4).I = I_thight;
robot.links(5).I = I_calf;
robot.links(6).I = I_foot;


robot.links(1).Jm=0;
robot.links(2).Jm=0;
robot.links(3).Jm=0;
robot.links(4).Jm=0;
robot.links(5).Jm=0;
robot.links(6).Jm=0;

q1 = sym('q1',[1,numJ]);
q2 = sym('q2',[1,numJ]);

syms sampT;
dyn = dynGen_discrete(robot,end_eff,sampT);


% the return values are the Discrete Lagrangian
% dL1 means dL(x_k,x_k+1)/dx_k
% dL2 means dL(x_k-1,x_k)/dx_k  (always takes derivative at x_k, 1,2 just
% menas the order of x_k)

% in other words, for dL1, it is taking derivative on q1, but dL2 is on q2
tasks = cell(1,12);
task_i = 1;
tasks{1,task_i} =@()matlabFunction(dyn.dL1(1,1),'file','dyn/dL1_1','vars',{q1,q2,sampT}); task_i = task_i+1;
tasks{1,task_i} =@()matlabFunction(dyn.dL1(2,1),'file','dyn/dL1_2','vars',{q1,q2,sampT}); task_i = task_i+1;
tasks{1,task_i} =@()matlabFunction(dyn.dL1(3,1),'file','dyn/dL1_3','vars',{q1,q2,sampT}); task_i = task_i+1;
tasks{1,task_i} =@()matlabFunction(dyn.dL1(4,1),'file','dyn/dL1_4','vars',{q1,q2,sampT}); task_i = task_i+1;
tasks{1,task_i} =@()matlabFunction(dyn.dL1(5,1),'file','dyn/dL1_5','vars',{q1,q2,sampT}); task_i = task_i+1;
tasks{1,task_i} =@()matlabFunction(dyn.dL1(6,1),'file','dyn/dL1_6','vars',{q1,q2,sampT}); task_i = task_i+1;


tasks{1,task_i} =@()matlabFunction(dyn.dL2(1,1),'file','dyn/dL2_1','vars',{q1,q2,sampT}); task_i = task_i+1;
tasks{1,task_i} =@()matlabFunction(dyn.dL2(2,1),'file','dyn/dL2_2','vars',{q1,q2,sampT}); task_i = task_i+1;
tasks{1,task_i} =@()matlabFunction(dyn.dL2(3,1),'file','dyn/dL2_3','vars',{q1,q2,sampT}); task_i = task_i+1;
tasks{1,task_i} =@()matlabFunction(dyn.dL2(4,1),'file','dyn/dL2_4','vars',{q1,q2,sampT}); task_i = task_i+1;
tasks{1,task_i} =@()matlabFunction(dyn.dL2(5,1),'file','dyn/dL2_5','vars',{q1,q2,sampT}); task_i = task_i+1;
tasks{1,task_i} =@()matlabFunction(dyn.dL2(6,1),'file','dyn/dL2_6','vars',{q1,q2,sampT}); task_i = task_i+1;

parfor i=1:length(tasks)
    tasks{1,i}();
end

% now we need to generate the derivative of it

dL1_1 = dyn.dL1(1,1);
dL1_2 = dyn.dL1(2,1);
dL1_3 = dyn.dL1(3,1);
dL1_4 = dyn.dL1(4,1);
dL1_5 = dyn.dL1(5,1);
dL1_6 = dyn.dL1(6,1);

dL2_1 = dyn.dL2(1,1);
dL2_2 = dyn.dL2(2,1);
dL2_3 = dyn.dL2(3,1);
dL2_4 = dyn.dL2(4,1);
dL2_5 = dyn.dL2(5,1);
dL2_6 = dyn.dL2(6,1);


dL11_dq1 = [diff(dL1_1,q1(1));diff(dL1_1,q1(2));diff(dL1_1,q1(3));diff(dL1_1,q1(4));diff(dL1_1,q1(5));diff(dL1_1,q1(6))];
dL12_dq1 = [diff(dL1_2,q1(1));diff(dL1_2,q1(2));diff(dL1_2,q1(3));diff(dL1_2,q1(4));diff(dL1_2,q1(5));diff(dL1_2,q1(6))];
dL13_dq1 = [diff(dL1_3,q1(1));diff(dL1_3,q1(2));diff(dL1_3,q1(3));diff(dL1_3,q1(4));diff(dL1_3,q1(5));diff(dL1_3,q1(6))];
dL14_dq1 = [diff(dL1_4,q1(1));diff(dL1_4,q1(2));diff(dL1_4,q1(3));diff(dL1_4,q1(4));diff(dL1_4,q1(5));diff(dL1_4,q1(6))];
dL15_dq1 = [diff(dL1_5,q1(1));diff(dL1_5,q1(2));diff(dL1_5,q1(3));diff(dL1_5,q1(4));diff(dL1_5,q1(5));diff(dL1_5,q1(6))];
dL16_dq1 = [diff(dL1_6,q1(1));diff(dL1_6,q1(2));diff(dL1_6,q1(3));diff(dL1_6,q1(4));diff(dL1_6,q1(5));diff(dL1_6,q1(6))];

dL11_dq2 = [diff(dL1_1,q2(1));diff(dL1_1,q2(2));diff(dL1_1,q2(3));diff(dL1_1,q2(4));diff(dL1_1,q2(5));diff(dL1_1,q2(6))];
dL12_dq2 = [diff(dL1_2,q2(1));diff(dL1_2,q2(2));diff(dL1_2,q2(3));diff(dL1_2,q2(4));diff(dL1_2,q2(5));diff(dL1_2,q2(6))];
dL13_dq2 = [diff(dL1_3,q2(1));diff(dL1_3,q2(2));diff(dL1_3,q2(3));diff(dL1_3,q2(4));diff(dL1_3,q2(5));diff(dL1_3,q2(6))];
dL14_dq2 = [diff(dL1_4,q2(1));diff(dL1_4,q2(2));diff(dL1_4,q2(3));diff(dL1_4,q2(4));diff(dL1_4,q2(5));diff(dL1_4,q2(6))];
dL15_dq2 = [diff(dL1_5,q2(1));diff(dL1_5,q2(2));diff(dL1_5,q2(3));diff(dL1_5,q2(4));diff(dL1_5,q2(5));diff(dL1_5,q2(6))];
dL16_dq2 = [diff(dL1_6,q2(1));diff(dL1_6,q2(2));diff(dL1_6,q2(3));diff(dL1_6,q2(4));diff(dL1_6,q2(5));diff(dL1_6,q2(6))];



dL21_dq1 = [diff(dL2_1,q1(1));diff(dL2_1,q1(2));diff(dL2_1,q1(3));diff(dL2_1,q1(4));diff(dL2_1,q1(5));diff(dL2_1,q1(6))];
dL22_dq1 = [diff(dL2_2,q1(1));diff(dL2_2,q1(2));diff(dL2_2,q1(3));diff(dL2_2,q1(4));diff(dL2_2,q1(5));diff(dL2_2,q1(6))];
dL23_dq1 = [diff(dL2_3,q1(1));diff(dL2_3,q1(2));diff(dL2_3,q1(3));diff(dL2_3,q1(4));diff(dL2_3,q1(5));diff(dL2_3,q1(6))];
dL24_dq1 = [diff(dL2_4,q1(1));diff(dL2_4,q1(2));diff(dL2_4,q1(3));diff(dL2_4,q1(4));diff(dL2_4,q1(5));diff(dL2_4,q1(6))];
dL25_dq1 = [diff(dL2_5,q1(1));diff(dL2_5,q1(2));diff(dL2_5,q1(3));diff(dL2_5,q1(4));diff(dL2_5,q1(5));diff(dL2_5,q1(6))];
dL26_dq1 = [diff(dL2_6,q1(1));diff(dL2_6,q1(2));diff(dL2_6,q1(3));diff(dL2_6,q1(4));diff(dL2_6,q1(5));diff(dL2_6,q1(6))];

dL21_dq2 = [diff(dL2_1,q2(1));diff(dL2_1,q2(2));diff(dL2_1,q2(3));diff(dL2_1,q2(4));diff(dL2_1,q2(5));diff(dL2_1,q2(6))];
dL22_dq2 = [diff(dL2_2,q2(1));diff(dL2_2,q2(2));diff(dL2_2,q2(3));diff(dL2_2,q2(4));diff(dL2_2,q2(5));diff(dL2_2,q2(6))];
dL23_dq2 = [diff(dL2_3,q2(1));diff(dL2_3,q2(2));diff(dL2_3,q2(3));diff(dL2_3,q2(4));diff(dL2_3,q2(5));diff(dL2_3,q2(6))];
dL24_dq2 = [diff(dL2_4,q2(1));diff(dL2_4,q2(2));diff(dL2_4,q2(3));diff(dL2_4,q2(4));diff(dL2_4,q2(5));diff(dL2_4,q2(6))];
dL25_dq2 = [diff(dL2_5,q2(1));diff(dL2_5,q2(2));diff(dL2_5,q2(3));diff(dL2_5,q2(4));diff(dL2_5,q2(5));diff(dL2_5,q2(6))];
dL26_dq2 = [diff(dL2_6,q2(1));diff(dL2_6,q2(2));diff(dL2_6,q2(3));diff(dL2_6,q2(4));diff(dL2_6,q2(5));diff(dL2_6,q2(6))];

tasks{1,task_i} =@()matlabFunction(dL11_dq1,'file','grad/dL11_dq1','vars',{q1,q2,sampT},'Optimize','false'); task_i = task_i+1;
tasks{1,task_i} =@()matlabFunction(dL12_dq1,'file','grad/dL12_dq1','vars',{q1,q2,sampT},'Optimize','false'); task_i = task_i+1;
tasks{1,task_i} =@()matlabFunction(dL13_dq1,'file','grad/dL13_dq1','vars',{q1,q2,sampT},'Optimize','false'); task_i = task_i+1;
tasks{1,task_i} =@()matlabFunction(dL14_dq1,'file','grad/dL14_dq1','vars',{q1,q2,sampT},'Optimize','false'); task_i = task_i+1;
tasks{1,task_i} =@()matlabFunction(dL15_dq1,'file','grad/dL15_dq1','vars',{q1,q2,sampT},'Optimize','false'); task_i = task_i+1;
tasks{1,task_i} =@()matlabFunction(dL16_dq1,'file','grad/dL16_dq1','vars',{q1,q2,sampT},'Optimize','false'); task_i = task_i+1;

tasks{1,task_i} =@()matlabFunction(dL11_dq2,'file','grad/dL11_dq2','vars',{q1,q2,sampT}); task_i = task_i+1;
tasks{1,task_i} =@()matlabFunction(dL12_dq2,'file','grad/dL12_dq2','vars',{q1,q2,sampT}); task_i = task_i+1;
tasks{1,task_i} =@()matlabFunction(dL13_dq2,'file','grad/dL13_dq2','vars',{q1,q2,sampT}); task_i = task_i+1;
tasks{1,task_i} =@()matlabFunction(dL14_dq2,'file','grad/dL14_dq2','vars',{q1,q2,sampT}); task_i = task_i+1;
tasks{1,task_i} =@()matlabFunction(dL15_dq2,'file','grad/dL15_dq2','vars',{q1,q2,sampT}); task_i = task_i+1;
tasks{1,task_i} =@()matlabFunction(dL16_dq2,'file','grad/dL16_dq2','vars',{q1,q2,sampT}); task_i = task_i+1;

tasks{1,task_i} =@()matlabFunction(dL21_dq1,'file','grad/dL21_dq1','vars',{q1,q2,sampT}); task_i = task_i+1;
tasks{1,task_i} =@()matlabFunction(dL22_dq1,'file','grad/dL22_dq1','vars',{q1,q2,sampT}); task_i = task_i+1;
tasks{1,task_i} =@()matlabFunction(dL23_dq1,'file','grad/dL23_dq1','vars',{q1,q2,sampT}); task_i = task_i+1;
tasks{1,task_i} =@()matlabFunction(dL24_dq1,'file','grad/dL24_dq1','vars',{q1,q2,sampT}); task_i = task_i+1;
tasks{1,task_i} =@()matlabFunction(dL25_dq1,'file','grad/dL25_dq1','vars',{q1,q2,sampT}); task_i = task_i+1;
tasks{1,task_i} =@()matlabFunction(dL26_dq1,'file','grad/dL26_dq1','vars',{q1,q2,sampT}); task_i = task_i+1;

tasks{1,task_i} =@()matlabFunction(dL21_dq2,'file','grad/dL21_dq2','vars',{q1,q2,sampT}); task_i = task_i+1;
tasks{1,task_i} =@()matlabFunction(dL22_dq2,'file','grad/dL22_dq2','vars',{q1,q2,sampT}); task_i = task_i+1;
tasks{1,task_i} =@()matlabFunction(dL23_dq2,'file','grad/dL23_dq2','vars',{q1,q2,sampT}); task_i = task_i+1;
tasks{1,task_i} =@()matlabFunction(dL24_dq2,'file','grad/dL24_dq2','vars',{q1,q2,sampT}); task_i = task_i+1;
tasks{1,task_i} =@()matlabFunction(dL25_dq2,'file','grad/dL25_dq2','vars',{q1,q2,sampT}); task_i = task_i+1;
tasks{1,task_i} =@()matlabFunction(dL26_dq2,'file','grad/dL26_dq2','vars',{q1,q2,sampT}); task_i = task_i+1;


% % conditions for ground reaction force
% % to make this symbolic calculation simplier, we use dq, although it is
% % actually (q2-q1)/sampT
% qd = sym('qd',[1,numJ]);
% q = sym('q',[1,numJ]);
% syms Fx_heel Fy_heel Fx_toe Fy_toe s_toe s_heel
% syms k cmax dmax us;
% J_toe = six_J(q(1),q(2),q(3),q(4),q(5),q(6));
% J_heel = six_J2(q(1),q(2),q(3),q(4),q(5),q(6));
% toe_vel = J_toe(1:2,:)*qd.';
% heel_vel = J_heel(1:2,:)*qd.';
% 
% toe_vel_x = toe_vel(1);
% toe_vel_y = toe_vel(2);
% heel_vel_x =heel_vel(1);
% heel_vel_y = heel_vel(2);
% 
% toe_pos = 
% 
% 
% % first constraint
% % Fy is equal to the spring-damper function when slack variables s_toe,
% % s_heel is not zero
% 
% Fn_toe_c1 = Fy_toe-s_toe*   (k*(th-toe_y_pos)^ks -(cmax*(0.5*tanh(2*(th-toe_y_pos)/dmax)+0.5))*toe_vel_y);
% Fn_heel_c1 = Fy_heel-s_heel*(k*(th-heel_y_pos)^ks-(cmax*(0.5*tanh(2*(th-heel_y_pos)/dmax)+0.5))*heel_vel_y);
% 
% tasks{1,task_i}=@()matlabFunction(Fn_toe_c1,'file','grf/discrete/Fn_toe_c1d','vars',{q,dq,[Fx_toe,Fy_toe],s_toe,th,k,cmax,dmax}); task_i=task_i+1;
% tasks{1,task_i}=@()matlabFunction(Fn_heel_c1,'file','grf/discrete/Fn_heel_c1d','vars',{q,dq,[Fx_heel,Fy_heel],s_heel,th,k,cmax,dmax}); task_i=task_i+1;
% 
% dFn_toe_c1q = [diff(Fn_toe_c1,q(1));
%                diff(Fn_toe_c1,q(2));
%                diff(Fn_toe_c1,q(3));
%                diff(Fn_toe_c1,q(4));
%                diff(Fn_toe_c1,q(5));
%                diff(Fn_toe_c1,q(6))];
% dFn_toe_c1dq=[diff(Fn_toe_c1,dq(1));
%               diff(Fn_toe_c1,dq(2));
%               diff(Fn_toe_c1,dq(3));
%               diff(Fn_toe_c1,dq(4));
%               diff(Fn_toe_c1,dq(5));
%               diff(Fn_toe_c1,dq(6))];
% dFn_toe_c1s = diff(Fn_toe_c1,s_toe);
% tasks{1,task_i}=@()matlabFunction(dFn_toe_c1q,'file','grf/discrete/dFn_toe_c1dq','vars',{q,dq,[Fx_toe,Fy_toe],s_toe,th,k,cmax,dmax}); task_i=task_i+1;
% tasks{1,task_i}=@()matlabFunction(dFn_toe_c1dq,'file','grf/discrete/Fn_toe_c1ddq','vars',{q,dq,[Fx_toe,Fy_toe],s_toe,th,k,cmax,dmax}); task_i=task_i+1;
% tasks{1,task_i}=@()matlabFunction(dFn_toe_c1s,'file','grf/discrete/Fn_toe_c1ds','vars',{q,dq,[Fx_toe,Fy_toe],s_toe,th,k,cmax,dmax}); task_i=task_i+1;
%           
%           
% dFn_heel_c1q =[diff(Fn_heel_c1,q(1));
%                diff(Fn_heel_c1,q(2));
%                diff(Fn_heel_c1,q(3));
%                diff(Fn_heel_c1,q(4));
%                diff(Fn_heel_c1,q(5));
%                diff(Fn_heel_c1,q(6))];
%            
% dFn_heel_c1dq=[diff(Fn_heel_c1,dq(1));
%                diff(Fn_heel_c1,dq(2));
%                diff(Fn_heel_c1,dq(3));
%                diff(Fn_heel_c1,dq(4));
%                diff(Fn_heel_c1,dq(5));
%                diff(Fn_heel_c1,dq(6))];

             
              
               


