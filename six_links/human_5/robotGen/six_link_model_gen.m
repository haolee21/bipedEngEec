%%% This script will generate a 2-d biped model
%   the model has 5 segments, calf, thigh, and torso
%%% 

%% This is a male model, height 180 cm, body weight 72 kg
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



syms q1 q2 q3 q4 q5 q6 th qd1 qd2 qd3 qd4 qd5 qd6 th real;
endPos = turnRTtoMatrix(robot.A([1,2,3,4,5,6],[q1 q2 q3 q4 q5 q6]))*[l_foot,l_heel,0,1].';
endPos = simplify(endPos(1:3,1));

%% end effector

% pos
end_x_pos = endPos(1,1);
end_y_pos = endPos(2,1);
matlabFunction(end_x_pos,'file','posCons/end_x_pos','Vars',[q1,q2,q3,q4,q5,q6]);
matlabFunction(end_y_pos,'file','posCons/end_y_pos','Vars',{[q1,q2,q3,q4,q5,q6]});

end_x_grad = [diff(end_x_pos,q1),diff(end_x_pos,q2),diff(end_x_pos,q3),diff(end_x_pos,q4),diff(end_x_pos,q5),diff(end_x_pos,q6)];
matlabFunction(end_x_grad,'file','posCons/end_x_grad','Vars',[q1,q2,q3,q4,q5,q6]);
end_y_grad = [diff(end_y_pos,q1),diff(end_y_pos,q2),diff(end_y_pos,q3),diff(end_y_pos,q4),diff(end_y_pos,q5),diff(end_y_pos,q6)];
matlabFunction(end_y_grad,'file','posCons/end_y_grad','Vars',{[q1,q2,q3,q4,q5,q6]});

% heel pos
heelPos = turnRTtoMatrix(robot.A([1,2,3,4,5,6],[q1 q2 q3 q4 q5 q6]))*[0,l_heel,0,1].'; 
heel_x_pos = heelPos(1,1);
heel_y_pos = heelPos(2,1);
matlabFunction(heel_x_pos,'file','posCons/heel_x_pos','vars',{[q1,q2,q3,q4,q5,q6]});
matlabFunction(heel_y_pos,'file','posCons/heel_y_pos','vars',{[q1,q2,q3,q4,q5,q6]});

heel_x_grad = [diff(heel_x_pos,q1),diff(heel_x_pos,q2),diff(heel_x_pos,q3),diff(heel_x_pos,q4),diff(heel_x_pos,q5),diff(heel_x_pos,q6)];
heel_y_grad = [diff(heel_y_pos,q1),diff(heel_y_pos,q2),diff(heel_y_pos,q3),diff(heel_y_pos,q4),diff(heel_y_pos,q5),diff(heel_y_pos,q6)];
matlabFunction(heel_x_grad,'file','posCons/heel_x_grad','vars',{[q1,q2,q3,q4,q5,q6]});
matlabFunction(heel_y_grad,'file','posCons/heel_y_grad','vars',{[q1,q2,q3,q4,q5,q6]});

% ankpos
% ankPos = turnRTtoMatrix(robot.A([1,2,3,4,5,6],[q1 q2 q3 q4 q5 q6]))*[0,0,0,1].';
% ankPos = simplify(ankPos(1:3,1));
% ank_y_pos = ankPos(2,1);
% matlabFunction(ank_y_pos,'file','posCons/ank_y_pos','Vars',{[q1,q2,q3,q4,q5,q6]});
% ank_y_grad = [diff(ank_y_pos,q1),diff(ank_y_pos,q2),diff(ank_y_pos,q3),diff(ank_y_pos,q4),diff(ank_y_pos,q5),diff(ank_y_pos,q6)];
% matlabFunction(ank_y_grad,'file','posCons/ank_y_grad','Vars',{[q1,q2,q3,q4,q5,q6]});
%% head

headPos = turnRTtoMatrix(robot.A([1,2,3],[q1,q2,q3]))*[l_torso,0,0,1].';
headPos = simplify(headPos(1:3,1));
head_y_pos = headPos(2);
matlabFunction(head_y_pos,'File','posCons/head_y_pos','vars',{[q1,q2,q3]});
% generate the gradient, but only to q1 q2 q3 (form the real gradient in
% the function since we may add more joints to the end
head_y_grad = [diff(head_y_pos,q1),diff(head_y_pos,q2),diff(head_y_pos,q3)];
matlabFunction(head_y_grad,'file','posCons/head_y_grad','vars',{[q1,q2,q3]});

%% hip
hipPos = turnRTtoMatrix(robot.A([1,2,3],[q1,q2,q3]))*[0,0,0,1].';
hipPos = simplify(hipPos(1:3,1));
hip_x_pos = hipPos(1);
matlabFunction(hip_x_pos,'File','posCons/hip_x_pos','vars',{[q1,q2,q3]});
hip_x_grad = [diff(hip_x_pos,q1),diff(hip_x_pos,q2),diff(hip_x_pos,q3)];
matlabFunction(hip_x_grad,'file','posCons/hip_x_grad','vars',{[q1,q2,q3]});


% %% activation function (sigma)
% sigma_toe = (0.5*tanh(400*(th-end_y_pos))+0.5);
% matlabFunction(sigma_toe,'file','sigma_toe','vars',{[q1,q2,q3,q4,q5,q6],th});
% dsigma_toe = [diff(sigma_toe,q1);diff(sigma_toe,q2);diff(sigma_toe,q3);diff(sigma_toe,q4);diff(sigma_toe,q5);diff(sigma_toe,q6)];
% matlabFunction(dsigma_toe,'file','dsigma_toe','vars',{[q1,q2,q3,q4,q5,q6],th});
% 
% sigma_heel =  (0.5*tanh(400*(th-heel_y_pos))+0.5);
% matlabFunction(sigma_heel,'file','sigma_heel','vars',{[q1,q2,q3,q4,q5,q6],th});
% dsigma_heel = [diff(sigma_heel,q1);diff(sigma_heel,q2);diff(sigma_heel,q3);diff(sigma_heel,q4);diff(sigma_heel,q5);diff(sigma_heel,q6)];
% matlabFunction(dsigma_heel,'file','dsigma_heel','vars',{[q1,q2,q3,q4,q5,q6],th});


%% generate the joint/CoM pos (base frame) for graphing
% Joint
knee_front = turnRTtoMatrix(robot.A(1,q1))*[l_calf,0,0,1].';
hip_front = turnRTtoMatrix(robot.A([1,2],[q1,q2]))*[l_thigh,0,0,1].';
head = turnRTtoMatrix(robot.A([1,2,3],[q1,q2,q3]))*[l_torso,0,0,1].';
knee_back = turnRTtoMatrix(robot.A([1,2,3,4],[q1,q2,q3,q4]))*[l_thigh,0,0,1].';
ankle_back =turnRTtoMatrix(robot.A([1,2,3,4,5],[q1,q2,q3,q4,q5]))*[l_calf,0,0,1].';
toe_back = turnRTtoMatrix(robot.A([1,2,3,4,5,6],[q1,q2,q3,q4,q5,q6]))*[l_foot,l_heel,0,1].';
heel_back = turnRTtoMatrix(robot.A([1,2,3,4,5,6],[q1,q2,q3,q4,q5,q6]))*[0,l_heel,0,1].';
% Com
calf_front = turnRTtoMatrix(robot.A(1,q1))*[lc_calf1,0,0,1].';
thigh_front =turnRTtoMatrix(robot.A([1,2],[q1,q2]))*[lc_thigh1,0,0,1].';
torso = turnRTtoMatrix(robot.A([1,2,3],[q1,q2,q3]))*[lc_torso,0,0,1].';
thigh_back = turnRTtoMatrix(robot.A([1,2,3,4],[q1,q2,q3,q4]))*[lc_thigh2,0,0,1].';
calf_back = turnRTtoMatrix(robot.A([1,2,3,4,5],[q1,q2,q3,q4,q5]))*[lc_calf2,0,0,1].';
foot_back = turnRTtoMatrix(robot.A([1,2,3,4,5,6],[q1,q2,q3,q4,q5,q6]))*[lc_foot,l_heel,0,1].';

draw_pos = [knee_front,hip_front,head,knee_back,ankle_back,toe_back,heel_back;...
            calf_front,thigh_front,torso,thigh_back,calf_back,foot_back,zeros(4,1)]; %2 points to draw for the feet (toe and heel)
draw_pos = [draw_pos(1:3,:);draw_pos(5:8,:)];
matlabFunction(draw_pos,'file','getRobotPos','vars',[q1,q2,q3,q4,q5,q6]);
%q6 is never important since 
% 1. we directly calculate the ankle torque with jacobian
% 2. we have no feet link, the mass and inertia are 0 for that joint

dyn = dynGen(robot,end_eff);
% 
% 
task_i = 1;
tasks{1,task_i}=@()matlabFunction(dyn.M,'File','dyn/six_M','vars',[q2,q3,q4,q5,q6]);task_i=task_i+1;
tasks{1,task_i}=@()matlabFunction(dyn.G,'File','dyn/six_G','vars',[q1 q2 q3 q4 q5 q6]);task_i=task_i+1;
tasks{1,task_i}=@()matlabFunction(dyn.J,'File','dyn/six_J','vars',[q1 q2 q3 q4 q5 q6]);task_i=task_i+1;
tasks{1,task_i}=@()matlabFunction(dyn.V,'File','dyn/six_V','vars',[q2 q3 q4 q5 q6 qd1 qd2 qd3 qd4 qd5 qd6]);task_i=task_i+1;
tasks{1,task_i}=@()matlabFunction(dyn.J2,'file','dyn/six_J2','vars',[q1,q2,q3,q4,q5,q6]);task_i=task_i+1;
tasks{1,task_i}=@()matlabFunction(dyn.J_hip,'file','dyn/six_J_hip','vars',[q1,q2,q3,q4,q5,q6]);task_i=task_i+1;
% 
% 
% %% generate gradient for the M,G,J,V
G=dyn.G;
M=dyn.M;
V=dyn.V;
J=dyn.J;
J2 = dyn.J2; %jacobian on heel
J_hip = dyn.J_hip;

dG_dx = [diff(G,q1);diff(G,q2);diff(G,q3);diff(G,q4);diff(G,q5);diff(G,q6)];

tasks{1,task_i}=@()matlabFunction(dG_dx,'file','grad/dG_dx','vars',[q1,q2,q3,q4,q5,q6]);task_i=task_i+1;       



dV_dx = [diff(V,q1);diff(V,q2);diff(V,q3);diff(V,q4);diff(V,q5);diff(V,q6);diff(V,qd1);diff(V,qd2);diff(V,qd3);diff(V,qd4);diff(V,qd5);diff(V,qd6)];

tasks{1,task_i}=@()matlabFunction(dV_dx,'file','grad/dV_dx','vars',[q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6]);task_i=task_i+1;


dM_dx2 = diff(M,q2);
dM_dx3 = diff(M,q3);
dM_dx4 = diff(M,q4);
dM_dx5 = diff(M,q5);
dM_dx6 = diff(M,q6);
tasks{1,task_i}=@()matlabFunction(dM_dx2,'file','grad/dM_dx2','vars',[q2,q3,q4,q5,q6]);task_i=task_i+1;
tasks{1,task_i}=@()matlabFunction(dM_dx3,'file','grad/dM_dx3','vars',[q2,q3,q4,q5,q6]);task_i=task_i+1;
tasks{1,task_i}=@()matlabFunction(dM_dx4,'file','grad/dM_dx4','vars',[q2,q3,q4,q5,q6]);task_i=task_i+1;
tasks{1,task_i}=@()matlabFunction(dM_dx5,'file','grad/dM_dx5','vars',[q2,q3,q4,q5,q6]);task_i=task_i+1;
tasks{1,task_i}=@()matlabFunction(dM_dx6,'file','grad/dM_dx6','vars',[q2,q3,q4,q5,q6]);task_i=task_i+1;      

% parfor i=1:length(tasks)
%     tasks{1,i}();
% end
%% generate the beta function for GRF

% f_ext on toe
% J = dyn.J;
% J2 =dyn.J2;

% J_f_toe = J(1:2,:).';
% 
% tasks = cell(1,32);
% 
% 
% 
% Mext_toe = (J_f_toe.'*J_f_toe)\J_f_toe.';
% beta_grf_toe = sigma_toe*J_f_toe;
% 
% dMext_toe_dx1 = diff(Mext_toe,q1);
% dMext_toe_dx2 = diff(Mext_toe,q2);
% dMext_toe_dx3 = diff(Mext_toe,q3);
% dMext_toe_dx4 = diff(Mext_toe,q4);
% dMext_toe_dx5 = diff(Mext_toe,q5);
% dMext_toe_dx6 = diff(Mext_toe,q6);
% 
% 
% tasks{1,1}= @()matlabFunction(Mext_toe,'file','grf/Mext_toe','vars',{[q1,q2,q3,q4,q5,q6]});
% tasks{1,2}= @()matlabFunction(dMext_toe_dx1,'file','grf/dMext_toe_dx1','vars',{[q1,q2,q3,q4,q5,q6]});
% tasks{1,3}= @()matlabFunction(dMext_toe_dx2,'file','grf/dMext_toe_dx2','vars',{[q1,q2,q3,q4,q5,q6]});
% tasks{1,4}= @()matlabFunction(dMext_toe_dx3,'file','grf/dMext_toe_dx3','vars',{[q1,q2,q3,q4,q5,q6]});
% tasks{1,5}= @()matlabFunction(dMext_toe_dx4,'file','grf/dMext_toe_dx4','vars',{[q1,q2,q3,q4,q5,q6]});
% tasks{1,6}= @()matlabFunction(dMext_toe_dx5,'file','grf/dMext_toe_dx5','vars',{[q1,q2,q3,q4,q5,q6]});
% tasks{1,7}= @()matlabFunction(dMext_toe_dx6,'file','grf/dMext_toe_dx6','vars',{[q1,q2,q3,q4,q5,q6]});
% 
% 
% 
% 
% % beta_grf = sigma*J_f/(J_f.'*J_f)*J_f.';
% dbeta_toe_dx1 = diff(beta_grf_toe,q1);
% dbeta_toe_dx2 = diff(beta_grf_toe,q2);
% dbeta_toe_dx3 = diff(beta_grf_toe,q3);
% dbeta_toe_dx4 = diff(beta_grf_toe,q4);
% dbeta_toe_dx5 = diff(beta_grf_toe,q5);
% dbeta_toe_dx6 = diff(beta_grf_toe,q6);
% 
% tasks{1,8}=@()matlabFunction(beta_grf_toe,'file','grf/beta_grf_toe','vars',{[q1,q2,q3,q4,q5,q6],th});
% tasks{1,9}=@()matlabFunction(dbeta_toe_dx1,'file','grf/dbeta_toe_dx1','vars',{[q1,q2,q3,q4,q5,q6],th});
% tasks{1,10}=@()matlabFunction(dbeta_toe_dx2,'file','grf/dbeta_toe_dx2','vars',{[q1,q2,q3,q4,q5,q6],th});
% tasks{1,11}=@()matlabFunction(dbeta_toe_dx3,'file','grf/dbeta_toe_dx3','vars',{[q1,q2,q3,q4,q5,q6],th});
% tasks{1,12}=@()matlabFunction(dbeta_toe_dx4,'file','grf/dbeta_toe_dx4','vars',{[q1,q2,q3,q4,q5,q6],th});
% tasks{1,13}=@()matlabFunction(dbeta_toe_dx5,'file','grf/dbeta_toe_dx5','vars',{[q1,q2,q3,q4,q5,q6],th});
% tasks{1,14}=@()matlabFunction(dbeta_toe_dx6,'file','grf/dbeta_toe_dx6','vars',{[q1,q2,q3,q4,q5,q6],th});
% 
% 
% % external force on heel
% J_f_heel = J2(1:2,:).';
% Mext_heel = (J_f_heel.'*J_f_heel)\J_f_heel.';
% beta_grf_heel = sigma_heel*J_f_heel;
% 
% dMext_heel_dx1 = diff(Mext_heel,q1);
% dMext_heel_dx2 = diff(Mext_heel,q2);
% dMext_heel_dx3 = diff(Mext_heel,q3);
% dMext_heel_dx4 = diff(Mext_heel,q4);
% dMext_heel_dx5 = diff(Mext_heel,q5);
% dMext_heel_dx6 = diff(Mext_heel,q6);
% 
% tasks{1,15}=@()matlabFunction(Mext_heel,'file','grf/Mext_heel','vars',{[q1,q2,q3,q4,q5,q6]});
% tasks{1,16}=@()matlabFunction(dMext_heel_dx1,'file','grf/dMext_heel_dx1','vars',{[q1,q2,q3,q4,q5,q6]});
% tasks{1,17}=@()matlabFunction(dMext_heel_dx2,'file','grf/dMext_heel_dx2','vars',{[q1,q2,q3,q4,q5,q6]});
% tasks{1,18}=@()matlabFunction(dMext_heel_dx3,'file','grf/dMext_heel_dx3','vars',{[q1,q2,q3,q4,q5,q6]});
% tasks{1,19}=@()matlabFunction(dMext_heel_dx4,'file','grf/dMext_heel_dx4','vars',{[q1,q2,q3,q4,q5,q6]});
% tasks{1,20}=@()matlabFunction(dMext_heel_dx5,'file','grf/dMext_heel_dx5','vars',{[q1,q2,q3,q4,q5,q6]});
% tasks{1,21}=@()matlabFunction(dMext_heel_dx6,'file','grf/dMext_heel_dx6','vars',{[q1,q2,q3,q4,q5,q6]});
% % beta_grf = sigma*J_f/(J_f.'*J_f)*J_f.';
% dbeta_heel_dx1 = diff(beta_grf_heel,q1);
% dbeta_heel_dx2 = diff(beta_grf_heel,q2);
% dbeta_heel_dx3 = diff(beta_grf_heel,q3);
% dbeta_heel_dx4 = diff(beta_grf_heel,q4);
% dbeta_heel_dx5 = diff(beta_grf_heel,q5);
% dbeta_heel_dx6 = diff(beta_grf_heel,q6);
% 
% tasks{1,22}=@()matlabFunction(beta_grf_heel,'file','grf/beta_grf_heel','vars',{[q1,q2,q3,q4,q5,q6],th});
% tasks{1,23}=@()matlabFunction(dbeta_heel_dx1,'file','grf/dbeta_heel_dx1','vars',{[q1,q2,q3,q4,q5,q6],th});
% tasks{1,24}=@()matlabFunction(dbeta_heel_dx2,'file','grf/dbeta_heel_dx2','vars',{[q1,q2,q3,q4,q5,q6],th});
% tasks{1,25}=@()matlabFunction(dbeta_heel_dx3,'file','grf/dbeta_heel_dx3','vars',{[q1,q2,q3,q4,q5,q6],th});
% tasks{1,26}=@()matlabFunction(dbeta_heel_dx4,'file','grf/dbeta_heel_dx4','vars',{[q1,q2,q3,q4,q5,q6],th});
% tasks{1,27}=@()matlabFunction(dbeta_heel_dx5,'file','grf/dbeta_heel_dx5','vars',{[q1,q2,q3,q4,q5,q6],th});
% tasks{1,28}=@()matlabFunction(dbeta_heel_dx6,'file','grf/dbeta_heel_dx6','vars',{[q1,q2,q3,q4,q5,q6],th});
% 
% %% generate friction for object function 
% % we only consider the horizontal friction, so we use x_dot = J(1,:)*q_dot
% 
% %toe
% x_vel1 = J(1,:)*[qd1;qd2;qd3;qd4;qd5;qd6];
% fri1 = sigma_toe*(x_vel1*x_vel1);
% fri1_dx = [diff(fri1,q1);
%           diff(fri1,q2);
%           diff(fri1,q3);
%           diff(fri1,q4);
%           diff(fri1,q5);
%           diff(fri1,q6);
%           diff(fri1,qd1);
%           diff(fri1,qd2);
%           diff(fri1,qd3);
%           diff(fri1,qd4);
%           diff(fri1,qd5);
%           diff(fri1,qd6);
%           zeros(numJ,1)];
% tasks{1,29}=@()matlabFunction(fri1,'file','obj/fri_toe','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],th});
% tasks{1,30}=@()matlabFunction(fri1_dx,'file','obj/fri_toe_dx','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],th});
% 
% % ankle
% 
% x_vel2 = J2(1,:)*[qd1;qd2;qd3;qd4;qd5;qd6];
% fri2 = sigma_heel*(x_vel2*x_vel2);
% fri2_dx = [diff(fri2,q1);
%           diff(fri2,q2);
%           diff(fri2,q3);
%           diff(fri2,q4);
%           diff(fri2,q5);
%           diff(fri2,q6);
%           diff(fri2,qd1);
%           diff(fri2,qd2);
%           diff(fri2,qd3);
%           diff(fri2,qd4);
%           diff(fri2,qd5);
%           diff(fri2,qd6);
%           zeros(numJ,1)];
% tasks{1,31}=@()matlabFunction(fri2,'file','obj/fri_heel','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],th});
% tasks{1,32}=@()matlabFunction(fri2_dx,'file','obj/fri_heel_dx','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],th});
% 
% 
% %% front/back knee locking
% % we will add a torsional spring at the front knee
% % it is triggered when q2<1 deg
% 
% syms kne_stop ank_stop qkne qank real ;
% sigma_knee = 0.5*tanh(250*(kne_stop-qkne))+0.5;
% dsigma_knee = diff(sigma_knee,qkne);
% sigma_ank = 0.5*tanh(250*(ank_stop-qank))+0.5;
% dsigma_ank = diff(sigma_ank,qank);
% 
% task{1,33}=@()matlabFunction(sigma_knee,'file','knee_spring/sigma_knee','vars',[qkne,kne_stop]);
% task{1,34}=@()matlabFunction(dsigma_knee,'file','knee_spring/dsigma_knee','vars',[qkne,kne_stop]);
% task{1,35}=@()matlabFunction(sigma_ank,'file','knee_spring/sigma_ank','vars',[qank,ank_stop]);
% task{1,36}=@()matlabFunction(dsigma_ank,'file','knee_spring/dsigma_ank','vars',[qank,ank_stop]);




%% GRF 
% since for grf we need x_vel, we generate it at last

syms k cmax dmax us ud;
% k = 2e5;
% k = 2e6;
ks = 2; %use ks to replace the spring constant e in the paper
% cmax = 1250;
% 
% dmax = 1e-2;

toe_y_pos = end_y_pos;


% here we introduce the nondifferentiable term into the model
% while it is not differentiable, it is at least continuous (relu like
% function)

% toe

sigma_Fn_toe = 0.5*tanh(400*(th-toe_y_pos))+0.5;  % although we want a discrete behavior here, but make it change too drastic will cause instability 

toe_vel = J(1:2,:)*[qd1;qd2;qd3;qd4;qd5;qd6];
y_vel_toe = toe_vel(2,1);
x_vel_toe = toe_vel(1,1);

sigma_cmax = cmax*(0.5*tanh(2*(th-toe_y_pos)/dmax)+0.5);
Fn_toe = sigma_Fn_toe*(k*(th-toe_y_pos)^ks-sigma_cmax*y_vel_toe);




dFn_toe = [diff(Fn_toe,q1);diff(Fn_toe,q2);diff(Fn_toe,q3);diff(Fn_toe,q4);diff(Fn_toe,q5);diff(Fn_toe,q6);...
            diff(Fn_toe,qd1);diff(Fn_toe,qd2);diff(Fn_toe,qd3);diff(Fn_toe,qd4);diff(Fn_toe,qd5);diff(Fn_toe,qd6)];

% I give up on having static/dynamic friction constant since this add up
% too much complicity, and us=0.8 while ud = 0.7
        
sigma_u = tanh(500*x_vel_toe)*(us-(us-ud)*((0.5*tanh(100*(x_vel_toe-0.05))+0.5)-(0.5*tanh(100*(x_vel_toe+0.05))-0.5))); %by xvel the direction will change

Fs_toe = sigma_u*Fn_toe;
        
dFs_toe = [diff(Fs_toe,q1);diff(Fs_toe,q2);diff(Fs_toe,q3);diff(Fs_toe,q4);diff(Fs_toe,q5);diff(Fs_toe,q6);...
              diff(Fs_toe,qd1);diff(Fs_toe,qd2);diff(Fs_toe,qd3);diff(Fs_toe,qd4);diff(Fs_toe,qd5);diff(Fs_toe,qd6)];   

          

tasks{1,task_i} = @()matlabFunction(Fn_toe,'file','grf/Fn_toe','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],th,k,cmax,dmax,us,ud}); task_i=task_i+1;
tasks{1,task_i} = @()matlabFunction(dFn_toe,'file','grf/dFn_toe','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],th,k,cmax,dmax,us,ud});  task_i=task_i+1;
tasks{1,task_i} = @()matlabFunction(Fs_toe,'file','grf/Fs_toe','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],th,k,cmax,dmax,us,ud}); task_i=task_i+1;
tasks{1,task_i} = @()matlabFunction(dFs_toe,'file','grf/dFs_toe','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],th,k,cmax,dmax,us,ud});  task_i=task_i+1;


% heel
sigma_Fn_heel = 0.5*tanh(400*(th-heel_y_pos))+0.5;

heel_vel = J2(1:2,:)*[qd1;qd2;qd3;qd4;qd5;qd6];
y_vel_heel = heel_vel(2,1);
x_vel_heel = heel_vel(1,1);


sigma_cmax = cmax*0.5*tanh(2*(th-heel_y_pos)/dmax)+0.5;
Fn_heel = sigma_Fn_heel*(k*(th-heel_y_pos)^ks-sigma_cmax*y_vel_heel);

dFn_heel1 = diff(Fn_heel,q1);
dFn_heel2 = diff(Fn_heel,q2);
dFn_heel3 = diff(Fn_heel,q3);
dFn_heel4 = diff(Fn_heel,q4);
dFn_heel5 = diff(Fn_heel,q5);
dFn_heel6 = diff(Fn_heel,q6);
dFn_heel7 = diff(Fn_heel,qd1);
dFn_heel8 = diff(Fn_heel,qd2);
dFn_heel9 = diff(Fn_heel,qd3);
dFn_heel10 = diff(Fn_heel,qd4);
dFn_heel11 = diff(Fn_heel,qd5);
dFn_heel12 = diff(Fn_heel,qd6);


dFn_heel = [diff(Fn_heel,q1);diff(Fn_heel,q2);diff(Fn_heel,q3);diff(Fn_heel,q4);diff(Fn_heel,q5);diff(Fn_heel,q6);...
            diff(Fn_heel,qd1);diff(Fn_heel,qd2);diff(Fn_heel,qd3);diff(Fn_heel,qd4);diff(Fn_heel,qd5);diff(Fn_heel,qd6)];

sigma_u = tanh(500*x_vel_heel)*(us-(us-ud)*((0.5*tanh(100*(x_vel_heel-0.05))+0.5)-(0.5*tanh(100*(x_vel_heel+0.05))-0.5))); %by xvel the direction will change


Fs_heel = sigma_u*Fn_heel;


dFs_heel1 = diff(Fs_heel,q1);
dFs_heel2 = diff(Fs_heel,q2);
dFs_heel3 = diff(Fs_heel,q3);
dFs_heel4 = diff(Fs_heel,q4);
dFs_heel5 = diff(Fs_heel,q5);
dFs_heel6 = diff(Fs_heel,q6);
dFs_heel7 = diff(Fs_heel,qd1);
dFs_heel8 = diff(Fs_heel,qd2);
dFs_heel9 = diff(Fs_heel,qd3);
dFs_heel10 = diff(Fs_heel,qd4);
dFs_heel11 = diff(Fs_heel,qd5);
dFs_heel12 = diff(Fs_heel,qd6);
% dFs_heel = [diff(Fs_heel,q1);diff(Fs_heel,q2);diff(Fs_heel,q3);diff(Fs_heel,q4);diff(Fs_heel,q5);diff(Fs_heel,q6);...
%               diff(Fs_heel,qd1);diff(Fs_heel,qd2);diff(Fs_heel,qd3);diff(Fs_heel,qd4);diff(Fs_heel,qd5);diff(Fs_heel,qd6)];

          
% I have to make dFn_heel, dFs_heel into seperate files to expedite the code, for some reason it is taking around 3 hours 
tasks{1,task_i} = @()matlabFunction(Fn_heel,'file','grf/Fn_heel','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],th,k,cmax,dmax,us,ud}); task_i=task_i+1;
tasks{1,task_i} = @()matlabFunction(dFn_heel1,'file','grf/dFn_heel1','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],th,k,cmax,dmax,us,ud}); task_i=task_i+1;
tasks{1,task_i} = @()matlabFunction(dFn_heel2,'file','grf/dFn_heel2','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],th,k,cmax,dmax,us,ud}); task_i=task_i+1;
tasks{1,task_i} = @()matlabFunction(dFn_heel3,'file','grf/dFn_heel3','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],th,k,cmax,dmax,us,ud}); task_i=task_i+1;
tasks{1,task_i} = @()matlabFunction(dFn_heel4,'file','grf/dFn_heel4','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],th,k,cmax,dmax,us,ud}); task_i=task_i+1;
tasks{1,task_i} = @()matlabFunction(dFn_heel5,'file','grf/dFn_heel5','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],th,k,cmax,dmax,us,ud}); task_i=task_i+1;
tasks{1,task_i} = @()matlabFunction(dFn_heel6,'file','grf/dFn_heel6','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],th,k,cmax,dmax,us,ud}); task_i=task_i+1;
tasks{1,task_i} = @()matlabFunction(dFn_heel7,'file','grf/dFn_heel7','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],th,k,cmax,dmax,us,ud}); task_i=task_i+1;
tasks{1,task_i} = @()matlabFunction(dFn_heel8,'file','grf/dFn_heel8','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],th,k,cmax,dmax,us,ud}); task_i=task_i+1;
tasks{1,task_i} = @()matlabFunction(dFn_heel9,'file','grf/dFn_heel9','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],th,k,cmax,dmax,us,ud}); task_i=task_i+1;
tasks{1,task_i} = @()matlabFunction(dFn_heel10,'file','grf/dFn_heel10','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],th,k,cmax,dmax,us,ud}); task_i=task_i+1;
tasks{1,task_i} = @()matlabFunction(dFn_heel11,'file','grf/dFn_heel11','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],th,k,cmax,dmax,us,ud}); task_i=task_i+1;
tasks{1,task_i} = @()matlabFunction(dFn_heel12,'file','grf/dFn_heel12','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],th,k,cmax,dmax,us,ud}); task_i=task_i+1;
tasks{1,task_i} = @()matlabFunction(Fs_heel,'file','grf/Fs_heel','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],th,k,cmax,dmax,us,ud}); task_i=task_i+1;

tasks{1,task_i} = @()matlabFunction(dFs_heel1,'file','grf/dFs_heel1','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],th,k,cmax,dmax,us,ud}); task_i=task_i+1;
tasks{1,task_i} = @()matlabFunction(dFs_heel2,'file','grf/dFs_heel2','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],th,k,cmax,dmax,us,ud}); task_i=task_i+1;
tasks{1,task_i} = @()matlabFunction(dFs_heel3,'file','grf/dFs_heel3','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],th,k,cmax,dmax,us,ud}); task_i=task_i+1;
tasks{1,task_i} = @()matlabFunction(dFs_heel4,'file','grf/dFs_heel4','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],th,k,cmax,dmax,us,ud}); task_i=task_i+1;
tasks{1,task_i} = @()matlabFunction(dFs_heel5,'file','grf/dFs_heel5','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],th,k,cmax,dmax,us,ud}); task_i=task_i+1;
tasks{1,task_i} = @()matlabFunction(dFs_heel6,'file','grf/dFs_heel6','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],th,k,cmax,dmax,us,ud}); task_i=task_i+1;
tasks{1,task_i} = @()matlabFunction(dFs_heel7,'file','grf/dFs_heel7','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],th,k,cmax,dmax,us,ud}); task_i=task_i+1;
tasks{1,task_i} = @()matlabFunction(dFs_heel8,'file','grf/dFs_heel8','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],th,k,cmax,dmax,us,ud}); task_i=task_i+1;
tasks{1,task_i} = @()matlabFunction(dFs_heel9,'file','grf/dFs_heel9','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],th,k,cmax,dmax,us,ud}); task_i=task_i+1;
tasks{1,task_i} = @()matlabFunction(dFs_heel10,'file','grf/dFs_heel10','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],th,k,cmax,dmax,us,ud}); task_i=task_i+1;
tasks{1,task_i} = @()matlabFunction(dFs_heel11,'file','grf/dFs_heel11','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],th,k,cmax,dmax,us,ud}); task_i=task_i+1;
tasks{1,task_i} = @()matlabFunction(dFs_heel12,'file','grf/dFs_heel12','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],th,k,cmax,dmax,us,ud}); task_i=task_i+1;
 
% we need to take derivative of the jacobian
dJ_q1 = diff(J,q1);
dJ_q2 = diff(J,q2);
dJ_q3 = diff(J,q3);
dJ_q4 = diff(J,q4);
dJ_q5 = diff(J,q5);
dJ_q6 = diff(J,q6);

tasks{1,task_i}=@()matlabFunction(dJ_q1,'file','grad/dJ_q1','vars',{[q1,q2,q3,q4,q5,q6]}); task_i=task_i+1;
tasks{1,task_i}=@()matlabFunction(dJ_q2,'file','grad/dJ_q2','vars',{[q1,q2,q3,q4,q5,q6]}); task_i=task_i+1;
tasks{1,task_i}=@()matlabFunction(dJ_q3,'file','grad/dJ_q3','vars',{[q1,q2,q3,q4,q5,q6]}); task_i=task_i+1;
tasks{1,task_i}=@()matlabFunction(dJ_q4,'file','grad/dJ_q4','vars',{[q1,q2,q3,q4,q5,q6]}); task_i=task_i+1;
tasks{1,task_i}=@()matlabFunction(dJ_q5,'file','grad/dJ_q5','vars',{[q1,q2,q3,q4,q5,q6]}); task_i=task_i+1;
tasks{1,task_i}=@()matlabFunction(dJ_q6,'file','grad/dJ_q6','vars',{[q1,q2,q3,q4,q5,q6]}); task_i=task_i+1;

dJ2_q1 = diff(J2,q1);
dJ2_q2 = diff(J2,q2);
dJ2_q3 = diff(J2,q3);
dJ2_q4 = diff(J2,q4);
dJ2_q5 = diff(J2,q5);
dJ2_q6 = diff(J2,q6);      

tasks{1,task_i}=@()matlabFunction(dJ2_q1,'file','grad/dJ2_q1','vars',{[q1,q2,q3,q4,q5,q6]}); task_i=task_i+1;
tasks{1,task_i}=@()matlabFunction(dJ2_q2,'file','grad/dJ2_q2','vars',{[q1,q2,q3,q4,q5,q6]}); task_i=task_i+1;
tasks{1,task_i}=@()matlabFunction(dJ2_q3,'file','grad/dJ2_q3','vars',{[q1,q2,q3,q4,q5,q6]}); task_i=task_i+1;
tasks{1,task_i}=@()matlabFunction(dJ2_q4,'file','grad/dJ2_q4','vars',{[q1,q2,q3,q4,q5,q6]}); task_i=task_i+1;
tasks{1,task_i}=@()matlabFunction(dJ2_q5,'file','grad/dJ2_q5','vars',{[q1,q2,q3,q4,q5,q6]}); task_i=task_i+1;
tasks{1,task_i}=@()matlabFunction(dJ2_q6,'file','grad/dJ2_q6','vars',{[q1,q2,q3,q4,q5,q6]}); task_i=task_i+1;



save('model','model');
parfor i=1:length(tasks)
    tasks{1,i}();
end


%% hip velocity constraint
% J_hip = six_J_hip(q1,q2,q3,q4,q5,q6);
hipvel = J_hip*[qd1;qd2;qd3;qd4;qd5;qd6];
hip_x_vel = hipvel(1,1);
dHip_x_vel = [diff(hip_x_vel,q1);diff(hip_x_vel,q2);diff(hip_x_vel,q3);diff(hip_x_vel,q4);diff(hip_x_vel,q5);diff(hip_x_vel,q6);
              diff(hip_x_vel,qd1);diff(hip_x_vel,qd2);diff(hip_x_vel,qd3);diff(hip_x_vel,qd4);diff(hip_x_vel,qd5);diff(hip_x_vel,qd6)];  
matlabFunction(hip_x_vel,'file','posCons/hip_x_vel','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6]});
matlabFunction(dHip_x_vel,'file','posCons/dHip_x_vel','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6]});


%% below here we generate the equations for solver that has Fx,Fy as variables

syms Fx_heel Fy_heel Fx_toe Fy_toe s_toe s_heel
toe_vel = J(1:2,:)*[qd1;qd2;qd3;qd4;qd5;qd6];
y_vel_toe = toe_vel(2,1);
x_vel_toe = toe_vel(1,1);
heel_vel = J2(1:2,:)*[qd1;qd2;qd3;qd4;qd5;qd6];
y_vel_heel = heel_vel(2,1);
x_vel_heel = heel_vel(1,1);
%% constraint1: Fy = k*(th-ypos)^ks-cmax*d/dmax*yvel % (if touches the ground)
% sigma_Fn_toe = 0.5*tanh(400*(th-toe_y_pos))+0.5;
% sigma_Fn_heel = 0.5*tanh(400*(th-heel_y_pos))+0.5;

Fn_toe_c2 = Fy_toe-s_toe*   (k*(th-toe_y_pos)^ks -(cmax*(0.5*tanh(2*(th-toe_y_pos)/dmax)+0.5))*y_vel_toe);
Fn_heel_c2 = Fy_heel-s_heel*(k*(th-heel_y_pos)^ks-(cmax*(0.5*tanh(2*(th-heel_y_pos)/dmax)+0.5))*y_vel_heel);

tasks{1,task_i}=@()matlabFunction(Fn_toe_c2,'file','grf/Fn_toe_c2','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],[Fx_toe,Fy_toe],s_toe,th,k,cmax,dmax}); task_i=task_i+1;
tasks{1,task_i}=@()matlabFunction(Fn_heel_c2,'file','grf/Fn_heel_c2','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],[Fx_heel,Fy_heel],s_heel,th,k,cmax,dmax}); task_i=task_i+1;


dFn_toe_c2 = [diff(Fn_toe_c2,q1);
                diff(Fn_toe_c2,q2);
                diff(Fn_toe_c2,q3);
                diff(Fn_toe_c2,q4);
                diff(Fn_toe_c2,q5);
                diff(Fn_toe_c2,q6);
                diff(Fn_toe_c2,qd1);
                diff(Fn_toe_c2,qd2);
                diff(Fn_toe_c2,qd3);
                diff(Fn_toe_c2,qd4);
                diff(Fn_toe_c2,qd5);
                diff(Fn_toe_c2,qd6);
                zeros(6,1);
                diff(Fn_toe_c2,Fx_toe);
                diff(Fn_toe_c2,Fy_toe);
                zeros(2,1)
                diff(Fn_toe_c2,s_toe);
                0];

dFn_heel_c2 = [diff(Fn_heel_c2,q1);
                diff(Fn_heel_c2,q2);
                diff(Fn_heel_c2,q3);
                diff(Fn_heel_c2,q4);
                diff(Fn_heel_c2,q5);
                diff(Fn_heel_c2,q6);
                diff(Fn_heel_c2,qd1);
                diff(Fn_heel_c2,qd2);
                diff(Fn_heel_c2,qd3);
                diff(Fn_heel_c2,qd4);
                diff(Fn_heel_c2,qd5);
                diff(Fn_heel_c2,qd6);
                zeros(6,1);
                zeros(2,1);
                diff(Fn_heel_c2,Fx_heel);
                diff(Fn_heel_c2,Fy_heel);
                0;
                diff(Fn_heel_c2,s_heel)];
tasks{1,task_i}=@()matlabFunction(dFn_toe_c2,'file','grf/dFn_toe_c2','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],[Fx_toe,Fy_toe],s_toe,th,k,cmax,dmax}); task_i=task_i+1;
tasks{1,task_i}=@()matlabFunction(dFn_heel_c2,'file','grf/dFn_heel_c2','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],[Fx_heel,Fy_heel],s_heel,th,k,cmax,dmax}); task_i=task_i+1;
            
% Fs < uFn

Fs_toe_c1 = Fx_toe^2-(us-(us-ud)*((0.5*tanh(100*(x_vel_toe-0.05))+0.5)-(0.5*tanh(100*(x_vel_toe+0.05))-0.5)))^2-Fy_toe^2;

Fs_heel_c1 = Fx_heel^2-(us-(us-ud)*((0.5*tanh(100*(x_vel_heel-0.05))+0.5)-(0.5*tanh(100*(x_vel_heel+0.05))-0.5)))^2-Fy_heel^2;


tasks{1,task_i}=@()matlabFunction(Fs_toe_c1,'file','grf/Fs_toe_c1','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],[Fx_toe,Fy_toe],s_toe,us,ud}); task_i=task_i+1;
tasks{1,task_i}=@()matlabFunction(Fs_heel_c1,'file','grf/Fs_heel_c1','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],[Fx_heel,Fy_heel],s_heel,us,ud}); task_i=task_i+1;

dFs_toe_c1_1= diff(Fs_toe_c1,q1);
dFs_toe_c1_2= diff(Fs_toe_c1,q2);
dFs_toe_c1_3= diff(Fs_toe_c1,q3);
dFs_toe_c1_4= diff(Fs_toe_c1,q4);
dFs_toe_c1_5= diff(Fs_toe_c1,q5);
dFs_toe_c1_6= diff(Fs_toe_c1,q6);
dFs_toe_c1_7= diff(Fs_toe_c1,qd1);
dFs_toe_c1_8= diff(Fs_toe_c1,qd2);
dFs_toe_c1_9= diff(Fs_toe_c1,qd3);
dFs_toe_c1_10=diff(Fs_toe_c1,qd4);
dFs_toe_c1_11=diff(Fs_toe_c1,qd5);
dFs_toe_c1_12=diff(Fs_toe_c1,qd6);
              
dFs_toe_c1_19=diff(Fs_toe_c1,Fx_toe);
dFs_toe_c1_20=diff(Fs_toe_c1,Fy_toe);

dFs_toe_c1_23=diff(Fs_toe_c1,s_toe);
              

dFs_heel_c1_1= diff(Fs_heel_c1,q1);
dFs_heel_c1_2= diff(Fs_heel_c1,q2);
dFs_heel_c1_3= diff(Fs_heel_c1,q3);
dFs_heel_c1_4= diff(Fs_heel_c1,q4);
dFs_heel_c1_5= diff(Fs_heel_c1,q5);
dFs_heel_c1_6= diff(Fs_heel_c1,q6);
dFs_heel_c1_7= diff(Fs_heel_c1,qd1);
dFs_heel_c1_8= diff(Fs_heel_c1,qd2);
dFs_heel_c1_9= diff(Fs_heel_c1,qd3);
dFs_heel_c1_10=diff(Fs_heel_c1,qd4);
dFs_heel_c1_11=diff(Fs_heel_c1,qd5);
dFs_heel_c1_12=diff(Fs_heel_c1,qd6);
         
dFs_heel_c1_21=diff(Fs_heel_c1,Fx_heel);
dFs_heel_c1_22=diff(Fs_heel_c1,Fy_heel);
dFs_heel_c1_24=diff(Fs_heel_c1,s_heel);
tasks{1,task_i}=@()matlabFunction(dFs_toe_c1_1,'file','grf/dFs_toe_c1_1','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],[Fx_toe,Fy_toe],s_toe,us,ud}); task_i=task_i+1;
tasks{1,task_i}=@()matlabFunction(dFs_toe_c1_2,'file','grf/dFs_toe_c1_2','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],[Fx_toe,Fy_toe],s_toe,us,ud}); task_i=task_i+1;
tasks{1,task_i}=@()matlabFunction(dFs_toe_c1_3,'file','grf/dFs_toe_c1_3','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],[Fx_toe,Fy_toe],s_toe,us,ud}); task_i=task_i+1;
tasks{1,task_i}=@()matlabFunction(dFs_toe_c1_4,'file','grf/dFs_toe_c1_4','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],[Fx_toe,Fy_toe],s_toe,us,ud}); task_i=task_i+1;
tasks{1,task_i}=@()matlabFunction(dFs_toe_c1_5,'file','grf/dFs_toe_c1_5','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],[Fx_toe,Fy_toe],s_toe,us,ud}); task_i=task_i+1;
tasks{1,task_i}=@()matlabFunction(dFs_toe_c1_6,'file','grf/dFs_toe_c1_6','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],[Fx_toe,Fy_toe],s_toe,us,ud}); task_i=task_i+1;
tasks{1,task_i}=@()matlabFunction(dFs_toe_c1_7,'file','grf/dFs_toe_c1_7','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],[Fx_toe,Fy_toe],s_toe,us,ud}); task_i=task_i+1;
tasks{1,task_i}=@()matlabFunction(dFs_toe_c1_8,'file','grf/dFs_toe_c1_8','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],[Fx_toe,Fy_toe],s_toe,us,ud}); task_i=task_i+1;
tasks{1,task_i}=@()matlabFunction(dFs_toe_c1_9,'file','grf/dFs_toe_c1_9','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],[Fx_toe,Fy_toe],s_toe,us,ud}); task_i=task_i+1;
tasks{1,task_i}=@()matlabFunction(dFs_toe_c1_10,'file','grf/dFs_toe_c1_10','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],[Fx_toe,Fy_toe],s_toe,us,ud}); task_i=task_i+1;
tasks{1,task_i}=@()matlabFunction(dFs_toe_c1_11,'file','grf/dFs_toe_c1_11','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],[Fx_toe,Fy_toe],s_toe,us,ud}); task_i=task_i+1;
tasks{1,task_i}=@()matlabFunction(dFs_toe_c1_12,'file','grf/dFs_toe_c1_12','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],[Fx_toe,Fy_toe],s_toe,us,ud}); task_i=task_i+1;
tasks{1,task_i}=@()matlabFunction(dFs_toe_c1_19,'file','grf/dFs_toe_c1_19','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],[Fx_toe,Fy_toe],s_toe,us,ud}); task_i=task_i+1;
tasks{1,task_i}=@()matlabFunction(dFs_toe_c1_20,'file','grf/dFs_toe_c1_20','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],[Fx_toe,Fy_toe],s_toe,us,ud}); task_i=task_i+1;
tasks{1,task_i}=@()matlabFunction(dFs_toe_c1_23,'file','grf/dFs_toe_c1_23','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],[Fx_toe,Fy_toe],s_toe,us,ud}); task_i=task_i+1;

tasks{1,task_i}=@()matlabFunction(dFs_heel_c1_1,'file','grf/dFs_heel_c1_1','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],[Fx_heel,Fy_heel],s_heel,us,ud}); task_i=task_i+1;
tasks{1,task_i}=@()matlabFunction(dFs_heel_c1_2,'file','grf/dFs_heel_c1_2','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],[Fx_heel,Fy_heel],s_heel,us,ud}); task_i=task_i+1;
tasks{1,task_i}=@()matlabFunction(dFs_heel_c1_3,'file','grf/dFs_heel_c1_3','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],[Fx_heel,Fy_heel],s_heel,us,ud}); task_i=task_i+1;
tasks{1,task_i}=@()matlabFunction(dFs_heel_c1_4,'file','grf/dFs_heel_c1_4','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],[Fx_heel,Fy_heel],s_heel,us,ud}); task_i=task_i+1;
tasks{1,task_i}=@()matlabFunction(dFs_heel_c1_5,'file','grf/dFs_heel_c1_5','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],[Fx_heel,Fy_heel],s_heel,us,ud}); task_i=task_i+1;
tasks{1,task_i}=@()matlabFunction(dFs_heel_c1_6,'file','grf/dFs_heel_c1_6','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],[Fx_heel,Fy_heel],s_heel,us,ud}); task_i=task_i+1;
tasks{1,task_i}=@()matlabFunction(dFs_heel_c1_7,'file','grf/dFs_heel_c1_7','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],[Fx_heel,Fy_heel],s_heel,us,ud}); task_i=task_i+1;
tasks{1,task_i}=@()matlabFunction(dFs_heel_c1_8,'file','grf/dFs_heel_c1_8','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],[Fx_heel,Fy_heel],s_heel,us,ud}); task_i=task_i+1;
tasks{1,task_i}=@()matlabFunction(dFs_heel_c1_9,'file','grf/dFs_heel_c1_9','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],[Fx_heel,Fy_heel],s_heel,us,ud}); task_i=task_i+1;
tasks{1,task_i}=@()matlabFunction(dFs_heel_c1_10,'file','grf/dFs_heel_c1_10','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],[Fx_heel,Fy_heel],s_heel,us,ud}); task_i=task_i+1;
tasks{1,task_i}=@()matlabFunction(dFs_heel_c1_11,'file','grf/dFs_heel_c1_11','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],[Fx_heel,Fy_heel],s_heel,us,ud}); task_i=task_i+1;
tasks{1,task_i}=@()matlabFunction(dFs_heel_c1_12,'file','grf/dFs_heel_c1_12','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],[Fx_heel,Fy_heel],s_heel,us,ud}); task_i=task_i+1;
tasks{1,task_i}=@()matlabFunction(dFs_heel_c1_21,'file','grf/dFs_heel_c1_21','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],[Fx_heel,Fy_heel],s_heel,us,ud}); task_i=task_i+1;
tasks{1,task_i}=@()matlabFunction(dFs_heel_c1_22,'file','grf/dFs_heel_c1_22','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],[Fx_heel,Fy_heel],s_heel,us,ud}); task_i=task_i+1;
tasks{1,task_i}=@()matlabFunction(dFs_heel_c1_24,'file','grf/dFs_heel_c1_24','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],[Fx_heel,Fy_heel],s_heel,us,ud}); task_i=task_i+1;



% constraint for no Fn above toe_th

Fn_toe_c3 = s_toe*Fy_toe;
Fn_heel_c3 = s_heel*Fy_heel;

tasks{1,task_i}=@()matlabFunction(Fn_toe_c3,'file','grf/Fn_toe_c3','vars',{Fy_toe,s_toe}); task_i=task_i+1;
tasks{1,task_i}=@()matlabFunction(Fn_heel_c3,'file','grf/Fn_heel_c3','vars',{Fy_heel,s_heel}); task_i=task_i+1;

dFn_toe_c3 = [zeros(19,1);
              diff(Fn_toe_c3,Fy_toe);
              zeros(2,1);
              diff(Fn_toe_c3,s_toe);
              0];
          
dFn_heel_c3 = [zeros(21,1)
               diff(Fn_toe_c3,Fy_heel);
               0;
               diff(Fn_heel_c3,s_heel)];           
                
tasks{1,task_i}=@()matlabFunction(dFn_toe_c3,'file','grf/dFn_toe_c3','vars',{Fy_toe,s_toe}); task_i=task_i+1;
tasks{1,task_i}=@()matlabFunction(dFn_heel_c3,'file','grf/dFn_heel_c3','vars',{Fy_heel,s_heel}); task_i=task_i+1;             


% constraint for s_toe, s_heel

s_toe_con = s_toe*(toe_y_pos-th);
s_heel_con = s_heel*(heel_y_pos-th);

ds_toe_con = [diff(s_toe_con,q1);
              diff(s_toe_con,q2);
              diff(s_toe_con,q3);
              diff(s_toe_con,q4);
              diff(s_toe_con,q5);
              diff(s_toe_con,q6);
              zeros(16,1);
              diff(s_toe_con,s_toe);
              diff(s_heel_con,s_heel)];
ds_heel_con = [diff(s_heel_con,q1);
              diff(s_heel_con,q2);
              diff(s_heel_con,q3);
              diff(s_heel_con,q4);
              diff(s_heel_con,q5);
              diff(s_heel_con,q6);
              zeros(16,1);
              diff(s_toe_con,s_toe);
              diff(s_heel_con,s_heel)];         
          
tasks{1,task_i}=@()matlabFunction(s_toe_con,'file','grf/s_toe_con','vars',{[q1,q2,q3,q4,q5,q6],s_toe,th}); task_i=task_i+1;
tasks{1,task_i}=@()matlabFunction(s_heel_con,'file','grf/s_heel_con','vars',{[q1,q2,q3,q4,q5,q6],s_heel,th}); task_i=task_i+1;             
tasks{1,task_i}=@()matlabFunction(ds_toe_con,'file','grf/ds_toe_con','vars',{[q1,q2,q3,q4,q5,q6],s_toe,th}); task_i=task_i+1;
tasks{1,task_i}=@()matlabFunction(ds_heel_con,'file','grf/ds_heel_con','vars',{[q1,q2,q3,q4,q5,q6],s_heel,th}); task_i=task_i+1;  


% fs must be on the opposite direction of x_vel
Fs_toe_c2 = Fx_toe*x_vel_toe;
Fs_heel_c2 = Fx_heel*x_vel_heel;

dFs_toe_c2 = [diff(Fs_toe_c2,q1);
              diff(Fs_toe_c2,q2);
              diff(Fs_toe_c2,q3);
              diff(Fs_toe_c2,q4);
              diff(Fs_toe_c2,q5);
              diff(Fs_toe_c2,q6);
              diff(Fs_toe_c2,qd1);
              diff(Fs_toe_c2,qd2);
              diff(Fs_toe_c2,qd3);
              diff(Fs_toe_c2,qd4);
              diff(Fs_toe_c2,qd5);
              diff(Fs_toe_c2,qd6);
              zeros(6,1);
              diff(Fs_toe_c2,Fx_toe);
              zeros(5,1)];
dFs_heel_c2 = [diff(Fs_heel_c2,q1);
               diff(Fs_heel_c2,q2);
               diff(Fs_heel_c2,q3);
               diff(Fs_heel_c2,q4);
               diff(Fs_heel_c2,q5);
               diff(Fs_heel_c2,q6);
               diff(Fs_heel_c2,qd1);
               diff(Fs_heel_c2,qd2);
               diff(Fs_heel_c2,qd3);
               diff(Fs_heel_c2,qd4);
               diff(Fs_heel_c2,qd5);
               diff(Fs_heel_c2,qd6);
               zeros(8,1);
               diff(Fs_heel_c2,Fx_heel);
               zeros(3,1)];
           
tasks{1,task_i}=@()matlabFunction(Fs_toe_c2,'file','grf/Fs_toe_c2','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],Fx_toe}); task_i=task_i+1;           
tasks{1,task_i}=@()matlabFunction(Fs_heel_c2,'file','grf/Fs_heel_c2','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],Fx_heel}); task_i=task_i+1; 
tasks{1,task_i}=@()matlabFunction(dFs_toe_c2,'file','grf/dFs_toe_c2','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],Fx_toe}); task_i=task_i+1;           
tasks{1,task_i}=@()matlabFunction(dFs_heel_c2,'file','grf/dFs_heel_c2','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],Fx_heel}); task_i=task_i+1;               
               
% friction term
% no x_vel when there is external force
fri_toe = (s_toe*x_vel_toe)^2;
fri_heel = (s_heel*x_vel_heel)^2;

dfri_toe = [diff(fri_toe,q1);
            diff(fri_toe,q2);
            diff(fri_toe,q3);
            diff(fri_toe,q4);
            diff(fri_toe,q5);
            diff(fri_toe,q6);
            diff(fri_toe,qd1);
            diff(fri_toe,qd2);
            diff(fri_toe,qd3);
            diff(fri_toe,qd4);
            diff(fri_toe,qd5);
            diff(fri_toe,qd6);
            zeros(10,1);
            diff(fri_toe,s_toe)
            0];
dfri_heel = [diff(fri_heel,q1);
             diff(fri_heel,q2);
             diff(fri_heel,q3);
             diff(fri_heel,q4);
             diff(fri_heel,q5);
             diff(fri_heel,q6);
             diff(fri_heel,qd1);
             diff(fri_heel,qd2);
             diff(fri_heel,qd3);
             diff(fri_heel,qd4);
             diff(fri_heel,qd5);
             diff(fri_heel,qd6);
             zeros(11,1);
             diff(fri_heel,s_heel)];
             
tasks{1,task_i}=@()matlabFunction(fri_toe,'file','grf/fri_toe','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],s_toe}); task_i=task_i+1;           
tasks{1,task_i}=@()matlabFunction(fri_heel,'file','grf/fri_heel','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],s_heel}); task_i=task_i+1; 
tasks{1,task_i}=@()matlabFunction(dfri_toe,'file','grf/dfri_toe','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],s_toe}); task_i=task_i+1;           
tasks{1,task_i}=@()matlabFunction(dfri_heel,'file','grf/dfri_heel','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],s_heel}); task_i=task_i+1;    




parfor i=1:length(tasks)
    tasks{1,i}();
end         
          
