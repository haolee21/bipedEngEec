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
m_thight = 0.1416*totM;
m_head = 0.0694*totM;
m_trunk = 0.4346*totM;
m_torso = m_trunk+m_head; %this total mass is not 100% though

model.totM = totM;
model.l_heel = l_heel;
model.l_foot = l_foot;
save('model','model');
%CoM pos
lc_thigh2 = 0.433*l_thigh;
lc_calf2 = 0.433*l_calf;
lc_calf1 = l_calf-lc_calf2;
lc_thigh1 = l_thigh-lc_thigh2;

lc_foot = 0.4415*l_foot; % de Leva 
% include head when calculating, suppose neck is fixed
lc_torso = (0.4486*l_torso*m_torso+(0.1166*totH*(1-0.5976)+l_torso)*m_head)/(m_head+m_trunk); 




I_calf = zeros(3);
% I_calf1(3,3) = m_calf*l_calf^2/12;



I_thight = zeros(3);
% I_thight(3,3) = m_thight*l_thight^2/12;
I_torso = zeros(3);
% I_torso(3,3) = m_torso*l_torso^2/12;
I_foot = zeros(3);


% we need extra struct to storage link length, since for the feet, we need
% heel and foot long
end_eff = [l_foot,l_heel,l_calf];


robot.links(1).m = m_calf;
robot.links(2).m = m_thight;
robot.links(3).m = m_torso;
robot.links(4).m = m_thight;
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



syms q1 q2 q3 q4 q5 q6 th qd1 qd2 qd3 qd4 qd5 qd6 th;
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
heelPos = turnRTtoMatrix(robot.A([1,2,3,4,5],[q1 q2 q3 q4 q5]))*[l_calf+l_heel,0,0,1].'; %we suppose heel is on calf, doesn't rotate
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


%% activation function (sigma)
sigma_toe = (0.5*tanh(400*(th-end_y_pos))+0.5);
matlabFunction(sigma_toe,'file','sigma_toe','vars',{[q1,q2,q3,q4,q5,q6],th});
dsigma_toe = [diff(sigma_toe,q1);diff(sigma_toe,q2);diff(sigma_toe,q3);diff(sigma_toe,q4);diff(sigma_toe,q5);diff(sigma_toe,q6)];
matlabFunction(dsigma_toe,'file','dsigma_toe','vars',{[q1,q2,q3,q4,q5,q6],th});

sigma_heel =  (0.5*tanh(400*(th-heel_y_pos))+0.5);
matlabFunction(sigma_heel,'file','sigma_heel','vars',{[q1,q2,q3,q4,q5,q6],th});
dsigma_heel = [diff(sigma_heel,q1);diff(sigma_heel,q2);diff(sigma_heel,q3);diff(sigma_heel,q4);diff(sigma_heel,q5);diff(sigma_heel,q6)];
matlabFunction(dsigma_heel,'file','dsigma_heel','vars',{[q1,q2,q3,q4,q5,q6],th});


%% generate the joint/CoM pos (base frame) for graphing
% Joint
knee_front = turnRTtoMatrix(robot.A(1,q1))*[l_calf,0,0,1].';
hip_front = turnRTtoMatrix(robot.A([1,2],[q1,q2]))*[l_thigh,0,0,1].';
head = turnRTtoMatrix(robot.A([1,2,3],[q1,q2,q3]))*[l_torso,0,0,1].';
knee_back = turnRTtoMatrix(robot.A([1,2,3,4],[q1,q2,q3,q4]))*[l_thigh,0,0,1].';
ankle_back =turnRTtoMatrix(robot.A([1,2,3,4,5],[q1,q2,q3,q4,q5]))*[l_calf,0,0,1].';
toe_back = turnRTtoMatrix(robot.A([1,2,3,4,5,6],[q1,q2,q3,q4,q5,q6]))*[l_foot,l_heel,0,1].';
heel_back = turnRTtoMatrix(robot.A([1,2,3,4,5],[q1,q2,q3,q4,q5]))*[l_calf+l_heel,0,0,1].';
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
tasks = cell(1,34);
tasks{1,1}=@()matlabFunction(dyn.M,'File','dyn/six_M','vars',[q2,q3,q4,q5,q6]);
tasks{1,2}=@()matlabFunction(dyn.G,'File','dyn/six_G','vars',[q1 q2 q3 q4 q5 q6]);
tasks{1,3}=@()matlabFunction(dyn.J,'File','dyn/six_J','vars',[q1 q2 q3 q4 q5 q6]);
tasks{1,4}=@()matlabFunction(dyn.V,'File','dyn/six_V','vars',[q2 q3 q4 q5 q6 qd1 qd2 qd3 qd4 qd5 qd6]);
tasks{1,5}=@()matlabFunction(dyn.J2,'file','dyn/six_J2','vars',[q1,q2,q3,q4,q5,q6]);
% 
% 
% %% generate gradient for the M,G,J,V
G=dyn.G;
M=dyn.M;
V=dyn.V;
J=dyn.J;
J2 = dyn.J2; %jacobian on heel

dG_dx = [diff(G,q1);diff(G,q2);diff(G,q3);diff(G,q4);diff(G,q5);diff(G,q6);zeros(2*numJ,numJ)];

tasks{1,6}=@()matlabFunction(dG_dx,'file','grad/dG_dx','vars',[q1,q2,q3,q4,q5,q6]);       



dV_dx = [diff(V,q1);diff(V,q2);diff(V,q3);diff(V,q4);diff(V,q5);diff(V,q6);diff(V,qd1);diff(V,qd2);diff(V,qd3);diff(V,qd4);diff(V,qd5);diff(V,qd6);zeros(numJ,numJ)];

tasks{1,7}=@()matlabFunction(dV_dx,'file','grad/dV_dx','vars',[qd1,qd2,qd3,qd4,qd5,qd6,q2,q3,q4,q5,q6]);


dM_dx2 = diff(M,q2);
dM_dx3 = diff(M,q3);
dM_dx4 = diff(M,q4);
dM_dx5 = diff(M,q5);
dM_dx6 = diff(M,q6);
tasks{1,8}=@()matlabFunction(dM_dx2,'file','grad/dM_dx2','vars',[q2,q3,q4,q5,q6]);
tasks{1,9}=@()matlabFunction(dM_dx3,'file','grad/dM_dx3','vars',[q2,q3,q4,q5,q6]);
tasks{1,10}=@()matlabFunction(dM_dx4,'file','grad/dM_dx4','vars',[q2,q3,q4,q5,q6]);
tasks{1,11}=@()matlabFunction(dM_dx5,'file','grad/dM_dx5','vars',[q2,q3,q4,q5,q6]);
tasks{1,12}=@()matlabFunction(dM_dx6,'file','grad/dM_dx6','vars',[q2,q3,q4,q5,q6]);        

parfor i=1:length(tasks)
    tasks{1,i}();
end
%% generate the beta function for GRF

% f_ext on toe
% J = dyn.J;
% J2 =dyn.J2;

J_f_toe = J(1:2,:).';

tasks = cell(1,32);



Mext_toe = (J_f_toe.'*J_f_toe)\J_f_toe.';
beta_grf_toe = sigma_toe*J_f_toe;

dMext_toe_dx1 = diff(Mext_toe,q1);
dMext_toe_dx2 = diff(Mext_toe,q2);
dMext_toe_dx3 = diff(Mext_toe,q3);
dMext_toe_dx4 = diff(Mext_toe,q4);
dMext_toe_dx5 = diff(Mext_toe,q5);
dMext_toe_dx6 = diff(Mext_toe,q6);


tasks{1,1}= @()matlabFunction(Mext_toe,'file','grf/Mext_toe','vars',{[q1,q2,q3,q4,q5,q6]});
tasks{1,2}= @()matlabFunction(dMext_toe_dx1,'file','grf/dMext_toe_dx1','vars',{[q1,q2,q3,q4,q5,q6]});
tasks{1,3}= @()matlabFunction(dMext_toe_dx2,'file','grf/dMext_toe_dx2','vars',{[q1,q2,q3,q4,q5,q6]});
tasks{1,4}= @()matlabFunction(dMext_toe_dx3,'file','grf/dMext_toe_dx3','vars',{[q1,q2,q3,q4,q5,q6]});
tasks{1,5}= @()matlabFunction(dMext_toe_dx4,'file','grf/dMext_toe_dx4','vars',{[q1,q2,q3,q4,q5,q6]});
tasks{1,6}= @()matlabFunction(dMext_toe_dx5,'file','grf/dMext_toe_dx5','vars',{[q1,q2,q3,q4,q5,q6]});
tasks{1,7}= @()matlabFunction(dMext_toe_dx6,'file','grf/dMext_toe_dx6','vars',{[q1,q2,q3,q4,q5,q6]});




% beta_grf = sigma*J_f/(J_f.'*J_f)*J_f.';
dbeta_toe_dx1 = diff(beta_grf_toe,q1);
dbeta_toe_dx2 = diff(beta_grf_toe,q2);
dbeta_toe_dx3 = diff(beta_grf_toe,q3);
dbeta_toe_dx4 = diff(beta_grf_toe,q4);
dbeta_toe_dx5 = diff(beta_grf_toe,q5);
dbeta_toe_dx6 = diff(beta_grf_toe,q6);

tasks{1,8}=@()matlabFunction(beta_grf_toe,'file','grf/beta_grf_toe','vars',{[q1,q2,q3,q4,q5,q6],th});
tasks{1,9}=@()matlabFunction(dbeta_toe_dx1,'file','grf/dbeta_toe_dx1','vars',{[q1,q2,q3,q4,q5,q6],th});
tasks{1,10}=@()matlabFunction(dbeta_toe_dx2,'file','grf/dbeta_toe_dx2','vars',{[q1,q2,q3,q4,q5,q6],th});
tasks{1,11}=@()matlabFunction(dbeta_toe_dx3,'file','grf/dbeta_toe_dx3','vars',{[q1,q2,q3,q4,q5,q6],th});
tasks{1,12}=@()matlabFunction(dbeta_toe_dx4,'file','grf/dbeta_toe_dx4','vars',{[q1,q2,q3,q4,q5,q6],th});
tasks{1,13}=@()matlabFunction(dbeta_toe_dx5,'file','grf/dbeta_toe_dx5','vars',{[q1,q2,q3,q4,q5,q6],th});
tasks{1,14}=@()matlabFunction(dbeta_toe_dx6,'file','grf/dbeta_toe_dx6','vars',{[q1,q2,q3,q4,q5,q6],th});


% external force on heel
J_f_heel = J2(1:2,:).';
Mext_heel = (J_f_heel.'*J_f_heel)\J_f_heel.';
beta_grf_heel = sigma_heel*J_f_heel;

dMext_heel_dx1 = diff(Mext_heel,q1);
dMext_heel_dx2 = diff(Mext_heel,q2);
dMext_heel_dx3 = diff(Mext_heel,q3);
dMext_heel_dx4 = diff(Mext_heel,q4);
dMext_heel_dx5 = diff(Mext_heel,q5);
dMext_heel_dx6 = diff(Mext_heel,q6);

tasks{1,15}=@()matlabFunction(Mext_heel,'file','grf/Mext_heel','vars',{[q1,q2,q3,q4,q5,q6]});
tasks{1,16}=@()matlabFunction(dMext_heel_dx1,'file','grf/dMext_heel_dx1','vars',{[q1,q2,q3,q4,q5,q6]});
tasks{1,17}=@()matlabFunction(dMext_heel_dx2,'file','grf/dMext_heel_dx2','vars',{[q1,q2,q3,q4,q5,q6]});
tasks{1,18}=@()matlabFunction(dMext_heel_dx3,'file','grf/dMext_heel_dx3','vars',{[q1,q2,q3,q4,q5,q6]});
tasks{1,19}=@()matlabFunction(dMext_heel_dx4,'file','grf/dMext_heel_dx4','vars',{[q1,q2,q3,q4,q5,q6]});
tasks{1,20}=@()matlabFunction(dMext_heel_dx5,'file','grf/dMext_heel_dx5','vars',{[q1,q2,q3,q4,q5,q6]});
tasks{1,21}=@()matlabFunction(dMext_heel_dx6,'file','grf/dMext_heel_dx6','vars',{[q1,q2,q3,q4,q5,q6]});
% beta_grf = sigma*J_f/(J_f.'*J_f)*J_f.';
dbeta_heel_dx1 = diff(beta_grf_heel,q1);
dbeta_heel_dx2 = diff(beta_grf_heel,q2);
dbeta_heel_dx3 = diff(beta_grf_heel,q3);
dbeta_heel_dx4 = diff(beta_grf_heel,q4);
dbeta_heel_dx5 = diff(beta_grf_heel,q5);
dbeta_heel_dx6 = diff(beta_grf_heel,q6);

tasks{1,22}=@()matlabFunction(beta_grf_heel,'file','grf/beta_grf_heel','vars',{[q1,q2,q3,q4,q5,q6],th});
tasks{1,23}=@()matlabFunction(dbeta_heel_dx1,'file','grf/dbeta_heel_dx1','vars',{[q1,q2,q3,q4,q5,q6],th});
tasks{1,24}=@()matlabFunction(dbeta_heel_dx2,'file','grf/dbeta_heel_dx2','vars',{[q1,q2,q3,q4,q5,q6],th});
tasks{1,25}=@()matlabFunction(dbeta_heel_dx3,'file','grf/dbeta_heel_dx3','vars',{[q1,q2,q3,q4,q5,q6],th});
tasks{1,26}=@()matlabFunction(dbeta_heel_dx4,'file','grf/dbeta_heel_dx4','vars',{[q1,q2,q3,q4,q5,q6],th});
tasks{1,27}=@()matlabFunction(dbeta_heel_dx5,'file','grf/dbeta_heel_dx5','vars',{[q1,q2,q3,q4,q5,q6],th});
tasks{1,28}=@()matlabFunction(dbeta_heel_dx6,'file','grf/dbeta_heel_dx6','vars',{[q1,q2,q3,q4,q5,q6],th});

%% generate friction for object function 
% we only consider the horizontal friction, so we use x_dot = J(1,:)*q_dot

%toe
x_vel1 = J(1,:)*[qd1;qd2;qd3;qd4;qd5;qd6];
fri1 = sigma_toe*(x_vel1*x_vel1);
fri1_dx = [diff(fri1,q1);
          diff(fri1,q2);
          diff(fri1,q3);
          diff(fri1,q4);
          diff(fri1,q5);
          diff(fri1,q6);
          diff(fri1,qd1);
          diff(fri1,qd2);
          diff(fri1,qd3);
          diff(fri1,qd4);
          diff(fri1,qd5);
          diff(fri1,qd6);
          zeros(numJ,1)];
tasks{1,29}=@()matlabFunction(fri1,'file','obj/fri_toe','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],th});
tasks{1,30}=@()matlabFunction(fri1_dx,'file','obj/fri_toe_dx','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],th});

% ankle

x_vel2 = J2(1,:)*[qd1;qd2;qd3;qd4;qd5;qd6];
fri2 = sigma_heel*(x_vel2*x_vel2);
fri2_dx = [diff(fri2,q1);
          diff(fri2,q2);
          diff(fri2,q3);
          diff(fri2,q4);
          diff(fri2,q5);
          diff(fri2,q6);
          diff(fri2,qd1);
          diff(fri2,qd2);
          diff(fri2,qd3);
          diff(fri2,qd4);
          diff(fri2,qd5);
          diff(fri2,qd6);
          zeros(numJ,1)];
tasks{1,31}=@()matlabFunction(fri2,'file','obj/fri_heel','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],th});
tasks{1,32}=@()matlabFunction(fri2_dx,'file','obj/fri_heel_dx','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],th});


%% front/back knee locking
% we will add a torsional spring at the front knee
% it is triggered when q2<1 deg

syms kne_stop ank_stop qkne qank real ;
sigma_knee = 0.5*tanh(250*(kne_stop-qkne))+0.5;
dsigma_knee = diff(sigma_knee,qkne);
sigma_ank = 0.5*tanh(250*(ank_stop-qank))+0.5;
dsigma_ank = diff(sigma_ank,qank);

task{1,33}=@()matlabFunction(sigma_knee,'file','knee_spring/sigma_knee','vars',[qkne,kne_stop]);
task{1,34}=@()matlabFunction(dsigma_knee,'file','knee_spring/dsigma_knee','vars',[qkne,kne_stop]);
task{1,35}=@()matlabFunction(sigma_ank,'file','knee_spring/sigma_ank','vars',[qank,ank_stop]);
task{1,36}=@()matlabFunction(dsigma_ank,'file','knee_spring/dsigma_ank','vars',[qank,ank_stop]);




%% GRF 
% since for grf we need x_vel, we generate it at last


k = 1e7;
ks = 2.2; %use ks to replace the spring constant e in the paper
cmax = 1500;

dmax = 1e-2;
us = 0.8;
ud=0.7;

toe_y_pos = end_y_pos;

model.dmax = dmax;
% here we introduce the nondifferentiable term into the model
% while it is not differentiable, it is at least continuous (relu like
% function)

% toe

toe_vel = J(1:2,:)*[qd1;qd2;qd3;qd4;qd5;qd6];
y_vel_toe = toe_vel(2,1);
x_vel_toe = toe_vel(1,1);
Fn_toe1 = k*(th-toe_y_pos)^ks-(th-toe_y_pos)/dmax*cmax*y_vel_toe;%y_vel_toe<0 when Fn>0
Fn_toe2 = k*(th-toe_y_pos)^ks-cmax*y_vel_toe;
dFn_toe1 = [diff(Fn_toe1,q1);diff(Fn_toe1,q2);diff(Fn_toe1,q3);diff(Fn_toe1,q4);diff(Fn_toe1,q5);diff(Fn_toe1,q6);...
            diff(Fn_toe1,qd1);diff(Fn_toe1,qd2);diff(Fn_toe1,qd3);diff(Fn_toe1,qd4);diff(Fn_toe1,qd5);diff(Fn_toe1,qd6)];
dFn_toe2 = [diff(Fn_toe2,q1);diff(Fn_toe2,q2);diff(Fn_toe2,q3);diff(Fn_toe2,q4);diff(Fn_toe2,q5);diff(Fn_toe2,q6);...
            diff(Fn_toe2,qd1);diff(Fn_toe2,qd2);diff(Fn_toe2,qd3);diff(Fn_toe2,qd4);diff(Fn_toe2,qd5);diff(Fn_toe2,qd6)];

Fs_toe1_s = us*Fn_toe1;
Fs_toe1_d = ud*Fn_toe1;

Fs_toe2_s = us*Fn_toe2;
Fs_toe2_d = ud*Fn_toe2;
        
dFs_toe1_s = [diff(Fs_toe1_s,q1);diff(Fs_toe1_s,q2);diff(Fs_toe1_s,q3);diff(Fs_toe1_s,q4);diff(Fs_toe1_s,q5);diff(Fs_toe1_s,q6);...
              diff(Fs_toe1_s,qd1);diff(Fs_toe1_s,qd2);diff(Fs_toe1_s,qd3);diff(Fs_toe1_s,qd4);diff(Fs_toe1_s,qd5);diff(Fs_toe1_s,qd6)];

dFs_toe1_d = [diff(Fs_toe1_d,q1);diff(Fs_toe1_d,q2);diff(Fs_toe1_d,q3);diff(Fs_toe1_d,q4);diff(Fs_toe1_d,q5);diff(Fs_toe1_d,q6);...
              diff(Fs_toe1_d,qd1);diff(Fs_toe1_d,qd2);diff(Fs_toe1_d,qd3);diff(Fs_toe1_d,qd4);diff(Fs_toe1_d,qd5);diff(Fs_toe1_d,qd6)];
          
dFs_toe2_s = [diff(Fs_toe2_s,q1);diff(Fs_toe2_s,q2);diff(Fs_toe2_s,q3);diff(Fs_toe2_s,q4);diff(Fs_toe2_s,q5);diff(Fs_toe2_s,q6);...
              diff(Fs_toe2_s,qd1);diff(Fs_toe2_s,qd2);diff(Fs_toe2_s,qd3);diff(Fs_toe2_s,qd4);diff(Fs_toe2_s,qd5);diff(Fs_toe2_s,qd6)];    
          
dFs_toe2_d = [diff(Fs_toe2_d,q1);diff(Fs_toe2_d,q2);diff(Fs_toe2_d,q3);diff(Fs_toe2_d,q4);diff(Fs_toe2_d,q5);diff(Fs_toe2_d,q6);...
              diff(Fs_toe2_d,qd1);diff(Fs_toe2_d,qd2);diff(Fs_toe2_d,qd3);diff(Fs_toe2_d,qd4);diff(Fs_toe2_d,qd5);diff(Fs_toe2_d,qd6)];    
          
tasks{1,37} = @()matlabFunction(Fn_toe1,'file','grf/Fn_toe1','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],th}); 
tasks{1,38} = @()matlabFunction(Fn_toe2,'file','grf/Fn_toe2','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],th}); 
tasks{1,39} = @()matlabFunction(dFn_toe1,'file','grf/dFn_toe1','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],th}); 
tasks{1,40} = @()matlabFunction(dFn_toe2,'file','grf/dFn_toe2','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],th}); 
tasks{1,41} = @()matlabFunction(Fs_toe1_s,'file','grf/Fs_toe1_s','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],th}); 
tasks{1,42} = @()matlabFunction(Fs_toe1_d,'file','grf/Fs_toe1_d','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],th}); 
tasks{1,43} = @()matlabFunction(Fs_toe2_s,'file','grf/Fs_toe2_s','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],th}); 
tasks{1,44} = @()matlabFunction(Fs_toe2_d,'file','grf/Fs_toe2_d','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],th}); 

tasks{1,45} = @()matlabFunction(dFs_toe1_s,'file','grf/dFs_toe1_s','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],th}); 
tasks{1,46} = @()matlabFunction(dFs_toe1_d,'file','grf/dFs_toe1_d','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],th}); 
tasks{1,47} = @()matlabFunction(dFs_toe2_s,'file','grf/dFs_toe2_s','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],th}); 
tasks{1,48} = @()matlabFunction(dFs_toe2_d,'file','grf/dFs_toe2_d','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],th}); 
% heel
heel_vel = J(1:2,:)*[qd1;qd2;qd3;qd4;qd5;qd6];
y_vel_heel = heel_vel(2,1);
x_vel_heel = heel_vel(1,1);
Fn_heel1 = k*(th-heel_y_pos)^ks-(th-heel_y_pos)/dmax*cmax*y_vel_heel;%y_vel_heel<0 when Fn>0
Fn_heel2 = k*(th-heel_y_pos)^ks-cmax*y_vel_heel;
dFn_heel1 = [diff(Fn_heel1,q1);diff(Fn_heel1,q2);diff(Fn_heel1,q3);diff(Fn_heel1,q4);diff(Fn_heel1,q5);diff(Fn_heel1,q6);...
            diff(Fn_heel1,qd1);diff(Fn_heel1,qd2);diff(Fn_heel1,qd3);diff(Fn_heel1,qd4);diff(Fn_heel1,qd5);diff(Fn_heel1,qd6)];
dFn_heel2 = [diff(Fn_heel2,q1);diff(Fn_heel2,q2);diff(Fn_heel2,q3);diff(Fn_heel2,q4);diff(Fn_heel2,q5);diff(Fn_heel2,q6);...
            diff(Fn_heel2,qd1);diff(Fn_heel2,qd2);diff(Fn_heel2,qd3);diff(Fn_heel2,qd4);diff(Fn_heel2,qd5);diff(Fn_heel2,qd6)];

Fs_heel1_s = us*Fn_heel1;
Fs_heel1_d = ud*Fn_heel1;

Fs_heel2_s = us*Fn_heel2;
Fs_heel2_d = ud*Fn_heel2;
        
dFs_heel1_s = [diff(Fs_heel1_s,q1);diff(Fs_heel1_s,q2);diff(Fs_heel1_s,q3);diff(Fs_heel1_s,q4);diff(Fs_heel1_s,q5);diff(Fs_heel1_s,q6);...
              diff(Fs_heel1_s,qd1);diff(Fs_heel1_s,qd2);diff(Fs_heel1_s,qd3);diff(Fs_heel1_s,qd4);diff(Fs_heel1_s,qd5);diff(Fs_heel1_s,qd6)];

dFs_heel1_d = [diff(Fs_heel1_d,q1);diff(Fs_heel1_d,q2);diff(Fs_heel1_d,q3);diff(Fs_heel1_d,q4);diff(Fs_heel1_d,q5);diff(Fs_heel1_d,q6);...
              diff(Fs_heel1_d,qd1);diff(Fs_heel1_d,qd2);diff(Fs_heel1_d,qd3);diff(Fs_heel1_d,qd4);diff(Fs_heel1_d,qd5);diff(Fs_heel1_d,qd6)];
          
dFs_heel2_s = [diff(Fs_heel2_s,q1);diff(Fs_heel2_s,q2);diff(Fs_heel2_s,q3);diff(Fs_heel2_s,q4);diff(Fs_heel2_s,q5);diff(Fs_heel2_s,q6);...
              diff(Fs_heel2_s,qd1);diff(Fs_heel2_s,qd2);diff(Fs_heel2_s,qd3);diff(Fs_heel2_s,qd4);diff(Fs_heel2_s,qd5);diff(Fs_heel2_s,qd6)];    
          
dFs_heel2_d = [diff(Fs_heel2_d,q1);diff(Fs_heel2_d,q2);diff(Fs_heel2_d,q3);diff(Fs_heel2_d,q4);diff(Fs_heel2_d,q5);diff(Fs_heel2_d,q6);...
              diff(Fs_heel2_d,qd1);diff(Fs_heel2_d,qd2);diff(Fs_heel2_d,qd3);diff(Fs_heel2_d,qd4);diff(Fs_heel2_d,qd5);diff(Fs_heel2_d,qd6)]; 
          
          
tasks{1,49} = @()matlabFunction(Fn_heel1,'file','grf/Fn_heel1','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],th}); 
tasks{1,50} = @()matlabFunction(Fn_heel2,'file','grf/Fn_heel2','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],th}); 
tasks{1,51} = @()matlabFunction(dFn_heel1,'file','grf/dFn_heel1','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],th}); 
tasks{1,52} = @()matlabFunction(dFn_heel2,'file','grf/dFn_heel2','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],th}); 
tasks{1,53} = @()matlabFunction(Fs_heel1_s,'file','grf/Fs_heel1_s','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],th}); 
tasks{1,54} = @()matlabFunction(Fs_heel1_d,'file','grf/Fs_heel1_d','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],th}); 
tasks{1,55} = @()matlabFunction(Fs_heel2_s,'file','grf/Fs_heel2_s','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],th}); 
tasks{1,56} = @()matlabFunction(Fs_heel2_d,'file','grf/Fs_heel2_d','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],th}); 

tasks{1,57} = @()matlabFunction(dFs_heel1_s,'file','grf/dFs_heel1_s','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],th}); 
tasks{1,58} = @()matlabFunction(dFs_heel1_d,'file','grf/dFs_heel1_d','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],th}); 
tasks{1,59} = @()matlabFunction(dFs_heel2_s,'file','grf/dFs_heel2_s','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],th}); 
tasks{1,60} = @()matlabFunction(dFs_heel2_d,'file','grf/dFs_heel2_d','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],th}); 
% we need to take derivative of the jacobian
dJ_q1 = diff(J,q1);
dJ_q2 = diff(J,q2);
dJ_q3 = diff(J,q3);
dJ_q4 = diff(J,q4);
dJ_q5 = diff(J,q5);
dJ_q6 = diff(J,q6);

tasks{1,61}=@()matlabFunction(dJ_q1,'file','grad/dJ_q1','vars',{[q1,q2,q3,q4,q5,q6]}); 
tasks{1,62}=@()matlabFunction(dJ_q2,'file','grad/dJ_q2','vars',{[q1,q2,q3,q4,q5,q6]}); 
tasks{1,63}=@()matlabFunction(dJ_q3,'file','grad/dJ_q3','vars',{[q1,q2,q3,q4,q5,q6]}); 
tasks{1,64}=@()matlabFunction(dJ_q4,'file','grad/dJ_q4','vars',{[q1,q2,q3,q4,q5,q6]}); 
tasks{1,65}=@()matlabFunction(dJ_q5,'file','grad/dJ_q5','vars',{[q1,q2,q3,q4,q5,q6]}); 
tasks{1,66}=@()matlabFunction(dJ_q6,'file','grad/dJ_q6','vars',{[q1,q2,q3,q4,q5,q6]}); 

dJ2_q1 = diff(J2,q1);
dJ2_q2 = diff(J2,q2);
dJ2_q3 = diff(J2,q3);
dJ2_q4 = diff(J2,q4);
dJ2_q5 = diff(J2,q5);
dJ2_q6 = diff(J2,q6);      

tasks{1,67}=@()matlabFunction(dJ2_q1,'file','grad/dJ2_q1','vars',{[q1,q2,q3,q4,q5,q6]}); 
tasks{1,68}=@()matlabFunction(dJ2_q2,'file','grad/dJ2_q2','vars',{[q1,q2,q3,q4,q5,q6]}); 
tasks{1,69}=@()matlabFunction(dJ2_q3,'file','grad/dJ2_q3','vars',{[q1,q2,q3,q4,q5,q6]}); 
tasks{1,70}=@()matlabFunction(dJ2_q4,'file','grad/dJ2_q4','vars',{[q1,q2,q3,q4,q5,q6]}); 
tasks{1,71}=@()matlabFunction(dJ2_q5,'file','grad/dJ2_q5','vars',{[q1,q2,q3,q4,q5,q6]}); 
tasks{1,72}=@()matlabFunction(dJ2_q6,'file','grad/dJ2_q6','vars',{[q1,q2,q3,q4,q5,q6]}); 



save('model','model');
parfor i=1:length(tasks)
    tasks{1,i}();
end
