%%% This script will generate a 2-d biped model
%   the model has 5 segments, calf, thigh, and torso
%%% 
clear;

addpath ../

numJ = 6;
l_foot = 0.26;
l_calf = 0.43;
l_thigh =0.4;
l_torso = 0.52+0.19;
DH = [0, 0, 0, 0, 0, 0;...
      0, l_calf, 0, 0, 0, 0;...
      0, l_thigh, 0, 0, 0, 0;...
      0,  0, 0, 0, 0, 0;...
      0,l_thigh, 0, 0, 0, 0;...
      0,l_calf,0,0,0,0];
robot=DH2Robot(DH,1);
robot.gravity=[0,9.81,0];

lc_calf1 = 0.43-0.19;
lc_thigh1 = 0.4-0.17;
lc_torso = 0.29;
lc_thigh2 = 0.17;
lc_calf2 = 0.19;
lc_foot = 0.13;


%kg
m_foot = 1.05;
m_calf = 3.52;
m_thight = 7.77;
m_torso = 34.66+6.07;


I_calf = zeros(3);
% I_calf1(3,3) = m_calf*l_calf^2/12;



I_thight = zeros(3);
% I_thight(3,3) = m_thight*l_thight^2/12;
I_torso = zeros(3);
% I_torso(3,3) = m_torso*l_torso^2/12;
I_foot = zeros(3);


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
endPos = turnRTtoMatrix(robot.A([1,2,3,4,5,6],[q1 q2 q3 q4 q5 q6]))*[l_foot,0,0,1].';
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

% ankpos
ankPos = turnRTtoMatrix(robot.A([1,2,3,4,5,6],[q1 q2 q3 q4 q5 q6]))*[0,0,0,1].';
ankPos = simplify(ankPos(1:3,1));
ank_y_pos = ankPos(2,1);
matlabFunction(ank_y_pos,'file','posCons/ank_y_pos','Vars',{[q1,q2,q3,q4,q5,q6]});
ank_y_grad = [diff(ank_y_pos,q1),diff(ank_y_pos,q2),diff(ank_y_pos,q3),diff(ank_y_pos,q4),diff(ank_y_pos,q5),diff(ank_y_pos,q6)];
matlabFunction(ank_y_grad,'file','posCons/ank_y_grad','Vars',{[q1,q2,q3,q4,q5,q6]});
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
sigma = 0.8*(0.5*tanh(400*(th-end_y_pos))+0.5);
matlabFunction(sigma,'file','sigma_out','vars',[q1,q2,q3,q4,q5,q6,th]);
dsigma_dq = [diff(sigma,q1);diff(sigma,q2);diff(sigma,q3);diff(sigma,q4);diff(sigma,q5);diff(sigma,q6)];
matlabFunction(dsigma_dq,'file','dsigma_dq','vars',[q1,q2,q3,q4,q5,q6,th]);


%% generate the joint/CoM pos (base frame) for graphing
% Joint
knee_front = turnRTtoMatrix(robot.A(1,q1))*[l_calf,0,0,1].';
hip_front = turnRTtoMatrix(robot.A([1,2],[q1,q2]))*[l_thigh,0,0,1].';
head = turnRTtoMatrix(robot.A([1,2,3],[q1,q2,q3]))*[l_torso,0,0,1].';
knee_back = turnRTtoMatrix(robot.A([1,2,3,4],[q1,q2,q3,q4]))*[l_thigh,0,0,1].';
ankle_back =turnRTtoMatrix(robot.A([1,2,3,4,5],[q1,q2,q3,q4,q5]))*[l_calf,0,0,1].';
toe_back = turnRTtoMatrix(robot.A([1,2,3,4,5,6],[q1,q2,q3,q4,q5,q6]))*[l_foot,0,0,1].';
% Com
calf_front = turnRTtoMatrix(robot.A(1,q1))*[lc_calf1,0,0,1].';
thigh_front =turnRTtoMatrix(robot.A([1,2],[q1,q2]))*[lc_thigh1,0,0,1].';
torso = turnRTtoMatrix(robot.A([1,2,3],[q1,q2,q3]))*[lc_torso,0,0,1].';
thigh_back = turnRTtoMatrix(robot.A([1,2,3,4],[q1,q2,q3,q4]))*[lc_thigh2,0,0,1].';
calf_back = turnRTtoMatrix(robot.A([1,2,3,4,5],[q1,q2,q3,q4,q5]))*[lc_calf2,0,0,1].';
foot_back = turnRTtoMatrix(robot.A([1,2,3,4,5,6],[q1,q2,q3,q4,q5,q6]))*[lc_foot,0,0,1].';

draw_pos = [knee_front,hip_front,head,knee_back,ankle_back,toe_back;...
            calf_front,thigh_front,torso,thigh_back,calf_back,foot_back];
draw_pos = [draw_pos(1:3,:);draw_pos(5:7,:)];
matlabFunction(draw_pos,'file','getRobotPos','vars',[q1,q2,q3,q4,q5,q6]);
%q6 is never important since 
% 1. we directly calculate the ankle torque with jacobian
% 2. we have no feet link, the mass and inertia are 0 for that joint

dyn = dynGen(robot);


matlabFunction(dyn.M,'File','dyn/six_M','vars',[q2,q3,q4,q5,q6]);
matlabFunction(dyn.G,'File','dyn/six_G','vars',[q1 q2 q3 q4 q5 q6]);
matlabFunction(dyn.J,'File','dyn/six_J','vars',[q1 q2 q3 q4 q5 q6]);
matlabFunction(dyn.V,'File','dyn/six_V','vars',[q2 q3 q4 q5 q6 qd1 qd2 qd3 qd4 qd5 qd6]);


%% generate gradient for the M,G,J,V
G=dyn.G;
M=dyn.M;
V=dyn.V;
J=dyn.J;


dG_dx = [diff(G,q1);diff(G,q2);diff(G,q3);diff(G,q4);diff(G,q5);diff(G,q6);zeros(2*numJ,numJ)];

matlabFunction(dG_dx,'file','grad/dG_dx','vars',[q1,q2,q3,q4,q5,q6]);       



dV_dx = [diff(V,q1);diff(V,q2);diff(V,q3);diff(V,q4);diff(V,q5);diff(V,q6);diff(V,qd1);diff(V,qd2);diff(V,qd3);diff(V,qd4);diff(V,qd5);diff(V,qd6);zeros(numJ,numJ)];

matlabFunction(dV_dx,'file','grad/dV_dx','vars',[qd1,qd2,qd3,qd4,qd5,qd6,q2,q3,q4,q5,q6]);


dM_dx2 = vpa(diff(M,q2),5);
dM_dx3 = vpa(diff(M,q3),5);
dM_dx4 = vpa(diff(M,q4),5);
dM_dx5 = vpa(diff(M,q5),5);
dM_dx6 = vpa(diff(M,q6),5);
matlabFunction(dM_dx2,'file','grad/dM_dx2','vars',[q2,q3,q4,q5,q6]);
matlabFunction(dM_dx3,'file','grad/dM_dx3','vars',[q2,q3,q4,q5,q6]);
matlabFunction(dM_dx4,'file','grad/dM_dx4','vars',[q2,q3,q4,q5,q6]);
matlabFunction(dM_dx5,'file','grad/dM_dx5','vars',[q2,q3,q4,q5,q6]);
matlabFunction(dM_dx6,'file','grad/dM_dx6','vars',[q2,q3,q4,q5,q6]);        


%% generate the beta function for GRF

J_f = J(1:2,:).';
Mext = (J_f.'*J_f)\J_f.';
beta_grf = sigma*J_f;

dMext_dx1 = diff(Mext,q1);
dMext_dx2 = diff(Mext,q2);
dMext_dx3 = diff(Mext,q3);
dMext_dx4 = diff(Mext,q4);
dMext_dx5 = diff(Mext,q5);
dMext_dx6 = diff(Mext,q6);

matlabFunction(Mext,'file','grf/Mext','vars',{[q1,q2,q3,q4,q5,q6]});
matlabFunction(dMext_dx1,'file','grf/dMext_dx1','vars',{[q1,q2,q3,q4,q5,q6]});
matlabFunction(dMext_dx2,'file','grf/dMext_dx2','vars',{[q1,q2,q3,q4,q5,q6]});
matlabFunction(dMext_dx3,'file','grf/dMext_dx3','vars',{[q1,q2,q3,q4,q5,q6]});
matlabFunction(dMext_dx4,'file','grf/dMext_dx4','vars',{[q1,q2,q3,q4,q5,q6]});
matlabFunction(dMext_dx5,'file','grf/dMext_dx5','vars',{[q1,q2,q3,q4,q5,q6]});
matlabFunction(dMext_dx6,'file','grf/dMext_dx6','vars',{[q1,q2,q3,q4,q5,q6]});
% beta_grf = sigma*J_f/(J_f.'*J_f)*J_f.';
dbeta_dx1 = diff(beta_grf,q1);
dbeta_dx2 = diff(beta_grf,q2);
dbeta_dx3 = diff(beta_grf,q3);
dbeta_dx4 = diff(beta_grf,q4);
dbeta_dx5 = diff(beta_grf,q5);
dbeta_dx6 = diff(beta_grf,q6);

matlabFunction(beta_grf,'file','grf/beta_grf','vars',{[q1,q2,q3,q4,q5,q6],th});
matlabFunction(dbeta_dx1,'file','grf/dbeta_dx1','vars',{[q1,q2,q3,q4,q5,q6],th});
matlabFunction(dbeta_dx2,'file','grf/dbeta_dx2','vars',{[q1,q2,q3,q4,q5,q6],th});
matlabFunction(dbeta_dx3,'file','grf/dbeta_dx3','vars',{[q1,q2,q3,q4,q5,q6],th});
matlabFunction(dbeta_dx4,'file','grf/dbeta_dx4','vars',{[q1,q2,q3,q4,q5,q6],th});
matlabFunction(dbeta_dx5,'file','grf/dbeta_dx5','vars',{[q1,q2,q3,q4,q5,q6],th});
matlabFunction(dbeta_dx6,'file','grf/dbeta_dx6','vars',{[q1,q2,q3,q4,q5,q6],th});

%% generate friction for object function 
% we only consider the horizontal friction, so we use x_dot = J(1,:)*q_dot
x_vel = J(1,:)*[qd1;qd2;qd3;qd4;qd5;qd6];
fri = sigma_out(q1,q2,q3,q4,q5,q6,th)*(x_vel*x_vel);
fri_dx = [diff(fri,q1);
          diff(fri,q2);
          diff(fri,q3);
          diff(fri,q4);
          diff(fri,q5);
          diff(fri,q6);
          diff(fri,qd1);
          diff(fri,qd2);
          diff(fri,qd3);
          diff(fri,qd4);
          diff(fri,qd5);
          diff(fri,qd6);
          zeros(numJ,1)];
matlabFunction(fri,'file','obj/fri','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],th});
matlabFunction(fri_dx,'file','obj/fri_dx','vars',{[q1,q2,q3,q4,q5,q6,qd1,qd2,qd3,qd4,qd5,qd6],th});



