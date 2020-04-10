% This script will generate a 2-d biped model's walking trajectories

% the calculated torque profile will be used in openai gym box2d model
clear;
clc;
syms q1 q2 q3 q4 q5 q6 q7
syms q1_dot q2_dot q3_dot q4_dot q5_dot q6_dot q7_dot
syms q1_dot2 q2_dot2 q3_dot2 q4_dot2 q5_dot2 q6_dot2 q7_dot2


l_feet = 0.3048;
l_calf = 0.45;
l_thight =0.45;
l_torso = 0.9;
DH = [0, 0, 0, 0, 0, 0;...
      0, l_feet, 0, 0, 0, 0;...
      0, l_calf, 0, 0, 0, 0;...
      0, l_thight, 0, 0, 0, 0;
      0,l_thight, 0, 0, 0, 0;...
      0,l_calf, 0, 0, 0, 0;...
      0,l_feet, 0, 0, 0, 0];

% DH = [0, 0, 0, 0, 0, 0;...
%       0, l_feet, 0, 0, 0, 0;...
%       0, l_calf, 0, 0, 0, 0;...
%       0, l_thight, 0, 0, 0, 0;...
%       pi, 0, 0 ,0,0,0;...
%       0,l_thight, 0, 0, 0, 0];

robot=DH2Robot(DH,1);

% robot.plot([10,70,-15,35,180,15,-90,0]/180*pi);

totM = 70; %kg
m_feet = 0.0143*totM;
m_calf = 0.0475*totM;
m_thight = 0.105*totM;
m_torso = totM - m_feet-m_calf-m_thight;


robot.links(1).m = m_feet;
robot.links(2).m = m_calf;
robot.links(3).m = m_thight;
% robot.links(4).m = m_torso;
robot.links(4).m = m_thight;
robot.links(5).m = m_calf;
robot.links(6).m = m_feet;
robot.links(7).m =0;


robot.links(1).r = [l_feet/2,0,0]';
robot.links(2).r = [l_calf/2,0,0]';
robot.links(3).r = [l_thight/2,0,0]';
% robot.links(4).r = [l_torso/2,0,0]';
robot.links(4).r = [l_thight/2,0,0]';
robot.links(5).r = [l_calf/2,0,0]';
robot.links(6).r = [l_feet/2,0,0]';
robot.links(7).r =[0,0,0]';

I_feet = zeros(3);
I_feet(3,3) =m_feet*l_feet^2/12;
I_calf = zeros(3);
I_calf(3,3) = m_calf*l_calf^2/12;
I_thight = zeros(3);
I_thight(3,3) = m_thight*l_thight^2/12;
I_torso = zeros(3);
I_torso(3,3) = m_torso*l_torso^2/12;
robot.links(1).I = I_feet;
robot.links(2).I = I_calf;
robot.links(3).I = I_thight;
% robot.links(4).I = I_torso;
robot.links(4).I = I_thight;
robot.links(5).I = I_calf;
robot.links(6).I = I_feet;
robot.links(7).I = zeros(3);

% robot.links.r = [l_feet/2 0 0;...
%           l_calf/2 0 0;...
%           l_thight/2 0 0;...
%           l_torso/2 0 0;...
%           -l_thight/2 0 0;...
%           -l_calf/2 0 0 ;...
%           -l_feet/2 0 0 ].';


robot.links(1).Jm=0;
robot.links(2).Jm=0;
robot.links(3).Jm=0;
robot.links(4).Jm=0;
robot.links(5).Jm=0;
robot.links(6).Jm=0;
robot.links(7).Jm=0;
% robot.links(8).Jm=0;


syms q1 q2 q3 q4 q5 q6 q7
syms q1_dot q2_dot q3_dot q4_dot q5_dot q6_dot q7_dot
syms q1_dot2 q2_dot2 q3_dot2 q4_dot2 q5_dot2 q6_dot2 q7_dot2
mod.var = [q1 q2 q3 q4 q5 q6 q7];
mod.var_dot = [q1_dot q2_dot q3_dot q4_dot q5_dot q6_dot q7_dot];
mod.var_dot2 = [q1_dot2 q2_dot2 q3_dot2 q4_dot2 q5_dot2 q6_dot2 q7_dot2];
syms fx fy fz
mod.fext =[fx fy fz];

syms g % it has to be a symbolic variable since we need to use it to solve G
mod.g = g;
mod.total_var = [q1 q2 q3 q4 q5 q6 q7 q1_dot q2_dot q3_dot q4_dot q5_dot q6_dot q7_dot q1_dot2 q2_dot2 q3_dot2 q4_dot2 q5_dot2 q6_dot2 q7_dot2 fx fy fz g];
T0_1 = turnRTtoMatrix(robot.A([1],[q1,q2,q3,q4,q5,q6,q7,0]));
T1_2 = turnRTtoMatrix(robot.A([2],[q1,q2,q3,q4,q5,q6,q7,0]));
T2_3 = turnRTtoMatrix(robot.A([3],[q1,q2,q3,q4,q5,q6,q7,0]));
T3_4 = turnRTtoMatrix(robot.A([4],[q1,q2,q3,q4,q5,q6,q7,0]));
T4_5 = turnRTtoMatrix(robot.A([5],[q1,q2,q3,q4,q5,q6,q7,0]));
T5_6 = turnRTtoMatrix(robot.A([6],[q1,q2,q3,q4,q5,q6,q7,0]));
T6_7 = turnRTtoMatrix(robot.A([7],[q1,q2,q3,q4,q5,q6,q7,0]));
T7_E = turnRTtoMatrix(robot.A([8],[q1,q2,q3,q4,q5,q6,q7,0]));

p=[T0_1(1:3,4).';T1_2(1:3,4).';T2_3(1:3,4).';T3_4(1:3,4).';T4_5(1:3,4).';T5_6(1:3,4).';T6_7(1:3,4).';T7_E(1:3,4).'].';
rMat = {T0_1(1:3,1:3),T1_2(1:3,1:3),T2_3(1:3,1:3),T3_4(1:3,1:3),T4_5(1:3,1:3),T5_6(1:3,1:3),T6_7(1:3,1:3),T7_E(1:3,1:3)};

for i=1:length(mod.p)
    p(i) = simplify(p(i));
    rMat{i} = simplify(rMat{i});
end
mod.p = p;
mod.rMat = rMat;




I_feet = zeros(3);
I_feet(3,3) =m_feet*l_feet^2/12;
I_calf = zeros(3);
I_calf(3,3) = m_calf*l_calf^2/12;
I_thight = zeros(3);
I_thight(3,3) = m_thight*l_thight^2/12;
I_torso = zeros(3);
I_torso(3,3) = m_torso*l_torso^2/12;
mod.m = [m_feet,m_calf,m_thight,m_torso,m_thight,m_calf,m_feet];
mod.Ic = {I_feet,I_calf,I_thight,I_torso,I_thight,I_calf,I_feet};


      
dynEqn = newtonEuler(mod);    
% matlabFunction(dynEqn.tau,'File','Tau');
matlabFunction(vpa(subs(dynEqn.M,g,9.81),5),'File','M_mat');
matlabFunction(vpa(subs(dynEqn.G,g,9.81),5),'File','G_mat');
matlabFunction(vpa(subs(dynEqn.F,g,9.81),5),'File','F_mat');
% matlabFunction(dynEqn.V,'File','V_mat','Optimize',false);
