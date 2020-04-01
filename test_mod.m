
% This script will create a model to validate the possibility of parameters
% adaptation during motion
% Model will be 2-d R-R-R minipulator
clear;

syms l1 l2 l3 
syms q1 q2 q3 q1_dot q2_dot q3_dot q1_dot2 q2_dot2 q3_dot2
syms lc1 lc2 lc3 m1 m2 m3 g
syms I1zz I2zz I3zz I1xx I1yy I2xx I2yy I3xx I3yy
syms fx fy fz
totalVar = [l1 l2 l3 q1 q2 q3 q1_dot q2_dot q3_dot q1_dot2 q2_dot2 q3_dot2 lc1 lc2 lc3 m1 m2 m3 g I1zz I2zz I3zz fx fy fz];
var = [q1 q2 q3];
var_dot = [q1_dot q2_dot q3_dot];
var_dot2 = [q1_dot2 q2_dot2 q3_dot2];
DH = [0,0,0,0,0,0;...
    0,l1,0,0,0,0;...
    0,l2,0,0,0,0;...
    0,l3,0,0,0,0];
robot=DH2Robot(DH,0);
T0_1 = turnRTtoMatrix(robot.A([1],[q1,q2,q3,0]));
T1_2 = turnRTtoMatrix(robot.A([2],[q1,q2,q3,0]));
T2_3 = turnRTtoMatrix(robot.A([3],[q1,q2,q3,0]));
T3_E = turnRTtoMatrix(robot.A([4],[q1,q2,q3,0]));

p=[T0_1(1:3,4).';T1_2(1:3,4).';T2_3(1:3,4).';T3_E(1:3,4).'].';
rotMat = {T0_1(1:3,1:3),T1_2(1:3,1:3),T2_3(1:3,1:3),T3_E(1:3,1:3)};
pc = [lc1 0 0;lc2 0 0;lc3 0 0].';
I1=[0 0 0;0 0 0;0 0 I1zz];
I2=[0 0 0;0 0 0;0 0 I2zz];
I3=[0 0 0;0 0 0;0 0 I3zz];
m = [m1 m2 m3];
Ic ={I1 I2 I3};
tou = newtonEuler(rotMat,p,pc,m,Ic,var,var_dot,var_dot2,totalVar);
M=[diff(tou(1,1),q1_dot2),diff(tou(1,1),q2_dot2),diff(tou(1,1),q3_dot2);...
    diff(tou(1,2),q1_dot2),diff(tou(1,2),q2_dot2),diff(tou(1,2),q3_dot2);...
    diff(tou(1,3),q1_dot2),diff(tou(1,3),q2_dot2),diff(tou(1,3),q3_dot2)];
G=[diff(tou(1,1),g);diff(tou(1,2),g);diff(tou(1,3),g)];
M=simplify(M);
G=simplify(G);
