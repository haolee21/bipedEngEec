%% Generating the plotting function (link length, com pos.....)

l_feet = 0.3048;
l_calf = 0.45;
l_thigh =0.45;
l_torso = 0.9;
DH = [0, 0, 0, 0, 0, 0;...
      0, l_calf, 0, 0, 0, 0;...
      0, l_thigh, 0, 0, 0, 0;...
      0,  0, 0, 0, 0, 0;...
      0,l_thigh, 0, 0, 0, 0];
robot=DH2Robot(DH,1);


robot.links(1).r = [l_calf/2,0,0]';
robot.links(2).r = [l_thigh/2,0,0]';
robot.links(3).r = [l_torso/2,0,0]';
robot.links(4).r = [l_thigh/2,0,0]';
robot.links(5).r = [l_calf/2,0,0]';


jointP2 = turnRTtoMatrix(robot.A([1],[q1 q2 q3 q4 q5]))*[l_calf,0,0,1].';
jointP3 = turnRTtoMatrix(robot.A([1,2],[q1 q2 q3 q4 q5]))*[l_thigh,0,0,1].';
joint_head = turnRTtoMatrix(robot.A([1,2,3],[q1 q2 q3 q4 q5]))*[l_torso,0,0,1].';
jointP5 = turnRTtoMatrix(robot.A([1,2,3,4],[q1 q2 q3 q4 q5]))*[l_thigh,0,0,1].';
jointP6 = turnRTtoMatrix(robot.A([1,2,3,4,5],[q1 q2 q3 q4 q5]))*[l_calf,0,0,1].';

jointPos = [jointP2(1:3,:),jointP3(1:3,:),joint_head(1:3,:),jointP5(1:3,:),jointP6(1:3)];

com1 = turnRTtoMatrix(robot.A([1],[q1 q2 q3 q4 q5]))*[l_calf/2,0,0,1].';
com2 = turnRTtoMatrix(robot.A([1,2],[q1 q2 q3 q4 q5]))*[l_thigh/2,0,0,1].';
com3 = turnRTtoMatrix(robot.A([1,2,3],[q1 q2 q3 q4 q5]))*[l_torso/2,0,0,1].';
com4 = turnRTtoMatrix(robot.A([1,2,3,4],[q1 q2 q3 q4 q5]))*[l_thigh/2,0,0,1].';
com5 = turnRTtoMatrix(robot.A([1,2,3,4,5],[q1 q2 q3 q4 q5]))*[l_calf/2,0,0,1].';

comPos = [com1(1:3,:),com2(1:3,:),com3(1:3,:),com4(1:3,:),com5(1:3,:)];

matlabFunction([jointPos;comPos],'file','getRobotPos');